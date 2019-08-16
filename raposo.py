'''
Implements the Raposo line simplification
algorithm, in which consecutive vertices in hexagonal 
grid cells are collapsed to a single
point. For details on the algorithm, see

Raposo (2013) Scale-specific automated line simplification by
vertex clustering on a hexagonal tessellation. CAGIS 40(5):427-443

Note: The 'untwisting' algorithm to remove self-intersections
described in the above paper has not been implemented

Boring license stuff:
=========================================
The MIT License (MIT)
Copyright (c) 2019 Barry Kronenfeld
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
================================
'''

class cHexGrid():
    '''Represents a hexagonal grid with the following properties:
        (1) the bottom left grid cell is centered on the bottom
             left corner of the bounding box,
        (2) hexagons in a column have the same x-value,
        (3) hexagons in odd-numbered columns are a half-step above
             hexagons in even-numbered columns in same row
        '''
    __sqrt_3 = 3.0**0.5
    def __init__(self, top, left, bottom, right, radius):
        '''Defines a hexagonal grid covering the bounding box defined by top, left, bottom and right'''
        self.__top=top
        self.__left=left
        self.__bottom=bottom
        self.__right=right
        self.__radius=radius
        self.__colWid = 1.5 * self.__radius
        self.__rowHt = self.__radius * self.__sqrt_3
    def hexRowCol(self,x,y):
        '''Returns the row & column of the hexagon containing the input point.'''
        # determine nearest column
        C = int((x+self.__colWid/2)/self.__colWid)
        # determine if point is left or right of column center
        nearX, nearY = self.hexCenter(0,C)
        if x < nearX:
            cOff = -1
        else:
            cOff = 1
        # determine nearest row
        # if nearest column is even, candidates will include row below
        # otherwise candidates will include row above
        if C % 2 == 0:
            R = int((y+self.__rowHt/2)/self.__rowHt)
            rOff = -1
        else:
            R = int(y/self.__rowHt)
            rOff = 1
        # make list of three candidate hexagons
        candidates=[[R,C],[R,C+cOff], [R+rOff,C+cOff]]
        # winner is candidate closest to input point
        mind2=self.__radius*self.__radius*2
        candidate=candidates[0]
        cX, cY = self.hexCenter(candidate[0],candidate[1])
        d2 = (x-cX)**2 + (y-cY)**2
        winner=candidate
        for candidate in candidates[1:]:
            cX, cY = self.hexCenter(candidate[0],candidate[1])
            d2 = (x-cX)**2 + (y-cY)**2
            if d2 <= mind2:
                mind2 = d2
                winner=candidate
        return [winner[0], winner[1]]

    def hexCenter(self,row,col):
        '''Returns the center of the hexagon at the input row & column.'''
        x = col*self.__colWid
        y=row*self.__rowHt
        if col%2==1: # odd column
            y+=self.__rowHt/2.0
        return x,y

def quant_avg(points,startid,nextstartid):
    '''Returns the mean center of all vertices in group.'''
    n=nextstartid-startid
    meanx=sum([pt[0] for pt in points[startid:nextstartid]])/n
    meany=sum([pt[1] for pt in points[startid:nextstartid]])/n
    return ((meanx,meany))

def quant_mid(points, startid, nextstartid):
    '''Returns the midpoint between the start and end vertices.'''
    p1=points[startid]
    p2=points[nextstartid-1]
    return ((p1[0]+p2[0])/2.0,(p1[1]+p2[1])/2.0)

def simplify_raposo(points,hexradius,quant_method=quant_avg, keep_first_last=True):
    '''Constructs a simplified line using the hexagon clustering method described in:

        PARAMETERS
        points: list/set of coordinate pairs
        ............e.g. [(16.3, 47.1), (-256.5, 132.5), (-780.5, 191.4), (-1791.7, 213.2)]
        hexradius: length of the radius/sides of hexagons in grid
        ............to remove more vertices, use larger radius
        quant_method: quantization function on consecutive vertices in the same hexagon
        ............built-in options: quant_avg, quant_mid
        keep_first_last: if true, first and last points from original line will be retained.'''

    # determine bounds of line and create hexagon grid
    x=[pt[0] for pt in points]
    y=[pt[1] for pt in points]
    hexgrid=cHexGrid(max(y),min(x),min(y),max(x),hexradius)
    # determine hexagon of each vertex
    hexRC = [hexgrid.hexRowCol(pt[0],pt[1]) for pt in points]
    # create simplified line and add first point (if keep_first_last)
    simplified_line=[]
    if keep_first_last and hexRC[0]==hexRC[1]:
        simplified_line.append(points[0])
    # identify consecutive points in each hexagon, quantize (aggregate)
    # and add to simplified line
    startid=0
    for i in range(len(points)):
        if hexRC[i] != hexRC[startid]:
            simplified_line.append(quant_method(points,startid,i))
            startid=i
    simplified_line.append(quant_method(points,startid,len(points)))
    # add last point (if keep_first_last)
    if keep_first_last and hexRC[-2]==hexRC[-1]:
        simplified_line.append(points[-1])
    # return result
    return simplified_line

def simplify_to_n_points_raposo(n,points,max_tries=25,quant_method=quant_avg, keep_first_last=True):
    '''Attempts to simplify the input line to exactly n vertices using
       the Raposo line simplification algorithm, by searching for a
       hexagon radius that will produce the desired number of points.

       PARAMETERS:
        n: target number of points in simplified line
        max_tries: maximum number of iterations to try before giving up
        other parameters are the same as function simplify_raposo'''

    # determine bounds of line
    x=[pt[0] for pt in points]
    y=[pt[1] for pt in points]
    dx = max(x) - min(x)
    dy = max(y) - min(y)
    # default low radius to zero, which would result in infinite vertices
    low_radius=0
    # default high radius to diagonal of bounding box, which would result in a single vertex
    high_radius=((dx**2)+(dy**2))**0.5
    # pivot radius is average of low and high radius
    pivot = (low_radius + high_radius)/2.0
    # try simplifying line with pivot radius and calc difference between num vertices & target
    simplified = simplify_raposo(points,pivot,quant_method,keep_first_last)
    dn=abs(len(simplified)-n) # difference between target and current n
    numtries=0
    # loop until we get the right number of points or make enough tries
    while dn > 0 and numtries <= max_tries:
        # move pivot
        if len(simplified) < n: # too few points, need to decrease radius
            high_radius=pivot
        else: # too many points, need to increase radius
            low_radius=pivot
        pivot = (low_radius + high_radius)/2.0
        # simplify again
        simplified = simplify_raposo(points,pivot,quant_method,keep_first_last)
        # see if we got closer
        new_dn = abs(len(simplified)-n)
        if new_dn < dn: # if so, reset try count
            numtries=0
        else: # otherwise, increase try count
            numtries +=1
        # reset dn
        dn=new_dn
    return simplified

if __name__ == '__main__':
    pass
