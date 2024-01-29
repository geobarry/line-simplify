#-------------------------------------------------------------------------------
# Name:        geom utils
# Purpose:     various geometry functions
# License:     MIT License
#-------------------------------------------------------------------------------

from .data_utils import min_id,sequence_info
import math as __m
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from decimal import Decimal as dec


def clockwise_sequence(edges):
    """
    Determines a clockwise sequence of the input edges around their shared 1st vertex.
    * Can handle edges that start in the same direction but then later split, but 
    no accomodations made for floating point errors or precision tolerance.

    Parameters
    ----------
    edges : list of lists of (x,y) tuples
        List of edges. Each edge must have the same first vertex. No edge should be 
        identical to or a subset of another edge

    Returns
    -------
    seq : list
        Sequence of edge indices, so that result[i] represents the edge in position i.

    """
    # HELPER FUNCTION
    def get_subedges(tied_edges):
        """Creates a set of sub edges from split vertex, including reverse leader segment
        Input edges should start from same pt and go in same direction"""
        origin = tied_edges[0][0] # get shared 1st vertex
        v2list = [e[1] for e in tied_edges] # get 2nd vertex of all edges
        v2set = set(v2list) # find unique pts
        # get closest to origin
        if len(v2set) != 1:
            dist_list = [distance_pts(origin,v2) for v2 in v2list]
            min_id = min_id(dist_list)
#            min_id, min_d = min(enumerate(dist_list), key=itemgetter(1))
            fulcrum = v2list[min_id]
        else:
            fulcrum = v2list[0]
        # create sub_edges by truncating
        def check_subedge(e,fulcrum):
            if e[0] == fulcrum:
                return e
            else:
                return [fulcrum] + e
        subedges=[check_subedge(e[1:], fulcrum) for e in tied_edges]
        subedges.append([fulcrum,origin]) # add reverse leader
        return subedges
    
    # MAIN CODE
    # error checking

    # *** need to handle case where "edge" has only one vertex
    
    
    # *** by simply placing it in any place in sequence
    
    
    # pull out degenerate edges
    degenerate_edge_ids = [i for i in range(len(edges)) if len(edges[i])<2]
    legitimate_edge_ids = [i for i in range(len(edges)) if len(edges[i])>=2]
    legitimate_edges = [edges[i] for i in legitimate_edge_ids]
    # work on legitimate edges
    n = len(legitimate_edge_ids)
    bearings = [-1*bearing(edge[0],edge[1]) for edge in legitimate_edges]   
    seq,ranks = sequence_info(bearings)
    r=0
    # obtain and handle sets of ties
    while r < n-1:
        # get list of identical bearings
        tie_ranks = [r]
        while r < n-1 and bearings[seq[r+1]] == bearings[seq[r]]:
            r += 1
            tie_ranks.append(r)
        # handle set of identical bearings
        if len(tie_ranks) > 1:
            # get subedges
            tied_edges = [legitimate_edges[seq[r]] for r in tie_ranks]
            subedges = get_subedges(tied_edges)
            # get subedge order
            subseq = clockwise_sequence(subedges)
            lead_pos = subseq.index(len(subedges)-1)
            subseq = subseq[lead_pos+1:] + subseq[:lead_pos]
            # make replacements
            replaces = [(tie_ranks[i],seq[tie_ranks[subseq[i]]]) for i in range(len(tie_ranks))]
            for repl in replaces:
                seq[repl[0]] = repl[1]
        # move to next rank
        r += 1    
    # convert to original index
    seq=[legitimate_edge_ids[i] for i in seq]
    # add back in degenerate edge ids
    seq = seq + degenerate_edge_ids
    return seq 

def subtract_vectors(a,b):
    """
    Computes the vector difference a-b

    Parameters
    ----------
    a : (x,y) tuple
    b : (x,y) tuple

    Returns
    -------
    (x,y) tuple
    """
    return (a[0]-b[0],a[1]-b[1])

def dot_product(a,b):
    """
    Computes the dot-product of two vectors.

    Parameters
    ----------
    a : (x,y) tuple
    b : (x,y) tuple

    Returns
    -------
    float
    """
    return a[0]*b[0] + a[1]*b[1]

def length(v):
    """
    Computes the length of a vector

    Parameters
    ----------
    v : (x,y) tuple

    Returns
    -------
    float

    """
    return __m.sqrt((v[0]**2)+(v[1]**2))

def bearing(p1,p2):
    """
    Computes the direction from p1 to p2

    Parameters
    ----------
    p1 : (x,y) tuple
    p2 : (x,y) tuple

    Returns
    -------
    float
        The bearing in radians counter-clockwise from horizontal-right (=0).

    """
    v = subtract_vectors(p2, p1)
    return __m.atan2(v[1], v[0])

def angle(fulcrum,a,b):
    """
    Computes the angle from a to b, measured in radians counterclockwise.

    Parameters
    ----------
    fulcrum :: (x,y) tuple
        fulcrum of angle
    a :: (x,y) tuple
        end of one segment from fulcrum 
    b :: (x,y) tuple
        end of other segment from fulcrum

    Returns
    -------
    float in range [0,pi]

    """
    # uses method in https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
    # which is measured to be numerically stable to ~10^-15
    # get distances of triangle formed
    c=distance_pts(a,b)
    a=distance_pts(fulcrum,a)
    b=distance_pts(fulcrum,b)
    if b<0 or c<0: # this should never happen
        return None
    if b >= c:
        mu = c-(a-b)
    else:
        mu = b-(a-c)
    numerator=((a-b)+c)*mu
    denominator=(a+(b+c))*((a-c)+b)
    half_tangent = __m.sqrt(numerator/denominator)
    theta = 2 * __m.atan(half_tangent)    
    return theta 

def angle_quick(fulcrum,a,b):
    # approx. 18% faster, but errors up to 10 x e-8 occur approx. 1/million times
    a=subtract_vectors(a,fulcrum)
    b=subtract_vectors(b, fulcrum)
    dotprod = dot_product(a, b)
    a = length(a)
    b = length(b)
    cos_theta = dotprod/(a*b)
    return __m.acos(cos_theta)

def area(pts,absolute=False):
    """
    Computes the clockwise area of the polygon defined by the points.
    Args:
        pts: list of (x,y) tuples. Assumes last and first vertex are different.
             Cannot handle if input points contain None coordinates.
        absolute: if true, returns the absolute value
    Returns:
        float representing the area
    """
    if pts[len(pts)-1] != pts[0]:
        pts.append(pts[0])
    a=[(pts[i+1][0]-pts[i][0])*(pts[i][1]+pts[i+1][1]) for i in range(len(pts)-1)]
    A=sum(a)/2
    if absolute:
        return abs(A)
    else:
        return A

def triangle_area(a,b,c, robustZero = False):
    '''returns the area of the triangle 
       If robustZero: returns zero if 
       zero is calculated in any one of three ways, or 
       if both negative and positive values are 
       calculated by the different methods'''
    # reduce numbers to avoid floating point precision
    minx = min(a[0],b[0],c[0])
    miny = min(a[1],b[1],c[1])
    a,b,c=(a[0]-minx,a[1]-miny),(b[0]-minx,b[1]-miny),(c[0]-minx,c[1]-miny)
    # calculate area three different ways
    a1=(b[0]-a[0])*(b[1]+a[1])
    a2=(c[0]-b[0])*(c[1]+b[1])
    a3=(a[0]-c[0])*(c[1]+a[1])
    areas = [a1+a2+a3,a2+a3+a1,a3+a1+a2]
    if robustZero:
        # if any are zero, return zero
        if min([abs(x) for x in areas])==0:
            return 0
        # if there is disagreement over sign, return zero
        if (areas[0] > 0) != (areas[1] > 0):
            return 0
        if (areas[0] > 0) != (areas[2] > 0):
            return 0
    # return area with minimum absolute value
    w=0
    if abs(areas[1])<abs(areas[0]):
        w=1
    if abs(areas[2])<abs(areas[w]):
        w=2
    area = areas[w]
    return area/2

def distance_pts(a,b):
    ''' computes the distance between two points.
        inputs should be lists of two coordinates.'''
    return  __m.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)

def cycle_next(id, cycle_length):
    """
    computes the next id in a cycle

    Parameters
    ----------
    id : integer
        The index of an item in the cycle.
    cycle_length : integer
        The length of the cycle.

    Returns
    -------
    integer
        The index of the next item in the cycle

    """
    if id==cycle_length-1:
        return 0
    else:
        return id+1

def cycle_prev(id,cycle_length):
    """
    computes the previous id in a cycle

    Parameters
    ----------
    id : integer
        The index of an item in the cycle.
    cycle_length : integer
        The length of the cycle.

    Returns
    -------
    integer
        The index of the previous item in the cycle
        """
    if id==0:
        return cycle_length-1
    else:
        return id-1

def nextID_in_poly(poly, id):
    """
    Determines the id of the next point on a standard polygon.

    Parameters
    ----------
    poly :: List of (x,y) tuples
        Representation of a polygon, with identical first and last vertices.
    id : Integer

    Returns
    -------
    Integer

    """
    if id==len(poly)-1:
        return 1
    elif id== len(poly)-2:
        return 0
    else:
        return id+1

def prevID_in_poly(poly, id):
    """
    Determines the id of the previous point on a standard polygon.

    Parameters
    ----------
    poly :: List of (x,y) tuples
        Representation of a polygon, with identical first and last vertices.
    id : Integer

    Returns
    -------
    Integer

    """
    if id==0:
        return len(poly)-2
    else:
        return id-1


def clockwise(A,B,C):
    """Determines if points A,B & C are sequenced clockwise around a triangle.
    Args:
        A,B,C: (x,y) tuples
    Returns:
        True if points are sequenced clockwise. Result is unpredictable if points
        are collinear."""
    return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])

def intersect(A,B,C,D):
    """Quickly determines if segment AB intersects segment CD.
    Args:
        A,B,C,D: (x,y) tuples
    Returns:
        True if segments intersect. Result is unpredictable if
        lines are parallel or segments meet at endpoin.
        """
    return clockwise(A,C,D) != clockwise(B,C,D) and clockwise(A,B,C) != clockwise(A,B,D)

def intersection(A,B,C,D,infinite=True):
    """ Returns the intersection point of two lines AB & CD.
        A,B,C,D and return value are all lists of two coordinates each.
        If lines are parallel or do not intersect, returns a pair of Nones.
        Code taken from Stephen Wise, pp. 48-9
        Args:
            A,B,C,D: (x,y) tuples
            infinite: If False, will return (None, None)
        Returns:
            (x,y) tuple of coordinates of intersection point"""

##    if ultra_precise:
##        A=(Decimal(A[0]),Decimal(A[1])) # dec_pt(A)
##        B=(Decimal(B[0]),Decimal(B[1])) # dec_pt(B)
##        C=(Decimal(C[0]),Decimal(C[1])) # dec_pt(C)
##        D=(Decimal(D[0]),Decimal(D[1])) # dec_pt(D)
    if A[0]==B[0]:
        if C[0]==D[0]:
            xp,yp=None,None # lines are parallel
        else: # first line vertical
            b2=(D[1]-C[1])/(D[0]-C[0])
            a2=C[1]-b2*C[0]
            xp = A[0]
            yp = a2+b2*xp
    else:
        if C[0]==D[0]: # second line vertical
            if B[1]==None or A[1]==None or B[0]==None or A[0]==None:
                hell="yeah!"
            b1=(B[1]-A[1])/(B[0]-A[0])
            a1=A[1]-b1*A[0]
            xp = C[0]
            yp = a1+b1*xp
        else: # neither line vertical
            b1=(B[1]-A[1])/(B[0]-A[0])
            b2=(D[1]-C[1])/(D[0]-C[0])
            a1=A[1]-b1*A[0]
            a2=C[1]-b2*C[0]
            if b1==b2:
                xp,yp = None,None # lines are parallel
            else:
                xp = -(a1-a2)/(b1-b2)
                yp = a1+b1*xp
    # test whether intersection point falls on either line
    if infinite == False and xp != None:
        if (A[0]-xp)*(xp-B[0]) < 0 or (C[0]-xp)*(xp-D[0]) < 0 or (A[1]-yp)*(yp-B[1]) < 0 or (C[1]-yp)*(yp-D[1]) < 0:
            xp,yp = None, None
    if xp != None:
        xp=float(xp)
    if yp != None:
        yp = float(yp)
    return (xp,yp)

def distance_pt_line(p,a,b):
    """Computes the distance from point p to the infinite line through ab"""
    trianglearea=abs(area([a,b,p]))
    line_length=distance_pts(a,b)
    if line_length==0:
        return distance_pts(p,a)
    else:
        return 2*trianglearea/line_length

def X_Right_of_AB(A,B,X,tolerance=0):
    """
    Returns True if the input point or the first point on the input line
    is on or right of the line from A to B.

    A: first point on baseline (tuple of two floats)
    B: second point on baseline
    X: either a point (tuple of two floats) or a list of points
    """
    # get first point from X if it is a line
    if isinstance(X[0],list):
        X=X[0]
    # determine area of ABC
    ax,ay,bx,by,xx,xy=A[0],A[1],B[0],B[1],X[0],X[1]
    minx,miny=min(ax,bx,xx),min(ay,by,xy)
    ax,bx,xx=ax-minx,bx-minx,xx-minx
    ay,by,xy=ay-miny,by-miny,xy-miny
    area=ax*xy+bx*ay+xx*by-ax*by-bx*xy-xx*ay
    # check precision - THIS MIGHT REALLY SLOW THINGS DOWN!
    #tolerance=0
    if tolerance > 0:
        d = distance_pt_line(X,A,B)
        if d < tolerance:
            A = (dec(A[0]),dec(A[1]))
            B = (dec(B[0]),dec(B[1]))
            X = (dec(X[0]),dec(X[1]))
            ax,ay,bx,by,xx,xy=A[0],A[1],B[0],B[1],X[0],X[1]
            minx,miny=min(ax,bx,xx),min(ay,by,xy)
            ax,bx,xx=ax-minx,bx-minx,xx-minx
            ay,by,xy=ay-miny,by-miny,xy-miny
            area=ax*xy+bx*ay+xx*by-ax*by-bx*xy-xx*ay
        
    # if area is positive, C is right of AB
    if area >= 0:
        return True
    else:
        return False

def X_Right_of_ABC(A,B,C,X, tolerance=0):
    """ Returns TRUE if X is right of combined line ABC
        Also returns TRUE if X is exactly on combined line ABC
    """
    # Check point C is on which side of line AB
    if X_Right_of_AB(A,B,C,tolerance) == False:   # C is left of AB so ABC turns to left
        # If X is right of (or on) either AB or BC, X is right of combined line ABC
        if X_Right_of_AB(A,B,X,tolerance) == True or X_Right_of_AB(B,C,X,tolerance) == True:
            return True
        else:
            return False
    else: # C is right of AB so ABC turns to right
        # If X is left of either AB or BC, X is left of combined line ABC
        if X_Right_of_AB(A,B,X,tolerance) == False or X_Right_of_AB(B,C,X,tolerance) == False:
            return False
        else:
            return True

def __isLeft(p0,p1,p2):
    """Returns a float indicating if p2 is left or right of line from p0 and p1
        p0,p1,p2: tuples representing point coordinates
        return: float > 0 if p2 is left
                      < 0 if p2 is right
                      = 0 if p2 is on line"""
    p0x = p0[0]
    p0y = p0[1]
    p1x = p1[0]
    p1y = p1[1]
    p2x = p2[0]
    p2y = p2[1]
    return ((p1x-p0x)*(p2y-p0y)-(p2x-p0x)*(p1y-p0y))

def winding_number(pt,polygon):
    """ Calculates the winding number of a point in polygon
        
    Args: 
    -----
            pt: (x,y) tuple 
            
            polygon: list of all the points in this polygon
    Returns:
    --------

    the winding number, equal to 0 when pt is outside polygon, and equal to 
    winding number just to the right if pt is on the polygon boundary (?)
    """
    # This works by counting the number of edges
    # that cross a ray extending right from pt
    # upwards: winding number decreases
    # downwards: winding number increases
    
    # clean polygon - make sure last point is same as first
    if polygon[-1] != polygon[0]:
        polygon.append(polygon[0])
    # Initialize winding number
    w = 0
    # Loop through all edges of the polygon 
    for i in range(len(polygon)-1):            
        # Does poly edge cross horizontal line through point?
        # An upward edge includes its starting endpoint and excludes its final endpoint
        if polygon[i][1] <= pt[1]:
            if polygon[i+1][1] > pt[1]:
                # Is pt strictly left of polygon edge?
                if __isLeft(polygon[i],polygon[i+1],pt) > 0:
                    w -= 1
        else:
            # A downward edge excludes its starting endpoint and includes its final endpoint
            try:
                if polygon[i+1][1] <= pt[1]:
                    # Check pt is strictly right of polygon edge
                    if __isLeft(polygon[i],polygon[i+1],pt) < 0:
                        w += 1
            except:
#                print("****************************")
#                print(polygon_pts)
#                print("i: {}".format(i))
#                print("pt count: {}".format(len(polygon_pts)))
#                print("polygon_pts: {}".format(type(polygon_pts)))
#                print("polygon_pts[i+1]: {}".format(type(polygon_pts[i+1])))
#                print("polygon_pts[i+1][0]: {}".format(type(polygon_pts[i+1][0])))
#                print("****************************")
                pass
    return w 

def segmentOnLine(seg,line):
    """Determines if the input line segment is a proper part of the input line,
    i.e. if it consists of two consecutive segments on the input line.
    Args:
        seg: list two points comprising the segment to check
        line: list of any number of points comprising the full line
    Returns:
        True if the two segment points are consecutive points on the input line,
        False otherwise."""
    for i in range(len(line)-1):
        start=line[i]
        # check forward direction
        if start[0]==seg[0][0] and start[1]==seg[0][1]:
            end=line[i+1]
            if end[0]==seg[1][0] and end[1]==seg[1][1]:
                return True
        # check backwards direction
        if start[0]==seg[1][0] and start[1]==seg[1][1]:
            end=line[i+1]
            if end[0]==seg[0][0] and end[1]==seg[0][1]:
                return True
    return False

def same_point(p1,p2,tolerance):
    '''Determines if p1 and p2 are the same point, i.e. x & y coordinates
       are both within tolerance of each other'''
    if abs(p1[0]-p2[0]) > tolerance:
        return False
    if abs(p1[1]-p2[1]) > tolerance:
        return False
    return True

def pt_on_seg(p,s1,s2,countEndPoint=True,tolerance=0.000001,robustZero = False):
    """
    Returns TRUE if p is between S1 and S2, inclusive
    Returns FALSE if p is not exactly on line segment (no tolerance)
    args:
        p: (x,y) tuple representing point
        s1: (x,y) tuple representing one end of line segment
        s2: (x,y) tuple representing other end of line segment
        countEndPoint: value returned if p is identical to s1 or s2 within tolerance
    """
    # check for exact matches first
    if same_point(p,s1,tolerance):
        return countEndPoint
    if same_point(p,s2,tolerance):
        return countEndPoint
#    # p is on line through s1 and s2 only if area is 0
#    if triangle_area(p,s1,s2,robustZero) != 0:
#        return False
    # check distance from p to line through s1 & s2
    if distance_pt_line(p,s1,s2) > tolerance:
        return False
    # p is on line only if p.x is between s1.x and s2.x
    # and same is true for Y coordinates
    return ptsMonotonic(s1,p,s2)

def pt_on_ray(pt, ray_start, ray_pt, tolerance=0.000001):
    """
    Returns TRUE if pt is on ray (within tolerance perpendicularly)

    args:
        pt: (x,y) tuple representing point
        ray_start: (x,y) tuple representing origin of ray
        ray_pt: (x,y) tuple representing another point on ray
    """  
    rs=ray_start
    rp=ray_pt
    # check easy case
    if pt == ray_pt:
        return True
    # check if pt is on correct side of ray origin
    if (pt[0] >= rs[0]) != (rp[0] >= rs[0]):
        return False
    if (pt[1] >= pt[1]) != (rp[1] >= rs[1]):
        return False
    # check distance from pt to ray
    if distance_pt_line(pt,rs,rp) > tolerance:
        return False
    else:
        return True

def isMonotonic(a, b, c, tolerance=0):
    """Returns True if {a,b,c} is a monotonic sequence """
    if (a <= b+tolerance) == (b <= c+tolerance):
        return True
    elif (a >= b-tolerance) == (b >= c-tolerance):
        return True
    else:
        return False

def ptsMonotonic(a,b,c):
    """returns True of both x- and y- coordinates of points a,b,c are monotonic sequences"""
    if not isMonotonic(a[0],b[0],c[0]):
        return False
    if not isMonotonic(a[1],b[1],c[1]):
        return False
    return True

def removeConsecutiveCoincident(pts, track_id=0, tol=0.000001):
    '''
    Removes consecutive points that are within tolerance in x or y direction.
    Args:
        pts: list of (x,y) tuples
        tol: distance below which points will be considered identical
        track_id: id of point to track
    Returns:
        (list of (x,y) tuples, integer id of pt tracked)
    '''
    # include first point
    new_pts = [pts[0]]
    track_id2=0
    for i in range(1,len(pts)):
        # include point if different from previous 
        if (pts[i-1] != pts[i]) and distance_pts(pts[i-1] , pts[i]) > tol:
            new_pts.append(pts[i])
        # update position of point to track
        if i==track_id:
            track_id2=len(new_pts)-1
    return (new_pts,track_id2)
    
def _whisker_tips(p):
    """Identifies the vertices that are tips of whiskers.
    Args:
        p  a list of xy tuples, 
           with first and last point identical.
    Returns:
        a list of IDs of whisker tip vertices        
   """
    r = []
    # we'll refer to each vertex as b, skipping the first point
    for b in range(1, len(p)):
        # and the previous and next vertices as a & c
        a = b - 1
        if b==len(p)-1:
            c = 1
        else:
            c = b + 1
        A,B,C=p[a],p[b],p[c]
        # look for whiskers
        if pt_on_seg(A,B,C, countEndPoint = True, robustZero = True) or pt_on_seg(C,B,A, countEndPoint = True, robustZero = True):
            r.append(b)
    return r

def _one_per_sequence(ids,n):
    """Returns a new list of ids, including
       one id in each sequential group of input ids to a list with n elements.
       Assumes ids are in ascending order, with 0 excluded, and looks for
       sequential groups that 'wrap around' """
    # get the first id that doesn't have something in front of it
    pos = 0
    if ids[len(ids)-1]==n-1:   # last id is end of chain
        dontuse = 1 # can't use point 1
        while pos < len(ids) and ids[pos]==dontuse:
            pos+=1
            dontuse +=1
    # go through ids
    r=[]
    while pos < len(ids):
        r.append(ids[pos])
        while pos+1 < len(ids) and ids[pos+1] == ids[pos]+1:
            pos+=1
        pos +=1
    return r

def _remove_whiskers(polypts,whiskers,track_id):
    # attempts to remove each whisker from the polygon
    # but if there are multiple whisker points in a row,
    # only removes the last one
    # assumes whiskers are in ascending order
    # and that polygon first and last points are identical
    # and if first/last point is a whisker, it is included as last index

    # get points to remove
    # don't remove a "whisker" if it immediately follows another whisker
    whiskers = _one_per_sequence(whiskers,len(polypts))
    # if last point is a whisker, update xy of first point
    if whiskers[len(whiskers)-1]==len(polypts)-1:
        polypts[0] = polypts[len(polypts)-2]
    # create new polygon by replacing points in sequence, avoiding whiskers
    p=polypts[:-1*len(whiskers)]
    w=0 # whiskers encountered so far
    if len(whiskers)<w+1:
        next_exclude = len(p)
    else:
        next_exclude = whiskers[w]
    # go through all points
    for i in range(len(p)):
        # determine if next point needs to be excluded
        if i+w == next_exclude:
            # if so, increment number of exclusions so far
            w += 1
            if len(whiskers)<w+1:
                next_exclude = len(polypts)
            else:
                next_exclude = whiskers[w]
            # reduce id of track_pt if track_pt occurs after excluded pt
            if track_id >= i+1:
                track_id -=1
        # set i-th point to i-th non-excluded point
        p[i] = polypts[i+w]
    p,track_id=removeConsecutiveCoincident(p,track_id)
    return p,track_id

def remove_whiskers(polypts, track_id, verbose=False):
    """Removes all whiskers from the input polygon. A whisker is a section where
       the polygon boundary goes to a point and then backtracks from that point.
       Args:
            polypts: a list of xy tuples representing the original polygon
       Returns:
            a list of xy tuples representing the same polygon sans whiskers
    """
    w = _whisker_tips(polypts)
    if verbose:
        plot(polypts,"init")
    iteration=0
    while len(w) > 0:
        if verbose:
            print("Whiskers")
            print(w)
            for id in w:
                print(polypts[id])
        polypts,track_id = _remove_whiskers(polypts,w,track_id)
        if verbose:
            print("Polygon after whisker removal")
            iteration += 1
            plot(polypts,"iteration {}".format(iteration))
            for i in range(len(polypts)):
                print("{}\t{}\t{}".format(i,polypts[i][0],polypts[i][1]))
        w = _whisker_tips(polypts)
    return (polypts,track_id)

def coincident_point_list(pts, tolerance):
    """
    Records indices of pts coincident with each pt.
    Args:
        pts: list of (x,y) tuples. If your first & last pts are identical and you 
        don't want duplicates, it's YOUR job to remove the last pt
    Returns:
        list of lists, where each inner list contains indices of coincident pts
    """
    r=[]
    kt = cKDTree(pts)
    for i in range(len(pts)):
        pt=pts[i]
        coincident=kt.query_ball_point(pt,tolerance)
        coincident.remove(i)
        r.append(coincident)
    return r

def possiblyIntersectingSegmentPairs(poly, tolerance=0.000001, report_complexity=False):
    # polygon will have n vertices, n-1 segments
    r=[]
    n=len(poly)-1
    def next(id):
        if id==n-1:
            return 0
        else:
            return id+1
    def minx(id):
        return min(poly[id][0],poly[next(id)][0])
    def maxx(id):
        return max(poly[id][0],poly[next(id)][0])
    in_range=set() # set of segments with bounds around current x
    # get segments in ascending order of min x coordinate
    sorted_ids=sorted(range(len(poly)-1),key=lambda id: minx(id))
    # loop
    operations=0
    for id in sorted_ids:
        operations += len(in_range)
        # get current x
        curx=minx(id)
        # remove segments from range
        out_of_range=[id for id in in_range if maxx(id)+tolerance < curx]
        for partner in out_of_range:
            in_range.remove(partner)
        # pair current segment with segments in range
        for partner in in_range:
            if partner != id+1 and partner != id-1:
                pair=(min(id,partner),max(id,partner))
                if pair != (0,len(poly)-1):
                    
                    r.append((min(id,partner),max(id,partner)))
        # add current segment to range
        in_range.add(id)
    if report_complexity:
        return (r,operations)
    else:
        return r

def boundingBoxesOverlap(geom1, geom2, tolerance=0.000001, verbose=False):
    '''Determines if the axis-aligned boxes around two
       input geometries overlap within the given tolerance.
       Args:
           geom1: list of (x,y) tuples
           geom2: list of (x,y) tuples
       Returns:
           boolean (true or false)'''
    # calculate mins & maxes
    x1=[p[0] for p in geom1]
    x2=[p[0] for p in geom2]
    y1=[p[1] for p in geom1]
    y2=[p[1] for p in geom2]
    minx1=min(x1)
    maxx1=max(x1)
    minx2=min(x2)
    maxx2=max(x2)
    miny1=min(y1)
    maxy1=max(y1)
    miny2=min(y2)
    maxy2=max(y2)
    # check x
    if minx2 > maxx1 + tolerance:
        return False
    if minx1 > maxx2 + tolerance:
        return False
    # check y
    if miny2 > maxy1 + tolerance:
        return False
    if miny1 > maxy2 + tolerance:
        return False
    # checks out
    return True
    
def plot(poly, label, tolerance=0.000001, pt_vals=None):
    """
    Plots the polygon, highlighting intersections.

    Parameters
    ----------
    pts :: list of (x,y) tuples
        The points to be plotted.
    label :: string
        The plot title.
    tolerance :: float 
        Points within this distance will be marked as coincident
    pt_vals :: values
        Will be added to point labels.

    Returns
    -------
    None.

    """
    
    oneCount=[]
    twoCount=[]
    threeCount=[]
    fourCount=[]
    moreCount=[]
    i=0
    kt = cKDTree(poly[:-1])
    for pt in poly[:-1]:
        coincident=kt.query_ball_point(pt,tolerance)
        n=len(coincident)
        i+=1
        if n==1:
            oneCount.append((pt[0],pt[1], str(i-1)))
        if n==2:
            twoCount.append(pt)
        elif n==3:
            threeCount.append(pt)
        elif n==4:
            fourCount.append(pt)
        elif n>4:
            moreCount.append(pt)
    plt.figure(figsize=(8,5))
    x,y=zip(*poly)
    plt.plot(x,y, c='gray', linewidth=2.0, zorder=1)
    plt.scatter(x,y, s=30, c='b', zorder=1)
    plt.title(label)
    if len(twoCount) > 0:
        x2,y2=zip(*twoCount)
        plt.scatter(x2,y2, s=100, zorder=2, facecolors='None', edgecolors='r')
        plt.scatter(x2,y2, s=300, zorder=2, facecolors='None', edgecolors='r')
    if len(threeCount) > 0:
        x3,y3=zip(*threeCount)
        plt.scatter(x3,y3, s=100, zorder=3, facecolors='None', edgecolors='r')
        plt.scatter(x3,y3, s=300, zorder=3, facecolors='None', edgecolors='r')
        plt.scatter(x3,y3, s=500, zorder=3, facecolors='None', edgecolors='r')
    if len(fourCount) > 0:
        x4,y4=zip(*fourCount)
        plt.scatter(x4,y4, s=100, zorder=4, facecolors='None', edgecolors='r')
        plt.scatter(x4,y4, s=300, zorder=4, facecolors='None', edgecolors='r')
        plt.scatter(x4,y4, s=500, zorder=4, facecolors='None', edgecolors='r')
        plt.scatter(x4,y4, s=700, zorder=4, facecolors='None', edgecolors='r')
    if len(moreCount) > 0:
        xmore,ymore=zip(*moreCount)
        plt.scatter(xmore,ymore, s=500, zorder=5)
    plt.scatter([x[0]],[y[0]], s=100, c='b', zorder=6)
    plt.scatter([x[1]],[y[1]], s=65, c='y', zorder=6)
    # labels
    ymin = min(y)
    ymax = max(y)
    yoffset = (ymax-ymin)/15
    plt.axis('equal')
    ptdict={}
    for i in range(len(poly)):
        index = ptdict.get(poly[i],-1)
        if index==-1:
            ptdict[poly[i]]=poly[i][1] #+yoffset
        else:
            ptdict[poly[i]]=ptdict[poly[i]] #+yoffset
        txt=str(i)
        if pt_vals !=None:
            txt += " ({})".format(pt_vals[i])
        plt.annotate(txt,(poly[i][0],ptdict[poly[i]]))

def bounding_box(p1,p2,buffer=0):
    """Returns left, bottom, right and top coordinate values"""
    if p1[0] < p2[0]:
        left=p1[0]
        right=p2[0]
    else:
        left=p2[0]
        right=p1[0]
    if p1[1] < p2[1]:
        bottom=p1[1]
        top=p2[1]
    else:
        bottom=p2[1]
        top=p1[1]
    return (left-buffer, bottom-buffer, right+buffer, top+buffer)

def segment_bounding_box_generator(pts, other_pt_lists = None):
    """Returns a generator that yields tuples of (id, bounding box, None)
       for all of the segments between consecutive points in pts.
       This is used to create the R-Tree index."""
    segs=list(zip(pts[:-1],pts[1:]))
    for i, obj in enumerate(list(zip(range(len(segs)),segs))):
        # Get left, bottom, right and top coordinate values
        seg=obj[1]
        p1=seg[0]
        p2=seg[1]
        left, bottom, right, top=bounding_box(p1,p2)
        # Yield tuple as desired by rtree (sequence: LBRT)
        yield (i,(left, bottom, right, top),None)
    # do the same for the other lines
    if other_pt_lists != None:
        cur_id = len(pts)*2+1
        for pt_list in other_pt_lists:
            segs=list(zip(pt_list[:-1],pt_list[1:]))
            for i, obj in enumerate(list(zip(range(len(segs)),segs))):
                # Get left, bottom, right and top coordinate values
                seg=obj[1]
                p1=seg[0]
                p2=seg[1]
                left, bottom, right, top=bounding_box(p1,p2)
                # Yield tuple as desired by rtree (sequence: LBRT)
                yield (cur_id,(left, bottom, right, top),None)
                cur_id +=1