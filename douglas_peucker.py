'''
This is a very simple implementation of 
the Douglas Peucker line simplification algorithm.
It is not fast, nor does it include any checks to avoid self-intersection.

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

import numpy as np

#from numpy import cross,array, array_equal
#from numpy.linalg import norm

def __distance_pt_to_line(pt,line_start, line_end):
    '''calculates the perpendicular distance from the point to the line'''
    # convert to numpy arrays
    a=np.array(line_start)
    b=np.array(line_end)
    x=np.array(pt)
    if np.array_equal(a,b): # coincident start/end points
        d=np.linalg.norm(a-x)
    else: # ordinary case
        d=np.linalg.norm(np.cross(b-a,a-x))/np.linalg.norm(b-a)
    return d

def __calc_errors_recursive(pts,s,e,errors):
    '''Finds the point furthest from the line from pt s to pt e
    and marks it with a tolerance equal to the distance, or the last tolerance
    if it is greater.
    '''
    # initialize
    maxd=0
    furthest_pt=-1
    # loop
    for i in range(s+1,e):
        # get distance from pt i to line from s to e
        d = __distance_pt_to_line(pts[i],pts[s],pts[e])
        # update
        if d > maxd:
            maxd=d
            furthest_pt=i
    # check if we have any results
    if furthest_pt != -1:
        # update tolerance of furthest point
        errors[furthest_pt]=maxd
        # continue recursively
        maxleft=__calc_errors_recursive(pts,s,furthest_pt,errors)
        maxright=__calc_errors_recursive(pts,furthest_pt,e,errors)
        # make sure tolerance is at least as high as all child points
        if maxd < maxleft:
            maxd = maxleft
        if maxd < maxright:
            maxd = maxright
        # mark tolerance
        errors[furthest_pt] = maxd
        # return tolerance
        return maxd
    else:
        return -1


def __calc_errors(pts):
    '''calculates the tolerance above which each point should be retained'''
    errors=[0]*len(pts)
    __calc_errors_recursive(pts,0,len(pts)-1,errors)
    maxe = max(errors)
    errors[0]=maxe+1.0
    errors[-1]=maxe+1.0
    return errors

def simplify_by_error_tolerance(pts, errors, tolerance):
    '''returns points above given tolerance'''
    r=[]
    for i in range(len(pts)):
        if errors[i]>=tolerance:
            r.append(pts[i])
    return r

def get_errors_sortedErrors(pts):
    errors=__calc_errors(pts)
    sorted_errors=sorted(errors,reverse=True)
    return errors, sorted_errors

def fast_simplify_by_numPts(pts,numpts,errors, sorted_errors):
    t=sorted_errors[numpts-1]
    return simplify_by_error_tolerance(pts,errors,t)

def simplify_by_numpts(pts,numpts):
    '''returns simplified point list with given number of points'''
    r=[]
    # determine error tolerance needed to include given # points
    errors=__calc_errors(pts)
    sorted_errors=sorted(errors,reverse=True)
    t=sorted_errors[numpts-1]
    # get line
    return simplify_by_error_tolerance(pts,errors,t)

def main():
    # test on sample line
    sample_line=[]
    sample_line.append([0.0,0.0])
    sample_line.append([0.0,1.0])
    sample_line.append([1.0,2.0])
    sample_line.append([2.0,2.0])
    sample_line.append([3.0,1.0])
    sample_line.append([2.0,-2.0])
    sample_line.append([3.0,-4.0])
    sample_line.append([-1.0,-1.0])
    sample_line.append([0.0,0.0])

    errors=__calc_errors(sample_line)
    for t in errors:
        print(t)
    simple1=simplify_by_error_tolerance(sample_line,errors,1.0)
    print("by error:")
    print (simple1)
    simple3=simplify_by_numpts(sample_line,errors,3)
    print("by number of points:")
    print (simple3)

if __name__ == '__main__':
    main()
