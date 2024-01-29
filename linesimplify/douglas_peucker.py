#-------------------------------------------------------------------------------
# Name:        Douglas Peucker
# Description: A very simple implementation of the Ramer-Douglas-Peucker 
#              line simplification algorithm. This does not check for 
#              topological errors. The implementation here is not particularly
#              fast, but it does support pre-processing to quickly simplify
#              to any given distance tolerance or number of points. 
#
#              References:
#              
#              Urs Ramer (1972). An iterative procedure for the polygonal approximation 
#              of plane curves. Computer Graphics and Image Processing, vol. 1, pp. 244-256.
#
#              David H. Douglas and Thomas K. Peucker (1973). Algorithms for the reduction of 
#              the number of points required to represent a digitized line or its caricature.
#              The Canadian Cartographer, vol. 10, no. 2, pp. 112-122.
#
# Author:      Barry Kronenfeld
# License:     MIT License
#-------------------------------------------------------------------------------

import numpy as np


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


def get_errors_sortedErrors(pts):
    errors=__calc_errors(pts)
    sorted_errors=sorted(errors,reverse=True)
    return errors, sorted_errors

def simplify_by_distance_tolerance(pts, errors, tolerance):
    '''returns points above given tolerance'''
    r=[]
    for i in range(len(pts)):
        if errors[i]>=tolerance:
            r.append(pts[i])
    return r

def simplify_by_numPts(pts,numpts,errors, sorted_errors):
    t=sorted_errors[numpts-1]
    return simplify_by_distance_tolerance(pts,errors,t)

