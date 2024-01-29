#-------------------------------------------------------------------------------
# Name:        line_intersection
# Purpose:     utilities for detecting and fixing line intersections
# Author:      Barry Kronenfeld
# Created:     17/01/2018
# Licence:     MIT License
#-------------------------------------------------------------------------------

from .geom_utils import isMonotonic,distance_pts,triangle_area,ptsMonotonic,area


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



def __overlaps(a1,a2,b1,b2):
    """Returns True if ranges of a and b overlap
       Returns False otherwise"""
    if isMonotonic(a1,b1,a2):
        return True
    elif isMonotonic(a1,b2,a2):
        return True
    elif isMonotonic(b1,a1,b2):
        return True
    elif isMonotonic(b1,a2,b2):
        return True
    else:
        return False



def __segmentsIdentical(s1,s2):
    a=s1[0]
    b=s1[1]
    c=s2[0]
    d=s2[1]
    """Determines if segments ab and cd are identical, in any sequence"""
    if a==c and b==d:
        return True
    if a==d and b==c:
        return True
    return False

def points_sorted(p):
    """Takes a list of input points, and returns 
    a list of the same points sorted along a line between
    the furthest pair.
    Args:
        p - a list of (x,y) tuples"""
    if len(p)<3:
        return p
    else:
        # find furthest pair of points
        fpair=(0,1)
        fdist=distance_pts(p[0],p[1])
        for i in range(len(p)):
            for j in range(i+1,len(p)):
                d=distance_pts(p[i],p[j])
                if d > fdist:
                    fdist=d
                    fpair=(i,j)
        # sort points along line from first to last
        start=p[fpair[0]]
        p_dist=[(x,distance_pts(x,start)) for x in p]
        p_dist=sorted(p_dist,key = lambda x:x[1])
        return [x[0] for x in p_dist]        
        
def __segmentsOverlap(s1,s2,t):
    """Determines if segments ab and cd overlap, with distance tolerance t.
    Returns None if no overlap, otherwise returns list of two endpoints of overlap"""
    a=s1[0]
    b=s1[1]
    c=s2[0]
    d=s2[1]
    # triangle areas must be close to zero to overlap
    if abs(triangle_area(a,b,c)/distance_pts(a,b)) > t:
        return None
    if abs(triangle_area(a,b,d)/distance_pts(c,d)) > t:
        return None
    if abs(triangle_area(a,c,d)/distance_pts(a,b)) > t:
        return None
    if abs(triangle_area(b,c,d)/distance_pts(c,d)) > t:
        return None
    # x- and y-coordinates must be monotonic to overlap
    if not __overlaps(a[0],b[0],c[0],d[0]):
        return None
    if not __overlaps(a[1],b[1],c[1],d[1]):
        return None
    # if vertices are identical, there may not be overlap
    if a==c and ptsMonotonic(b,a,d):
        return None
    if a==d and ptsMonotonic(b,a,c):
        return None
    if b==c and ptsMonotonic(a,b,d):
        return None
    if b==d and ptsMonotonic(a,d,c):
        return None
    # can't sort just by x coordinates
    sorted_pts = points_sorted([a,b,c,d])
    overlap_endpoints = sorted_pts[1:3]
    return overlap_endpoints

def overlappingNonIdenticalSegments(polyline, t):
    """Searches for pairs of line segments that overlap within distance tolerance t. Returns two lists:
    (1) list of intersecting segment id pairs
    (2) list of lists of intersection points. most inner lists will have
        just one point. an inner list with two points means a pair of
        segments that overlaps."""
    # INITIALIZE RESULT LISTS
    intPairs=[]
    intPts=[]
    # FIND INTERSECTIONS BETWEEN NON-ADJACENT LINE SEGMENTS
    # loop through all segments
    for i in range(len(polyline)-1):
        seg1=polyline[i:i+2]
        # loop through all higher segments except the next one
        for j in range(i+2, len(polyline)-1):            
            seg2=polyline[j:j+2]
            # search for overlap          
            if not __segmentsIdentical(seg1,seg2):
                overlap = __segmentsOverlap(seg1,seg2,t)
                if not overlap == None:                   
                    # we have liftoff!
                    intPairs.append([i,j])
                    intPts.append(overlap)
    return list(zip(intPairs,intPts))

def segmentsIntersect_triangleMethod(seg1,seg2):
    """Attempts to check if triangles intersect by looking at how many 
    triangles have positive and negative areas.
    """
    # get variables
    a = seg1[0]
    b = seg1[1]
    c = seg2[0]
    d = seg2[1]
    ABC = triangle_area(a,b,c)
    ABD = triangle_area(a,b,d)
    ACD = triangle_area(a,c,d)
    BCD = triangle_area(b,c,d)
    ABC_ABD_opposite = (ABC >= 0) != (ABD >= 0)
    ACD_BCD_opposite = (ACD >= 0) != (BCD >= 0)
    if ABC_ABD_opposite and ACD_BCD_opposite:
        return True
    else:
        return False
    
def lineIntersection(seg1, seg2):
    """Returns a list containing a single intersection point between lines going through segments
    If segments do not intersect, too bad we'll return a point anyway!
    (based on Stephen Wise, GIS Basics, CRC Press, pp. 48-9)"""

    # get variables
    xs1=seg1[0][0]
    ys1=seg1[0][1]
    xe1=seg1[1][0]
    ye1=seg1[1][1]
    xs2=seg2[0][0]
    ys2=seg2[0][1]
    xe2=seg2[1][0]
    ye2=seg2[1][1]

    # work through cases
    if xs1==xe1:
        if xs2==xe2: # both lines vertical
            if xs1==xs2:
                # intersection will be middle two points in sequency of y-coordinates
                y = [ys1,ys2,ye1,ye2]
                y.sort()
                return[(xs1,y[1])]
        else: # first line vertical
            # determine y-coordinate of intersection
            b2=(ye2-ys2)/(xe2-xs2)
            a2=ys2-b2*xs2
            py=a2+b2*xs1
            # test if intersection is on vertical line
            return [(xs1,py)]
    elif xs2==xe2: # second line vertical
        # determine y-coordinate of intersection
        b1=(ye1-ys1)/(xe1-xs1)
        a1=ys1-b1*xs1
        py=a1+b1*xs2
        # test if intersection is on vertical line
        return [(xs2,py)]
    else: # neither line vertical
        b1=(ye1-ys1)/(xe1-xs1)
        b2=(ye2-ys2)/(xe2-xs2)
        a1=ys1-b1*xs1
        a2=ys2-b2*xs2
        if b1==b2: # lines are parallel                
            # intersection will be middle two points in sequence of either variable
            ptseq = sorted(list(zip([xs1,xe1,xs2,xe2],[ys1,ye1,ys2,ye2])),key=lambda pair: pair[0])
            return ptseq[1:2]
        else:
            # get possible intersection point
            px=-(a1-a2)/(b1-b2)
            py=a1+b1*px       
            return [(px,py)]

def segmentIntersection(seg1, seg2):
    """Returns list of one intersection point if line segments intersect.
    Returns list of two points if line segments overlap.
    Otherwise returns None.
    (based on Stephen Wise, GIS Basics, CRC Press, pp. 48-9)"""

    # identical segments don't intersect
    if __segmentsIdentical(seg1,seg2):
        return None
    
    # segments that share a vertex do intersect
    a=seg1[0]
    b=seg1[1]
    c=seg2[0]
    d=seg2[1]
    if a==c or a==d:
        return[a]
    if b==c or b==d:
        return[b]

    # get variables
    xs1=seg1[0][0]
    ys1=seg1[0][1]
    xe1=seg1[1][0]
    ye1=seg1[1][1]
    xs2=seg2[0][0]
    ys2=seg2[0][1]
    xe2=seg2[1][0]
    ye2=seg2[1][1]

    # pretest - bounding boxes (improves speed by ~27%)
    if xs1 < xe1:
        x1min=xs1
        x1max=xe1
    else:
        x1min=xe1
        x1max=xs1
    if xs2 < xe2:
        x2min=xs2
        x2max=xe2
    else:
        x2min=xe2
        x2max=xs2
    if ys1 < ye1:
        y1min=ys1
        y1max=ye1
    else:
        y1min=ye1
        y1max=ys1
    if ys2 < ye2:
        y2min=ys2
        y2max=ye2
    else:
        y2min=ye2
        y2max=ys2
    if x1min > x2max:
        return None
    if x2min > x1max:
        return None
    if y1min > y2max:
        return None
    if y2min > y1max:
        return None

    # work through cases
    if xs1==xe1:
        if xs2==xe2: # both lines vertical
            if xs1==xs2:
                if __overlaps(ys1,ye1,ys2,ye2):
                    # intersection will be middle two points in sequency of y-coordinates
                    y = [ys1,ys2,ye1,ye2]
                    y.sort()
                    if y[1]==y[2]:
                        return[(xs1,y[1])]
                    else:
                        return [(xs1,y[1]),(xs1,y[2])]
                else:
                    return None
        else: # first line vertical
            # test if intersection can possibly be on vertical line
            if isMonotonic(xs2,xs1,xe2):
                # determine y-coordinate of intersection
                b2=(ye2-ys2)/(xe2-xs2)
                a2=ys2-b2*xs2
                py=a2+b2*xs1
                # test if intersection is on vertical line
                if isMonotonic(ys1,py,ye1):
                    return [(xs1,py)]
                else:
                    return None
            else:
                return None
    elif xs2==xe2: # second line vertical
        # test if intersection can possibly be on vertical line
        if isMonotonic(xs1,xs2,xe1):
            # determine y-coordinate of intersection
            b1=(ye1-ys1)/(xe1-xs1)
            a1=ys1-b1*xs1
            py=a1+b1*xs2
            # test if intersection is on vertical line
            if isMonotonic(ys2,py,ye2):
                return [(xs2,py)]
            else:
                return None
        else:
            return None
    else: # neither line vertical
        b1=(ye1-ys1)/(xe1-xs1)
        b2=(ye2-ys2)/(xe2-xs2)
        a1=ys1-b1*xs1
        a2=ys2-b2*xs2
        if b1==b2: # lines are parallel
            if a1==a2: # intercepts are equal
                if __overlaps(xs1,xe1,xs2,xe2):
                    # intersection will be middle two points in sequence of either variable
                    ptseq = sorted(list(zip([xs1,xe1,xs2,xe2],[ys1,ye1,ys2,ye2])),key=lambda pair: pair[0])
                    if ptseq[1]==ptseq[2]:
                        return ptseq[1:2]
                    else:
                        return ptseq[1:3]
                else:
                    return None
            else:
                return None
        else:
            # get possible intersection point
            px=-(a1-a2)/(b1-b2)
            py=a1+b1*px
            # test if it is on both lines
            if isMonotonic(xs1,px,xe1) and isMonotonic(xs2, px, xe2):
                return [(px,py)]
            else:
                return None


def __adjacent_segment_intersections(polyline):
    """Searches for overlap between adjacent line segments in
    a polyline. Returns two lists:
    (1) list of intersecting segment id pairs
    (2) list of lists of intersection points. most inner lists will have
        just one point. an inner list with two points means a pair of
        segments that overlaps.
        """
    pl=polyline
    # init
    intPairs = []
    intPts = []
    # loop through vertex triplets
    for i in range(len(polyline)-2):
        # test for collinearity
        a=area(polyline[i:i+3])
        if a==0:
            # test for sequence
            if isMonotonic(pl[i][0],pl[i+2][0],pl[i+1][0]):
                if isMonotonic(pl[i][1],pl[i+2][1],pl[i+1][1]):
                    # 3rd point is between 1st and 2nd
                    intPairs.append([i,i+1])
                    intPts.append([pl[i+1],pl[i+2]])
            if isMonotonic(pl[i+1][0],pl[i][0],pl[i+2][0]):
                if isMonotonic(pl[i+1][1],pl[i][1],pl[i+2][1]):
                    # 1st point is between 2nd and 3rd
                    intPairs.append([i,i+1])
                    intPts.append([pl[i],pl[i+1]])
    return intPairs, intPts


def self_intersections_brute_force_triangle_method(polyline):
    """Searches every pair of line segments to search
    for intersections. Returns two lists:
    (1) list of intersecting segment id pairs
    (2) list of lists of intersection points. inner list will have
        just one point. 
    Note that intersections that occur at a vertex will be double- or quadruple-counted.
    """
    # INITIALIZE RESULT LISTS
    intPairs=[]
    intPts=[]
    # FIND INTERSECTIONS BETWEEN NON-ADJACENT LINE SEGMENTS
    # loop through all segments
    for i in range(len(polyline)-1):
        seg1=polyline[i:i+2]
        # loop through all higher segments except the next one
        for j in range(i+2, len(polyline)-1):
            seg2=polyline[j:j+2]
            # first check identity
            if not __segmentsIdentical(seg1,seg2):
                # now check overlap
                if segmentsIntersect_triangleMethod(seg1,seg2):
                    # get intersection point
                    intPt = lineIntersection(seg1,seg2)
                    intPairs.append([i,j])
                    intPts.append(intPt)
    return list(zip(intPairs,intPts))

    
def self_intersections_brute_force(polyline):
    """Searches every pair of line segments to search
    for intersections. Returns two lists:
    (1) list of intersecting segment id pairs
    (2) list of lists of intersection points. most inner lists will have
        just one point. an inner list with two points means a pair of
        segments that overlaps.
    Note that intersections that occur at a vertex will be double- or quadruple-counted.
    """
    # INITIALIZE RESULT LISTS
    intPairs=[]
    intPts=[]
    # FIND INTERSECTIONS BETWEEN NON-ADJACENT LINE SEGMENTS
    # loop through all segments
    for i in range(len(polyline)-1):
        seg1=polyline[i:i+2]
        # loop through all higher segments except the next one
        for j in range(i+2, len(polyline)-1):
            seg2=polyline[j:j+2]
            # *** debugging
            if i==6:
                if j==10:
                    dummy=True
            
            
            
            # *** first check overlap
            overlap = __segmentsOverlap(seg1,seg2,0.000001)
            if not overlap == None:
                # we have liftoff!
                intPairs.append([i,j])
                intPts.append(overlap)
            else:
            
                # get intersection
                thisInt = segmentIntersection(seg1,seg2)
                # if we caught something, add to list:
                if thisInt != None:
                    intPairs.append([i,j])
                    intPts.append(thisInt)
    # FIND OVERLAP BETWEEN ADJACENT LINE SEGMENTS
    adjIntPairs, adjIntPts = __adjacent_segment_intersections(polyline)
    intPairs = intPairs + adjIntPairs
    intPts = intPts + adjIntPts
    return list(zip(intPairs,intPts))

def self_intersections(polyline):
    """Searches for pairs of line segments that intersect using a line
    sweep algorithm. Returns list of tuples:
    (1) list of intersecting segment id pairs, each with two items:
    (2) list of lists of intersection points. most inner lists will have
        just one point. an inner list with two points means a pair of
        segments that overlaps.
    Note that intersections that occur at a vertex will be double- or quadruple-counted.
    Note that this is NOT the Bentley-Ottman algorithm,
    and is worst-case O(n^2) but it should usually be much faster (e.g. O(n^1.5))
    """
    # init
    intPairs=[]
    intPts=[]
    # SWEEP LINE SEARCH
    # sort vertices by x values and create double-lookup
    # v   : vertex at position i
    # pos : position of vertex i
    ids=range(len(polyline)) # just a list of IDs from 0 to n-1
    plids = list(zip(polyline,ids)) # a list of (point, id) tuples, in order of the input polyline
    sortedplids = sorted(plids, key = lambda pair: pair[0][0]) # a list of (point, id) tuples, in order of x-coord
    v = list(zip(*sortedplids)[1]) # list of ids, sorted by x-coord
    pos=[-1]*len(v) # list of indexes in v; if v[3]=7 then pos[7]=3
    for i in range(len(v)):
        pos[v[i]]=i
    # Local Set: 0, 1 or 2 segments extending right or vertically up from current vertex
    # Wide Set:  all segments with earlier-sequence vertices that
    #            might intersect a segment in local set

    # initialize wide set
    wide_set=set()
    # loop through vertices in order of x
    for cur_v in v:
        # get segments right of current vertex
        local_set=[]
        # add segments extending right to local set
        # check 2 adjacent segments, entering and exiting current vertex
        # (entering segment)
        if cur_v > 0:
            if pos[cur_v-1] > pos[cur_v]:
                local_set.append(cur_v-1)
        # (exiting segment)
        if cur_v < len(polyline)-1:
            if pos[cur_v+1] > pos[cur_v]:
                local_set.append(cur_v)
        # check for intersections between segs in local & wide sets
        for seg1 in local_set:
            for seg2 in wide_set:
                # make sure they are not adjacent
                if abs(seg1-seg2) >=2:
                    # get intersection
                    thisInt = segmentIntersection(polyline[seg1:seg1+2],polyline[seg2:seg2+2])
                    # if we caught something, add to list:
                    if thisInt != None:
                        intPairs.append(sorted([seg1,seg2]))
                        intPts.append(thisInt)
        # remove segments extending left from wide set
        # (entering segment)
        if cur_v > 0:
            if pos[cur_v-1] < pos[cur_v]:
                wide_set.remove(cur_v-1)
        # (exiting segment)
        if cur_v < len(polyline)-1:
            if pos[cur_v+1] < pos[cur_v]:
                wide_set.remove(cur_v)
        # add local set to wide set
        wide_set.update(local_set)
    # ADJACENT SEGMENT INTERSECTIONS
    adjIntPairs, adjIntPts = __adjacent_segment_intersections(polyline)
    intPairs = intPairs + adjIntPairs
    intPts = intPts + adjIntPts
    return list(zip(intPairs,intPts))


def intersections(polyline1, polyline2):
    """Searches for pairs of line segments that intersect using a line
    sweep algorithm. Returns a list of tuples, where each tuple contains three values:
    (0) ID of start vertex of segment from polyline 1
    (1) ID of start vertex of segment from polyline 2
    (2) point of intersection, as (x,y) tuple
    Note that intersections that occur at a vertex will be double- or quadruple-counted.
    Note that this is NOT the Bentley-Ottman algorithm,
    and is worst-case O(n^2) but it should usually be much faster (e.g. O(n^1.5))
    """
##    # init
##    intPairs=[]
##    merged=polyline1[:]
##    merged.extend(polyline2)
##    L1=len(polyline1)
##    # SWEEP LINE SEARCH
##    # sort vertices by x values and create double-lookup
##    # v   : vertex at position i
##    # pos : position of vertex i
##    ids=range(len(merged))
##    plids = list(zip(merged,ids))
##    sortedplids = sorted(plids, key = lambda pair: pair[0][0])
##    v = list(zip(*sortedplids)[1])
##    pos=[-1]*len(v) # initialize all positions to -1
##    for i in range(len(v)):
##        pos[v[i]]=i
##    # Local Set: 0, 1 or 2 segments extending right or vertically up from current vertex
##    # Wide Set:  all segments with earlier-sequence vertices that might
##    #            might intersect a segment in local set
##
##    # initialize local set
##    wide_set=set()
##    # loop through vertices in order of x
##    for cur_v in v:
##        # get segments right of current vertex
##        local_set=[]
##        # add segments extending right to local set
##        # (entering segment)
##        if cur_v > 0:
##            if pos[cur_v-1] > pos[cur_v]:
##                local_set.append(cur_v-1)
##        # (exiting segment)
##        if cur_v < len(merged)-1:
##            if pos[cur_v+1] > pos[cur_v]:
##                local_set.append(cur_v)
##        # check for intersections between segs in local & wide sets
##        for seg1 in local_set:
##            for seg2 in wide_set:
##                # make sure they are not adjacent
##                if abs(seg1-seg2) >=2:
##                    # make sure they are not the false connection between polylines 1 & 2
##                    if seg1 != L1-1:
##                        if seg2 != L1-1:
##                            # make sure they are on different polylines
##                            if (seg1 < L1) != (seg2 < L1):
##                                # get intersection
##                                thisInt = segmentIntersection(merged[seg1:seg1+2],merged[seg2:seg2+2])
##                                # if we caught something, add to list:
##                                if thisInt != None:
##                                    segIDs = sorted([seg1,seg2])
##                                    intPairs.append((segIDs[0],segIDs[1]-L1,thisInt))
##        # remove segments extending left from wide set
##        # (entering segment)
##        if cur_v > 0:
##            if pos[cur_v-1] < pos[cur_v]:
##                wide_set.remove(cur_v-1)
##        # (exiting segment)
##        if cur_v < len(merged)-1:
##            if pos[cur_v+1] < pos[cur_v]:
##                wide_set.remove(cur_v)
##        # add local set to wide set
##        wide_set.update(local_set)
##    # return result (no need to worry about adjacent segments
##    return intPairs

def checkInsertIntersectionPoints(polyline, s1, s2):
    '''
    Checks to see if one or more points should be inserted 
    on segments 1 or 2 to mark intersections
    args:
        polyline: list of (x,y) tuples
        s1: index of start point of first segment
        s2: index of start point of second segment
    returns:
        list of (id, x, y) tuples denoting new points to insert into polyline
        in reverse sequence of order they should be inserted
    '''
    # determine if point can be considered to be coincident with other segment
    
    pass
