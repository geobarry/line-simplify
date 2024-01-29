#-------------------------------------------------------------------------------
# This code implements the area preserving segment collapse (APSC) algorithm, 
# which simplifies a geographic line feature while preserving area and 
# minimizing areal displacement. The APSC algorithm is described in detail in:
#
# Kronenfeld, Stanislawski, Buttenfield and Brockmeyer (2020). Simplification
# of polylines by segment collapse: minimizing areal displacement while 
# preserving area. International Journal of Cartography, vol. 6, issue 1, 
# https://doi.org/10.1080/23729333.2019.1631535
#
# Code by Barry Kronenfeld (last updated June 2020)
# MIT License
#-------------------------------------------------------------------------------

"""
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import sys

from . import line_intersection as lineint

# import linesimplify.line_intersection as lineint
from . import geom_utils as g

try:
    import sortedcontainers as __sc
except ImportError:
    sys.exit("""You need sortedcontainers to run apsc.py.
                Install it from https://pypi.python.org/pypi/sortedcontainers
                or get pip (https://pip.pypa.io/en/stable/installing/) and run the following from a command prompt:
                      pip install sortedcontainers""")
try:
    from rtree import index
except ImportError:
    sys.exit("""You need Rtree to run apsc.py.
                To install it:
                # 1) go to https://www.lfd.uci.edu/~gohlke/pythonlibs/#rtree
                # 2) download the appropriate one of the following to the python scripts folder:
                #       Rtree-0.8.3-cp27-cp27m-win32.whl
                #       Rtree-0.8.3-cp27-cp27m-win_amd64.whl
                # 3) open a command-prompt as an administrator
                # 4) navigate to python scripts folder
                # 5) uninstall any previous version, just in case, by entering:
                #      pip uninstall rtree
                # 6) install by entering:
                #       pip install Rtree-0.8.3-cp27-cp27m-win32.whl""")

def __packUnpackSimplificationTableExample():
    # This function is not used anywhere
    # it merely provides examples to illustrate the structure of
    # the data in the APSC simplification table
    messages=[]
    messages.append('#The simplification "table" is a list of lists,best conceptualized as a table with 9 columns.')
    messages.append('i.e. s_tab= [ID,XY,error,parent,LC,RC,LS,RS,current]')
    messages.append('#Each inner list contains a column of data with 2n-3 elements.')
    messages.append('The following code obtaines these lists:')
    messages.append('ID,XY,error,parent,LC,RC,LS,RS,currents=s_tab[0],s_tab[1],s_tab[2],s_tab[3],s_tab[4],s_tab[5],s_tab[6],s_tab[7],s_tab[8]')
    messages.append('#The ID column contains the ID of each vertex')
    messages.append("#This column isn't really necessary, as it's always the same as the list index. So, for example, the following is always true:")
    messages.append('    ID[0]==0')
    messages.append('#The XY column contains the location of each vertex, as a tuple. For example:')
    messages.append('    XY[0]=(21.2,36.7)')
    messages.append("#The error column contains the displacement area, or other error metric, associated with collapsing a vertex's children.")
    messages.append('#Since leaf nodes have no children, their error will always be 0.0, as in:')
    messages.append('    error[0]=0')
    messages.append('# If the original line contains n vertices, the first n errors will be 0, so the following will generally be true:')
    messages.append('    error[n-1]==0')
    messages.append('    error[n] > 0')
    messages.append('# The parent column contains the parent ID of each vertex, e.g.')
    messages.append('    parent[1]=27')
    messages.append('# The LC and RC columns contain the indices of the left and right children, as in:')
    messages.append('    LC[27]=1')
    messages.append('# The LS and RS columns contain the indices of the left and right siblings')
    messages.append('# The first n vertices {i=0...n-1} are siblings with their neighbors, i.e. the following are true:')
    messages.append('    if i < n: LS[i]==i-1')
    messages.append('    if i < n-1: RS[i]==i+1')
    messages.append('# The "current" column is used during tree construction, and contains True for vertices that have not yet been collapsed.')
    messages.append('# after tree construction has finished, if every internal segment (excluding the two end segments) was able to be collapsed then only three vertices will be current.')
    return messages


def __sideshift_area(pA,pB,pC,pD,pE,overlap_endpt):
    """Returns the sideshift displacement area from ABCD to ACE.
       Args:
           pA, pB, pC, pD, pE: (x,y) tuples
           overlap_endpt:
                0      E is placed on AB)
                1      E is placed on CD)
                other  E is not placed on AB or CD)
       Returns:
            Total sideshift displacement area between ABCD and AED"""
    if g.intersect(pA,pB,pC,pD):  # original line self-intersects
        x=g.intersection(pA,pB,pC,pD)
        return g.area([x,pC,pB],True) + g.area([x,pA,pE,pD],True)
    elif overlap_endpt==0: # E is on line through AB
        if g.intersect(pB,pC,pE,pD):
            x=g.intersection(pB,pC,pE,pD)
            return g.area([pB,pE,x],True) + g.area([x,pD,pC],True)
        else:
            return g.area([pD,pC,pB,pE],True)
    elif overlap_endpt==1:  # E is on line through CD
        if g.intersect(pC,pB,pE,pA):
            x=g.intersection(pC,pB,pE,pA)
            return g.area([pC,pE,x],True) + g.area([x,pA,pB],True)
        else:
            return g.area([pA,pB,pC,pE],True)
    else: # possibilities: (1) AE^BC, (2) AE^BC & ED^BC, (3) AE^CD, (4) DE^AB, (5) DE^BC, (6) no intersections
        # check if AE intersects BC or CD
        if g.intersect(pA,pE,pB,pC):   # AE intersects BC
            x1=g.intersection(pA,pE,pB,pC)
            if g.intersect([pD, pE, pB, pC]):  # (2)
                x2=g.intersection(pD,pE,pB,pC)
                return g.area([pA,x1,pB],True)+g.area([x1,pE,x2],True)+g.area([x2,pD,pC],True)
            else:                             # (1)
                return g.area([pA,x1,pB],True) + g.area([x1,pE,pC,pD],True)
        elif g.intersect(pA,pE,pC,pD): # AE intersects CD  (3)
            x=g.intersection(pA,pE,pC,pD)
            return g.area([pA,x,pC,pB],True) + g.area([x,pE,pD],True)
        elif g.intersect(pD,pE,pA,pB): # DE intersects AB  (4)
            x=g.intersection(pD,pE,pA,pB)
            return g.area([pA,pE,x],True) + g.area([x,pD,pC,pB],True)
        elif g.intersect(pD,pE,pB,pC): # DE intersects BC  (5)
            x=g.intersection(pD,pE,pB,pC)
            return g.area([pA,pE,x,pB],True) + g.area([x,pD,pC],True)
        else: # no intersection (6)
            return g.area([pA,pE,pD,pC,pB],True)

def __placement_AP_EAmin(pA,pB,pC,pD, m_func):
    """Determines point E to replace BC that minimizes effective area of displacement while preserving area on each side.
    Args:
        pA,pB,pC,pD: xy tuples representing sequence of four vertices
        m_func: function to calculate displacement metric
    Returns:
        pE: location of new point E
        displacement: (half) area of displacement between ABCD and AED
        overlap_endpt: 0  (E is placed on AB)
                       1  (E is placed on CD)
                       -1 (E is not placed on AB or CD)
    Simultaneously returns error metric (1/2 of displacement area)."""
    if pA==pD:
        return None, -1, -1
    # get symmetry line
    S=__equalAreaLine(pA,pB,pC,pD)
    # determine situation
    if S==None:
        return None, -1, -1
    elif g.distance_pt_line(pA,S[0],S[1])==0:
        return None, -1, -1
    else:
        # default: intersect S with CD
        C2=pC
        D2=pD
        overlap_endpt=1
        if g.X_Right_of_AB(pA,pD,pB)==g.X_Right_of_AB(pA,pD,pC):
            B_further_from_AD = g.distance_pt_line(pB,pA,pD) - g.distance_pt_line(pC,pA,pD)
            if B_further_from_AD > 0: # intersect S with AB
                D2=pA
                C2=pB
                overlap_endpt=0
##            if B_further_from_AD == 0: # intersect S with midpoints (arbitrary)
##                D2=[(pA[0]+pD[0])/2,(pA[1]+pD[1])/2]
##                C2=[(pB[0]+pC[0])/2,(pB[1]+pC[1])/2]
##                overlap_endpt=-1
        else:
            if g.X_Right_of_AB(pA,pD,pB) == g.X_Right_of_AB(pA,pD,S): # intersect S with AB
                D2=pA
                C2=pB
                overlap_endpt=0
        # calculate intersection
        pE=g.intersection(D2,C2,S[0],S[1])
        # handle case of colinear points
        if pE[0]==None:
            pE=pC
        # calculate displacement area
        if overlap_endpt==1: # intersect S with DC
            area=m_func(pA,pB,pC,pD,pE, overlap_endpt)
        else: # intersect S with AB
            area=m_func(pD,pC,pB,pA,pE, overlap_endpt)
        return pE,area, overlap_endpt

def __causesSelfIntersection(A,B,C,D,pE, overlap_endpt,LS_A, vtree, idx):
    """Determines if the collapse of BC to E would create a self-intersection,
    assuming that E is placed on CD and that idx is an index of vtree.
    Input:
        A,B,C,D,E: IDs of vertices involved in collapse
        LS_A: left sibling of A
        vtree: vertex simplification tree
        idx: spatial index
    Returns:
        True if it would, False otherwise."""
    # get XY coordinate array, right siblings
    XY=vtree[1]
    RS=vtree[7]
    # get points
    pA, pB, pC, pD=(XY[A],XY[B],XY[C],XY[D])
    # create segments to check
    newSegs=[(pA, pE),(pE, pD)]
    # get intersecting segments
    for newSeg in newSegs:
        int_seg_IDs=list(idx.intersection(__bounding_box(newSeg[0],newSeg[1])))
        # check segments in list to see if they really intersect newSeg
        for int_seg_ID in int_seg_IDs:
            # make sure current segment is not on the collapse sequence
            if not int_seg_ID in [A,B,C]:
                # get segment vertices
                int_seg =(XY[int_seg_ID],XY[RS[int_seg_ID]])
                # look for intersection
                ints=lineint.segmentIntersection(newSeg,int_seg)
                if ints != None: # we might have an intersection
                    # see if we are looking at segment prior to A or just after D
                    if RS[int_seg_ID]==A or int_seg_ID==D:
                        # if so, we have intersection only if there's an overlap ,
                        # indicated by two intersection points
                        if len(ints) > 1:
                            return True
                    else:
                        # if we have an intersection with any other segment, it is real
                        return True
    # if we've found no intersections by now, then we have none to report
    return False


def __collapse(B,C,E,thisxy,thiserror,vtree):
    """
    records the collapse of BC to E
    by adding a new row to the vertex relation table
    representing the new Steiner point E
    and adjusting the relations of B & C
    """
    # obtain relations as lists
    ID=vtree[0]
    XY=vtree[1]
    error=vtree[2]
    parent=vtree[3]
    LC=vtree[4]
    RC=vtree[5]
    LS=vtree[6]
    RS=vtree[7]
    current=vtree[8]
    # add new row
    # id, xy, error, parent, LC, RC, LS, RS, current
    ID.append(E)
    XY.append(thisxy)
    error.append(thiserror)
    parent.append(-1)
    LC.append(B)
    RC.append(C)
    LS.append(LS[B])
    RS.append(RS[C])
    current.append(True)

    # reset sibling relations to E
    RS[LS[B]] = E
    LS[RS[C]] = E
    # reset properties of B
    LS[B] = -1
    RS[B] = C
    parent[B] = E
    current[B] = False
    # reset properties of C
    LS[C] = B
    RS[C] = -1
    parent[C] = E
    current[C] = False

def __addToPriorityList(A,B,C,D,vtree,pl, p_func, m_func, do_topo_check=False):
    """
    calculates the collapse of BC to E (using the placement function p_func)
    and the associated error metric (using the metric function m_func)
    and adds the collapse event to the priority list pl
    in sorted ascending order of error.
    Args:
        A,B,C,D: Four vertex sequence
        vtree: list of vertex relations
        pl: priority list
        p_func: placement function
        m_func: metric funciton
        do_topo_check: If true, will avoid self-intersections (less efficient)
    Adds a tuple to the priority list with the following values:
        (displacement,A,B,C,D,pE,overlap_endpt)
    """
    XY=vtree[1]
    # make sure all vertices are current
    if A != -1 and B != -1 and C != -1 and D != -1:
        # compute collapse point E and associated displacement error
        pE,displacement, overlap_endpt = p_func(XY[A],XY[B],XY[C],XY[D],m_func)
        # make sure collapse point is valid
        if pE != None:
            # get index of overlapping endpoint
            if overlap_endpt==0:
                overlap_endpt=A
            if overlap_endpt==1:
                overlap_endpt=D
            pl.add((displacement,A,B,C,D,pE,overlap_endpt))

def __init_priority_list(vtree, p_func, m_func, do_topo_check=False):
    """
    Initializes the priority list of potential collapses to perform
    and populates the list with possible collapses on the original
    polyline.
    """
    ID=vtree[0]
    pl= __sc.SortedListWithKey(key=lambda x:x[0]) # empty sorted list
    for i in range(0,len(ID)-3): # loop through collapsible line segments
        __addToPriorityList(i,i+1,i+2,i+3,vtree,pl,p_func,m_func,do_topo_check)
    return pl # return sorted priority list

def __equalAreaLine(A,B,C,D):
    """ Returns two points on the line of symmetrical areal displacement.
    Replacing points B & C with any point E on this line
    will result in a new line feature AED whose displacement from ABCD
    is symmetrical, i.e. equal areal displacement on either side of the original line."""
    # obtain line in standard form: aX+bY+C=0
    a=D[1]-A[1]
    b=A[0]-D[0]
    c=-A[0]*B[1]+(A[1]-C[1])*B[0]+(B[1]-D[1])*C[0]+C[1]*D[0]
    # create two points on line near AD
    if a==0:
        if b==0:
            return None
        else:   # b>a so use X from A and D
            liney=-c/b
            p1=[A[0],liney]
            p2=[D[0],liney]
    elif b==0:
        linex = -c/a
        p1=[linex,A[1]]
        p2=[linex,D[1]]
        pass
    else:
        if abs(a)>abs(b): # use y
            p1=[(-c-b*A[1])/a,A[1]]
            p2=[(-c-b*D[1])/a,D[1]]
        else:             # use x
            p1=[A[0],(-c-a*A[0])/b]
            p2=[D[0],(-c-a*D[0])/b]
            useX=True
##        p1=[-c/a,0]
##        p2=[0,-c/b]
    # create line as list of points (as list of coords)
    sl=[p1,p2]
    # return to sender
    return sl

def __smallestLine(vtree):
    """Extracts the smallest possible line from the current vertex tree
    by tracing from the start to each right sibling. When the full simplification tree
    has been computed, the smallest possible line should consist of three points.
    Args:
        vtree: list of vertexRelations representing the current simplification tree.
    Returns:
        List of points (tuples of floats) representing the smallest line.
    """
    XY,RS=vtree[1],vtree[7]
    # initialize output as empty list
    outLine=[]
    # get first vertex relations and add point to output
    v=0
    outLine.append(XY[v])
    # loop through remaining points and repeat
    while RS[v] > -1:
        v=RS[v]
        outLine.append(XY[v])
    # return to sender
    return outLine

def __initVertexTree(pts):
    """Creates a list vertex relations associated with each point in input list.
    Each item in returned list is a vertex relation object with:
    <id,xypt,e, parent, LC, RC, LS, RS, current>"""
    ID=[]
    XY=[]
    error=[]
    parent=[]
    LC=[]
    RC=[]
    LS=[]
    RS=[]
    current=[]
    n=len(pts)
    # initialize default values for each vertex
    for i in range(n):
        # id, xy, error, parent, LC, RC, LS, RS, current
        ID.append(i)
        XY.append(pts[i])
        error.append(0)
        parent.append(-1)
        LC.append(-1)
        RC.append(-1)
        LS.append(i-1)
        RS.append(i+1)
        current.append(True)
    # if not a closed loop, set start/end left/right siblings to -1
    RS[n-1]=-1
    LS[0]=-1
    # return packed list
    return [ID,XY,error, parent, LC, RC, LS, RS, current]

def __bounding_box(p1,p2):
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
    return (left, bottom, right, top)

def __segment_bounding_box_generator(pts, other_pt_lists = None):
    """Returns a generator that yields tuples of (id, bounding box, None)
       for all of the segments between consecutive points in pts.
       This is used to create the R-Tree index."""
    segs=list(zip(pts[:-1],pts[1:]))
    for i, obj in enumerate(list(zip(range(len(segs)),segs))):
        # Get left, bottom, right and top coordinate values
        seg=obj[1]
        p1=seg[0]
        p2=seg[1]
        left, bottom, right, top=__bounding_box(p1,p2)
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
                left, bottom, right, top=__bounding_box(p1,p2)
                # Yield tuple as desired by rtree (sequence: LBRT)
                yield (cur_id,(left, bottom, right, top),None)
                cur_id +=1


def displacementArea(simp_table,n):
    """
    Extracts the last (largest) diplacement area associated with simplification to n points

    Parameters
    ----------
    simp_table : list
        derived from the simplificationTable function.
    n : integer
        number of points in simplified line
    
    Returns
    -------
    disp_area : float
        the maximum displacement area associated with any single collapse event.

    """
    disp_areas=simp_table[2]
    orig_n = int(int(len(disp_areas)+3)/2)
    index_of_last_added_point = orig_n + (orig_n-n-1)
    
    return disp_areas[index_of_last_added_point]

def simplificationTable(pts, p_func=__placement_AP_EAmin, m_func=__sideshift_area, do_topo_check = False, other_pt_lists = None):
    """Collapses segments and builds simplification tree.
    from which simplified versions of line can be constructed.
    Args:
       pts: a list of (x,y) tuples making up a line
       optional p_func: f(A,B,C,D,m_func) returning E, displacement_metric
       optional m_func: f(A,B,C,D,E) returning displacement_metric.
       optional do_topo_check: if true, segments will not be collapsed if the produce a topological error
       optional other_pt_lists: list of lists of (x,y) tuples of other features to be included in the topo check
    Returns:
       List of vertexRelation objects, each of which contains the following
       accessible properties:
        id, xy, error, parent, LC, RC, LS, RS, current"""
    # build vertex relations on initial points
    vtree = __initVertexTree(pts)
    ID,XY,E,parent,LC,RC,LS,RS,current=vtree[0],vtree[1],vtree[2],vtree[3],vtree[4],vtree[5],vtree[6],vtree[7],vtree[8]
    # build priority list of collapse events
    pl = __init_priority_list(vtree,p_func, m_func,do_topo_check)
    # build spatial index for topology check
    # segment indices are indices of first vertex
    if do_topo_check:
        # create index of segments; use stream loading
        idx=index.Index(__segment_bounding_box_generator(XY,other_pt_lists))
    # initialize values for loop
    E = len(pts)-1
    maxerror=0
    collapsed=0
    report_interval=2000
    while len(pl)>0:
        error,A,B,C,D,pE, overlap_endpt = pl.pop(0)
        if current[A] and current[B] and current[C] and current[D]:
            topo_error=False
            if do_topo_check:
                topo_error=__causesSelfIntersection(A,B,C,D,pE,overlap_endpt,LS[A],vtree,idx)
            if not topo_error:
                E+=1
                maxerror=max(error,maxerror)
                if do_topo_check:
                    # remove AB, BC & CD from index
                    idx.delete(A,__bounding_box(XY[A],XY[B]))
                    idx.delete(B,__bounding_box(XY[B],XY[C]))
                    idx.delete(C,__bounding_box(XY[C],XY[D]))
                # perform collapse
                __collapse(B,C,E,pE,maxerror,vtree)
                if do_topo_check:
                    # add AE and ED to index
                    idx.insert(A,__bounding_box(XY[A],XY[E]))
                    idx.insert(E,__bounding_box(XY[E],XY[D]))

                collapsed += 1
                if collapsed % report_interval==0:
                    print("collapsed {} of {}".format(collapsed, len(pts)))
                try:
                    __addToPriorityList(LS[LS[A]],LS[A],A,E,vtree, pl, p_func, m_func,do_topo_check)
                    __addToPriorityList(LS[A],A,E,D, vtree, pl, p_func, m_func, do_topo_check)
                    __addToPriorityList(A,E,D,RS[D], vtree, pl, p_func, m_func, do_topo_check)
                    __addToPriorityList(E,D,RS[D],RS[RS[D]], vtree, pl, p_func, m_func, do_topo_check)
                except:
                    print('Error at {}'.format(E))
    return vtree


def simplifiedLine(simp_table,max_disp=-1, min_pts=-1):
    """ Returns a list of (x,y) tuples representing vertices of simplified line.
    Args:
        simp_table: simplification table derived from simplificationTable function
        max_disp:   maximum allowable displacement at each step
        min_pts:    output must include at least this many points
    (Either neither max_disp nor min_pts are supplied...
    """
    # vtree lists:
    # id, xy, error, parent, LC, RC, LS, RS, current
    ID,XY,error,parent,LC,RC,LS,RS,currents=simp_table[0],simp_table[1],simp_table[2],simp_table[3],simp_table[4],simp_table[5],simp_table[6],simp_table[7],simp_table[8]
    # determine original number of points,
    # maximum ID to allow based on displacement error criterion
    max_ID=len(ID)-1
    orig_count=0
    for i in range(len(ID)):
        if LC[i] == -1:
            orig_count +=1
        if max_disp >= 0:
            if error[i] <= max_disp:
                  max_ID=i
    # determine max ID to allow based on minimum number of points
    if min_pts > -1:
        max_ID_by_min_pts = orig_count -1 + orig_count - min_pts
        if max_disp > -1:
            max_ID = max([max_ID,max_ID_by_min_pts])
        else:
            max_ID = max_ID_by_min_pts
    # initialize output as empty list
    out_line=[]
    # start at beginning
    v=0
    # work to end
    while v > -1:
        if v <= max_ID:  # vertex is good
            out_line.append(XY[v]) # add to result
            # if we can't move right, move up
            while v > -1 and RS[v] == -1:
                v = parent[v]
            # move right if we're not at the end already
            if v > -1:
                v = RS[v]
        else: # vertex needs to be expanded; move to left child
            v = LC[v]
    # that's it!
    return out_line

def save(tree,tofile):
    import pickle
    output=open(tofile,'wb')
    pickle.dump(tree,output,-1)
    output.close()

def load(treefile):
    import pickle
    pkl_file = open(treefile,'rb')
    tree=pickle.load(pkl_file)
    pkl_file.close()
    return tree

if __name__ == '__main__':
    pass