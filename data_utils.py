#-------------------------------------------------------------------------------
# Name:        line_intersection
# Purpose:     utilities for sorting and inserting in sequence
# Author:      Barry Kronenfeld
# Created:     Nov. 2019
# Licence:     MIT License
#-------------------------------------------------------------------------------

from operator import itemgetter


def shift_list(a_list,new_start):
    return a_list[new_start:] + a_list[:new_start]

def min_id(val_list):
    minid, minval = min(enumerate(val_list), key=itemgetter(1))
    return minid

def sequence_info(a_list):
    """
    Provides a reverse lookup for any sortable list.

    Parameters
    ----------
    a_list : list
        The items to be put in sequence. A list, or any indexable collection of sortable values.

    Returns
    -------
    seq : list of integers
        The indices of the original list, in ascending order of value.
    ranks : TYPE
        The ranks of each value in the original list.

    """
    seq = [(i,a_list[i]) for i in range(len(a_list))]
    seq.sort(key = lambda x:x[1])
    seq = [x[0] for x in seq]
    ranks = [-1]*len(seq)
    for i in range(len(seq)):
        ranks[seq[i]]=i
    return seq,ranks

def bulk_insert(orig, insertion_list, track_id=0, do_sort = True):
    """
    Duplicates a list with new values inserted. 
    
    Parameters
    __________
    orig :: list of values
      The original list to be inserted into
    inserts :: list of (id, value) tuples
      Each id should denote the place to insert the new value
    track_id :: integer
      return value will include the id of the same item after insertion
    do_sort ::  boolean
      If False, the inserts list should be properly sorted
    
    Returns
    _______
    list, integer
      a new list containing the original items and inserted items, and the
      id of the tracked item 
      
    """
    # copy insertion list
    inserts = [x for x in insertion_list]
    # sort
    if do_sort:
        inserts.sort(key=itemgetter(0))
    #initialize
    r=[]
    new_track_id = -1
    origID,insertID = 0,0
    n = len(orig)
    inserts.append((n,None)) # to avoid out of range error
    # loop
    while origID < n:
        while inserts[insertID][0] == origID:
            r.append(inserts[insertID][1])
            insertID +=1
        r.append(orig[origID])
        if origID == track_id:
            new_track_id = len(r)-1
        origID += 1
    return r, new_track_id
    

