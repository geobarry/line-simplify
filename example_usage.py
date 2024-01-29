#-------------------------------------------------------------------------------
# line simplification usage examples
#-------------------------------------------------------------------------------

# input data should be a list of (x,y) tuples
pts=[(47.2,86.4),(63.9,103.3),(63.9,120.0),(63.9,136.8),(63.8,153.6),(79.6,170.4),(79.6,187.2),(96.3,202.9),(96.2,219.7),(112.9,236.5),(129.5,253.3),(145.3,270.2),(145.3,286.9),(162.0,302.6),(178.7,319.4),(194.5,336.3),(211.1,353.1),(227.8,369.9),(243.6,386.7),(243.5,403.5)]

# AREA-PRESERVING SEGMENT COLLAPSE (APSC)
import linesimplify
from linesimplify import apsc
simp_table=linesimplify.apsc.simplificationTable(pts) # create the simplification table 
simplifiedline = apsc.simplifiedLine(simp_table,min_pts=10) # simplify to 10 points
simplifiedline = apsc.simplifiedLine(simp_table,max_disp=50) # simplify to areal displacement tolerance


# OTHER ALGORITHMS
# Visvalingam & Whyatt (1993)
from linesimplify import visvalingam as v
import numpy as np
simplifier=v.VWSimplifier(np.array(pts)) # create simplifier object
simplified_line = simplifier.from_number(10).tolist() # simplify to 10 pts
simplified_line = simplifier.from_threshold(50).tolist() # simplify to areal displacement tolerance

# Ramer/Douglas Peucker (1972-3)
from linesimplify import douglas_peucker as dp
errors, sorted_errors = dp.get_errors_sortedErrors(pts) # get displacement distance for each point
simplified_line = dp.simplify_by_numPts(pts, 10, errors, sorted_errors) # simplify to 10 pts
simplified_line = dp.simplify_by_distance_tolerance(pts, errors, 10) # simplify to distance displacement tolerance

# Raposo (2013)
from linesimplify import raposo
simplified_line = raposo.simplify_raposo(pts, 5.0) # simplify by hexagon radius
simplified_line = raposo.simplify_to_n_points_raposo(pts,10) # simplify by target num pts

print("FINISHED.")