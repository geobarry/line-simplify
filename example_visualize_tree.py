#-------------------------------------------------------------------------------
# This script illustrates how to write the APSC simplification table to a csv 
# file, and how to visualize the simplification table as a tree-like graph
# using the tree_grapher.py module. The table/graph contains all original 
# vertices and new vertices, and links them in collapse events, so that the 
# parent of two vertices is the new vertex formed by their collapse. If the 
# original line has n vertices, the table will have 2n-3 rows numbered 0 to 
# 2n-4, and the graph will similarly have 2n-3 nodes numbered 0 to 2n-4.
#
# The tree-like structure shown in the output *.png file can be traversed  
# visually to reconstruct the line at any simplification level. Each version of 
# the line is obtained by traversing the graph from left to right. The original 
# line is obtained by traversing the leaf nodes (0,1,2,3,...,n-1), while 
# simplified versions are obtained by designating a cutoff id x (n <= x <= 2n-4). 
# Then traverse from left to right again, but at each node move up the tree 
# until the first ancestor vertex with id >= x is reached. The original line's first and 
# last vertex are always included.
#-------------------------------------------------------------------------------

# parameters
out_csv = r'C:\temp\APSC_table.csv'
out_png = r'C:\temp\APSC_table_tree.png'
check_topo = False

# imports
import apsc
import missouri_river
import tree_grapher as tg

# utility to write csv
def write_simplification_table_to_csv(simp_table,csvfile):
    # simp_table contains 9 lists, one each of (id,xy,error, parent, LC, RC, LS, RS, current)
    with open(csvfile,'w') as f:
        f.write(','.join('id,x,y,error,parent,left_child,right_child,left_sibling,right_sibling'.split(','))+'\n')
        for i in range(len(simp_table[0])):
            # get id & separate x & y coords
            csv_line = "{},{},{},".format(simp_table[0][i],simp_table[1][i][0],simp_table[1][i][1])
            # use loop to get remaining values, omitting 'current' value
            csv_line += ','.join([str(x[i]) for x in simp_table[2:-1]])+'\n' 
            # write line to csv file
            f.write(csv_line)

# get feature and simplify
pts=missouri_river.missouri_river()[:100] # limit to first 100 pts
simp_table = apsc.simplificationTable(pts,do_topo_check = check_topo)

# write simplification table to csv file
write_simplification_table_to_csv(simp_table,out_csv)

# create png image of tree structure of table
tg.createSimplificationTreeGraph(simp_table,out_png,False)

print("finished")