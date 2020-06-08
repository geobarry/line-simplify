# to make this work:
# pip/conda install pydot (MIT License)
# if necessary, install graphviz from https://graphviz.gitlab.io/_pages/Download/Download_windows.html

import os, sys
try:
    import pydot
except ImportError:
    sys.exit("""You need pydot to run tree_grapher.
                Refer to https://pypi.org/project/pydot/ for more information.
                """)
# if os.path.exists('C:/Program Files (x86)/Graphviz2.38/bin/'):
#     os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
# else:
#     sys.exit("""You need GraphViz to run tree_grapher.
#     Refer to https://graphviz.gitlab.io/_pages/Download/Download_windows.html 
#     and make sure the path to the Graphviz bin folder is correct.
#                 """)

def __to_digits(x,numdigits):
    if x==0:
        return "0"
    else:
        y=int(x*10)
        s=str(y)
        s=s[:-1] + "." + s[-1]
        return s

def __getedge(parentid, childid, simp_table,include_error=True,numdigits=1):
    ''''''
    def label(node_id):
        explicit_id=simp_table[0][node_id]
        if include_error:
            node_error = __to_digits(simp_table[2][node_id])
            text= "{}\n{}".format(explicit_id,node_error)
        else:
            text= str(explicit_id)
        #text='<b>{}</b>'.format(text)
        return text
    return pydot.Edge(label(parentid),label(childid))


def createSimplificationTreeGraph(simp_table,outfile, include_error=True, numdigits=1):
    # Each row of vtree is a list of [id,xy,error, parent, LC, RC, LS, RS, current]
    # create treegraph object
    graph = pydot.Dot(graph_type='graph', dpi="300")
    graph.set_node_defaults(fontsize='24', penwidth='3', style="filled", fillcolor="0.000 0.000 0.900")
    graph.set_edge_defaults(penwidth='6')
    worklist=[]
    # add edges
    # tree root will be last row in tree table
    rootid = len(simp_table[0])-1
    rootleft = rootid
    rootright= rootid
    worklist.append(rootid)
    # work through list
    while len(worklist) > 0:
        curid = worklist.pop() # retrieve next vertex in tree
        curnode=[col[curid] for col in simp_table]
        if curid == rootleft: # add edge for left sibling of root
            rootleft=curnode[6]
            if rootleft > -1:
                graph.add_edge(__getedge(curid,rootleft,simp_table,include_error,numdigits))
                worklist.append(rootleft)
        if curnode[4] >= 0: # add edges for children
            graph.add_edge(__getedge(curid,curnode[4],simp_table,include_error,numdigits))
            graph.add_edge(__getedge(curid,curnode[5],simp_table,include_error,numdigits))
        if curid==rootright: # add edge for right sibling of root
            rootright=curnode[7]
            if rootright > -1:
                graph.add_edge(__getedge(curid,rootright,simp_table,include_error,numdigits))
                worklist.append(rootright)
        if curnode[4] >= 0: # add children to worklist
            worklist.append(curnode[4])
            worklist.append(curnode[5])
    # write graph image to file
    if outfile[-4:] != ".png":
        outfile="{}.png".format(outfile)
    graph.write_png(outfile)
