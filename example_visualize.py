#----------------------------------------------------------------------------
# Name:        visualize_simplified_lines
# Purpose:     provide interactive plot to experiment with
#              (area-preserving) line simplification
# Author:      bjkronenfeld
# Created:     13/07/2017
# Copyright:   (c) bjkronenfeld 2017
# License:     MIT License
#----------------------------------------------------------------------------
# If using spyder, set options to use automatic graphics backend. See here:
#   https://stackoverflow.com/questions/23585126/how-do-i-get-interactive-plots-again-in-spyder-ipython-matplotlib
#----------------------------------------------------------------------------

import time
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import missouri_river as mo_riv
import apsc

# input parameters
check_topo=False
orig_line = mo_riv.missouri_river()

# perform simplification and time how long it takes
startTime=time.time()
simp_table = apsc.simplificationTable(orig_line,do_topo_check=check_topo)
constructionTime = time.time() - startTime

# show some metadata
print('Feature has {} vertices.'.format(len(orig_line)))
print("Construction Time: {} seconds".format(constructionTime))

# show original line
x,y=zip(*orig_line)
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.45, bottom=0.15)
t0 = 0
orig, = plt.plot(x,y,lw=2, color='black')
l, = plt.plot(x,y, lw=2, color='red')
maxx=max(x)
minx=min(x)
maxy=max(y)
miny=min(y)
dx=maxx-minx
dy=maxy-miny
plt.axis([float(a) for a in [minx-0.1*dx, maxx+0.1*dx, miny-0.1*dy,maxy+0.1*dy]])
plt.axis('equal')

# create slider
ax_slider = plt.axes([0.15, 0.03, 0.45, 0.03])
n=len(orig_line)

slider_vertex_count = Slider(ax_slider, 'Vertices:', 3, n, valinit=n,valstep=1,valfmt="%1.0f")

# create information box
def get_info_text(n_target, simple_pts):
    n_simple=len(simple_pts)
    info_text = "Simplified Line:\n   {} vertices\n   {}% of original".format(n_simple,round(100.0*n_simple/n,4))
    if n_simple != n_target:
        info_text += "\nTarget:\n   {} pts".format(n_target)
    info_text += "\nMax Displacement:\n   {:.1f}".format(apsc.displacementArea(simp_table, n_simple))
    return info_text

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
info_box = fig.text(0.02, 0.98, get_info_text(n,orig_line), transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

# update information box when new value is selected in slider
def update(val):
    n_target = int(slider_vertex_count.val)
    simple_pts=apsc.simplifiedLine(simp_table,-1,n_target)
    x,y=zip(*simple_pts)
    l.set_xdata(x)
    l.set_ydata(y)
    info_box.set_text(get_info_text(n_target,simple_pts))
    ax_slider.clear()
    slider_vertex_count.__init__(ax_slider, 'Vertices:', 3, n, valinit=val,valstep=1,valfmt="%1.0f")
    slider_vertex_count.on_changed(update)
    
slider_vertex_count.on_changed(update)

plt.show()