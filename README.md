# line-simplify
A collection of line simplification algorithms in pure python code

# Introduction
This repository contains basic implementations of several line simplification algorithms, for use in research and teaching. Currently the following algorithms are implemented:
- Area-Preserving Segment Collapse (APSC)
- the Ramer/Douglas-Peucker algorithm
- the Visvalingam-Whyatt algorithm
- Raposo's algorithm  

The Visvaling-Whyatt algorithm code is taken from Elliot Hallmark with only minor edits for consistency with other algorithms - the original project can be found at https://pypi.org/project/visvalingamwyatt/. All other algorithms were implemented by Barry Kronenfeld based on their descriptions in published papers. The algorithms run fairly quickly, but are not optimized to be blazingly fast. The Ramer/Douglas-Peucker algorithm in particular is implemented in a straightforward manner and does not utilize faster techniques that have been developed.

# Getting Started
1. Download the files
2. Install sortedcontainers (http://www.grantjenks.com/docs/sortedcontainers/)
3. Install Rtree (https://www.lfd.uci.edu/~gohlke/pythonlibs/#rtree)
   For windows users:
   a. download the 32- or 64-bit wheel:
        Rtree-0.8.3-cp27-cp27m-win32.whl
        Rtree-0.8.3-cp27-cp27m-win_amd64.whl
   b. open a command-prompt as an administrator
   c. use 'cd' to navigate to the python scripts folder
   d. uninstall any previous version, e.g. "pip uninstall rtree" 
   e. install from the wheel, e.g. "pip install Rtree-0.8.3-cp27-cp27m-win32.whl"
4. Run the 'example_usage.py' script for examples of how to use each line simplification algorithm.
   Run the 'example_visualize.py' script to view an example of the APSC algorithm in action.
   Run the 'example_visualize_tree.py' script to view the tree-structure of the APSC 
      simplification table

# References
Descriptions of each line simplification algorithm are given in the following papers:

Kronenfeld, Stanislawski, Buttenfield & Brockmeyer (2019) Simplification of polylines by segment collapse: minimizing areal displacement while preserving area. International Journal of Cartography, https://doi.org/10.1080/23729333.2019.1631535

Ramer (1972) An iterative procedure for the polygonal approximation of plane curves. Computer Grpahics and Image Processing, vol. 1, pp. 244-256.

Douglas and Peucker (1973) Algorithms for the reduction of the number of points required to represent a digitized line or its caricature. The Canadian Cartographer, vol. 10, no. 2, pp. 112-123.

Visvalingam and Whyatt (1993) Line generalisation by repeated elimination of points. The Cartographic Journal, vol. 30, no. 1, pp. 46-51.

Raposo (2013) Scale-specific automated line simplification by vertex clustering on a hexagonal tessellation. Cartography and Geographic Information Science, vol. 40, no. 5, pp. 427-443, http://dx.doi.org/10.1080/15230406.2013.803707
