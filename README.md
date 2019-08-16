# line-simplify
A collection of line simplification algorithms in pure python code

# Introduction
This repository contains basic implementations of several line simplification algorithms, including:
- Area-Preserving Segment Collapse (APSC)
- the Ramer/Douglas-Peucker algorithm
- the Visvalingam-Whyatt algorithm
- Raposo's algorithm  

The Visvaling-Whyatt algorithm was implemented by Elliot Hallmark - the original project can be found at https://pypi.org/project/visvalingamwyatt/. All other algorithms were implemented by Barry Kronenfeld based on their descriptions in published papers. The algorithms run fairly quickly, but are not optimized to be blazingly fast. The Ramer/Douglas-Peucker algorithm in particular is implemented in a straightforward manner and does not utilize faster techniques that have been developed.

# Getting Started
1. Download the files
2. Install sortedcontainers (http://www.grantjenks.com/docs/sortedcontainers/)
3. Run the `usage_examples.py` script for examples of how to use each module, or follow the instructions below

# Usage
Basic usage of each module is illustrated in the `usage_examples.py` file. All algorithms take a list of (x,y) tuples as input:

```
input_line = [(97.1,32.1),(138.6,200.1),(138.4,367.9),(177.8,535.8),(256.7,703.8),(337.8,869.2),(420.9,1031.8),(504.0,1197.0),(628.8,1365.0)]
```

To apply the APSC algorithm, first create a simplification table and then simplify to a given tolerance or target number of vertices:
```
import apsc
apsc_table = apsc.simplificationTable(input_line) 
apsc_line = apsc.simplifiedLine(apsc_table,500,-1) # simplify to 500 sq unit tolerance
apsc_line = apsc.simplifiedLine(apsc_table,-1,5) # simplify to 5 vertices

```

To apply the Ramer/Douglas-Peucker algorithm, first obtain error values and then simplify to a given number of points:
```
# Ramer -- Douglas-Peucker:
import douglas_peucker as dp
errors,sorted_errors = dp.get_errors_sortedErrors(input_line) 
rdp_line = dp.fast_simplify_by_numPts(input_line,25,errors,sorted_errors)  # simplify to 50 vertices

```

To apply the Visvalingam-Whyatt algorithm, first create a Simplifier object and then simplify to a given number of points or areal displacement threshold:
```
import visvalingam as vw 
vw_simplifier=vw.VWSimplifier(input_line) 
vw_line=vw_simplifier.from_number(25)
```

To apply the Raposo algorithm, specify a hexagon radius or else a target number of points:
```
import raposo
raposo_line = simplify_raposo(input_line,500) # simplify using 500 unit radius hexagons
raposo_line = raposo.simplify_to_n_points_raposo(25,input_line) # simply to n ~ 50 pts (sometimes will end up with n+1 or n-1)

```
More options are available for each function/module, refer to the documentation in the code.

# References
Descriptions of each line simplification algorithm are given in the following papers:

Kronenfeld, Stanislawski, Buttenfield & Brockmeyer (2019) Simplification of polylines by segment collapse: minimizing areal displacement while preserving area. International Journal of Cartography, https://doi.org/10.1080/23729333.2019.1631535

Ramer (1972) An iterative procedure for the polygonal approximation of plane curves. Computer Grpahics and Image Processing, vol. 1, pp. 244-256.

Douglas and Peucker (1973) Algorithms for the reduction of the number of points required to represent a digitized line or its caricature. The Canadian Cartographer, vol. 10, no. 2, pp. 112-123.

Visvalingam and Whyatt (1993) Line generalisation by repeated elimination of points. The Cartographic Journal, vol. 30, no. 1, pp. 46-51.

Raposo (2013) Scale-specific automated line simplification by vertex clustering on a hexagonal tessellation. Cartography and Geographic Information Science, vol. 40, no. 5, pp. 427-443, http://dx.doi.org/10.1080/15230406.2013.803707
