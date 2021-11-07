# GreedyColoringbyNeighbors
Greedy graph coloring a finite element mesh using element connectivity. This is an improved version of the original subroutines developed by Tomas Zegard, see http://paulino.ce.gatech.edu/educational_GreedyGraphCol.html.

These helper functions avoid computing the communication matrix, which can be a waste of memory. Instead, I use a [max_neb, NEL] array to describe the neighbors of each element, which is drastically more efficient in memory and also faster for searching. The number max_neb must be larger than the maximum number of neighboring elements that share the same element or the same node. I set max_neb = 50 by default.

I also add a fortran version of the program that takes in the connectivity array written in a text file and output the colors in text file. 

To compile:
`gfortran fortran/GreedyColoringbyNeighbors.f90 -o color_exe`

To run:
`./color_exe ELEM_FILE nthreads COLOR_FILE [max_neb] [NEB_FILE]`

Both `max_neb=50` and `BEB_FILE=''` are optional. 
