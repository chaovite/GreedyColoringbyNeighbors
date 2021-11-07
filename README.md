# GreedyColoringbyNeighbors
Greedy graph coloring a finite element mesh using element connectivity. This is an improved version of the original subroutines developed by Tomas Zegard.

These helper functions avoid computing the communication matrix, which can be a waste of memory. Instead, I use a [8, NEL] array to describe the neighboring elements for each element, drastically more efficient in memory storage and searching.
