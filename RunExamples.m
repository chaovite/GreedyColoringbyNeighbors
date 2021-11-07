% =========================================================================
%   *-------------------------------------------------------------------*
%   |              ~~* Greedy Graph Coloring Examples *~~               |
%   |              Written by Tomas Zegard, Oct 24th,2010               |
%   |                revised by Chao Liang,  Sep 1st, 2021           |
%   *-------------------------------------------------------------------*
% =========================================================================
clc, clear;
%% Example 1: Structured Q4 orthogonal mesh
close all
Nx=2; Ny=2;
disp('Coloring for a regular Quad mesh:')
[NODE,ELEM]=SimpleMesher('regquad',Nx,Ny);

% Do the coloring once more
[C,ne,NumberOfColors]=GreedyColoringbyNeighbors(ELEM',NODE',5);
%% Example 2: Structured T3 orthogonal mesh
close all
Nx=5; Ny=4;
disp('Coloring for a regular Tria mesh:')
[NODE,ELEM]=SimpleMesher('regtria',Nx,Ny);
for threads=2:2:6
    [C,ne,NumberOfColors]=GreedyColoringbyNeighbors(ELEM',NODE',threads);
end

disp('Hit any key to continue...')
pause

%% Example 3: Circular mesh with T3
close all
R=1.0; Dmin=0.2;
disp('Coloring for a regular Tria mesh:')
[NODE,ELEM]=SimpleMesher('circle',R,Dmin);
for threads=8:8:16
    [C,ne,NumberOfColors]=GreedyColoringbyNeighbors(ELEM',NODE',threads);
end