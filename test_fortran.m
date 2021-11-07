% test fortran
% this script is used to test the fortran program give exactly the same
% results as the matlab program.
%
% Must compile the fortran code first:
%
% type 'gfortran fortran/GreedyColoringbyNeighbors.f90 -o color_exe' in the terminal
%

%% Example 1: Structured Q4 orthogonal mesh
close all
Nx=5; Ny=4;

disp('Coloring for a regular Quad mesh:')
[NODE,ELEM]=SimpleMesher('regquad',Nx,Ny);

% write ELEM file for color_exe
ELEM_FILE = 'ELEM_regquad.dat';
COLOR_FILE = 'COLOR_regquad.dat';
NEB_FILE  = 'NEB_regquad.dat';

writeELEM(ELEM', ELEM_FILE);
nthreads = 4;

% Do the coloring using the matlab script
[C, ne, NumberOfColors]=GreedyColoringbyNeighbors(ELEM', NODE', nthreads);
% Do the coloring using the fortran script
cmd = sprintf('./color_exe  %s %d %s %s > log.txt',  ELEM_FILE, nthreads, COLOR_FILE, NEB_FILE);
system(cmd);

% read COLOR_FILE, and NEB_FILE and compare.
C_Fortran = importdata(COLOR_FILE);
ne_Fortran = importdata(NEB_FILE);
assert(norm(C_Fortran-C')==0 & norm(ne_Fortran-ne')==0);

fprintf('\n\n');
fprintf('Example 1, pass test!!\n');
fprintf('\n\n');

%% Example 2: Structured T3 orthogonal mesh
Nx=5; Ny=4;
disp('Coloring for a regular Tria mesh:')
[NODE,ELEM]=SimpleMesher('regtria',Nx,Ny);

% write ELEM file for color_exe
ELEM_FILE = 'ELEM_regtria.dat';
COLOR_FILE = 'COLOR_regtria.dat';
NEB_FILE  = 'NEB_regtria.dat';

writeELEM(ELEM', ELEM_FILE);
for nthreads=2:2:6
    [C, ne, NumberOfColors]=GreedyColoringbyNeighbors(ELEM',NODE', nthreads);
    
    % Do the coloring using the fortran script
    cmd = sprintf('./color_exe  %s %d %s %s > log.txt',  ELEM_FILE, nthreads, COLOR_FILE, NEB_FILE);
    system(cmd);

    % read COLOR_FILE, and NEB_FILE and compare.
    C_Fortran = importdata(COLOR_FILE);
    ne_Fortran = importdata(NEB_FILE);
    assert(norm(C_Fortran-C')==0 & norm(ne_Fortran-ne')==0);
end
disp('Hit any key to continue...')
fprintf('\n\n');
fprintf('Example 2, pass test!!\n');
fprintf('\n\n');

%%
R=1.0; Dmin=0.2;
disp('Coloring for a regular Tria mesh:');
[NODE,ELEM]=SimpleMesher('circle',R, Dmin);

% write ELEM file for color_exe
ELEM_FILE = 'ELEM_circletria.dat';
COLOR_FILE = 'COLOR_circletria.dat';
NEB_FILE  = 'NEB_circletria.dat';

writeELEM(ELEM', ELEM_FILE);

for nthreads=8:8:16
    [C,ne,NumberOfColors]=GreedyColoringbyNeighbors(ELEM',NODE',nthreads);
    
    % Do the coloring using the fortran script
    cmd = sprintf('./color_exe  %s %d %s %s > log.txt',  ELEM_FILE, nthreads, COLOR_FILE, NEB_FILE);
    system(cmd);

    % read COLOR_FILE, and NEB_FILE and compare.
    C_Fortran = importdata(COLOR_FILE);
    ne_Fortran = importdata(NEB_FILE);
    assert(norm(C_Fortran-C')==0 & norm(ne_Fortran-ne')==0);
end

fprintf('\n\n');
fprintf('Example 3, pass test!!\n');
fprintf('\n\n');

system('rm *.dat log.txt');