% =========================================================================
%   *-------------------------------------------------------------------*
%   |                    ~~* Mesh Generator *~~                         |
%   |              Written by Tomas Zegard, Oct 24th,2010               |
%   *-------------------------------------------------------------------*
% =========================================================================

function [NODE,ELEM]=SimpleMesher(type,param1,param2)

% Simple mesher for orthogonal meshes of quad, tria and circular with tria
% --For 'regquad' or 'regtria, param1 and param2 are the number of elements'
%   in X and Y respectively.
% --For 'circle', param1 is the radius, and param2 is the minimum distance
%   between the points.
switch type  
    case 'regquad'
        [NODE,ELEM]=QuadMesh(param1,param2);
    case 'regtria'
        [NODE,ELEM]=TriaMesh(param1,param2);
    case 'circle'
        [NODE,ELEM]=CircleMesh(param1,param2);
    otherwise
        disp(['ERROR: Mesh type must be either ''regquad'', ''regtria'''...
              ' or ''circle''.'])
        NODE=[];
        ELEM=[];
end

%% Structured mesher with quads function
function [NODE,ELEM]=QuadMesh(Nx,Ny)
        
ELEM=zeros(Nx*Ny,4);
for i=1:Nx
    for j=1:Ny
        nodes=[(j-1)*(Nx+1)+i (j-1)*(Nx+1)+i+1 j*(Nx+1)+i+1 j*(Nx+1)+i];
        ELEM((j-1)*Nx+i,:)=nodes;
    end
end

NODE=zeros((Nx+1)*(Ny+1),2);
for i=1:Nx+1
    for j=1:Ny+1
        NODE((j-1)*(Nx+1)+i,:)=[i-1 j-1];
    end
end

%% Structured mesher with tria function
function [NODE,ELEM]=TriaMesh(Nx,Ny)

ELEM=zeros(2*Nx*Ny,3);
for i=1:Nx
    for j=1:Ny
        nodes=[(j-1)*(Nx+1)+i (j-1)*(Nx+1)+i+1 j*(Nx+1)+i+1];
        ELEM((j-1)*Nx+i,:)=nodes;
        nodes=[(j-1)*(Nx+1)+i j*(Nx+1)+i+1 j*(Nx+1)+i];
        ELEM((j-1)*Nx+i+Nx*Ny,:)=nodes;
    end
end

NODE=zeros((Nx+1)*(Ny+1),2);
for i=1:Nx+1
    for j=1:Ny+1
        NODE((j-1)*(Nx+1)+i,:)=[i-1 j-1];
    end
end

%% Circular mesher with tria function
function [NODE,ELEM]=CircleMesh(R,Dmin)

% Seed the border
theta=linspace(0,2*pi,round(2*pi*R/Dmin));
NODE=R*[cos(theta(1:end-1))' sin(theta(1:end-1))'];

% Fill in the interior
[theta,rho]=meshgrid(linspace(0,2*pi,round(2*pi/Dmin)),...
    linspace(0,R-Dmin/3,round(R/Dmin)));
X=rho.*cos(theta); Y=rho.*sin(theta);
X=X(:,1:end-1); Y=Y(:,1:end-1);
Pcandidate=[reshape(X,size(X,1)*size(X,2),1) ...
    reshape(Y,size(Y,1)*size(Y,2),1)];

% Retain the nodes than meet the specified min distance parameter
for i=1:size(Pcandidate,1)
    TooClose=0;
    for j=1:size(NODE,1)
        if norm(Pcandidate(i,:)-NODE(j,:))<Dmin
            TooClose=1;
            break;
        end
    end
    if TooClose==0
        NODE=[NODE; Pcandidate(i,:)];
    end
end
% Use delaunay triangulation to generate the elements given the points
ELEM=delaunay(NODE(:,1),NODE(:,2));