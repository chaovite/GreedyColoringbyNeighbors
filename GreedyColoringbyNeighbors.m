% =========================================================================
%   *-------------------------------------------------------------------*
%   |                  ~~* Greedy Graph Coloring *~~                   |
%   |              Written by Tomas Zegard, Oct 24th,2010        |
%   |                revised by Chao Liang, Sep 1st, 2021           |
%   *-------------------------------------------------------------------*
% =========================================================================
%
% Note the dimension of ELEM and NODE are the transpose of those in the
% original scripts by Tomas Zegard to be consistent with faster access in
% fortran arrays.
%
% Also, note that, for coloring purposes, only ELEM array is necessary.
% NODE is used for visualizing the colored mesh.
%
% ELEM : element connectivity [NNODE_PER_ELEM, NEL]
% NODE: nodal coordinates [DIM, NNODE]
% threads: number of threads.
% max_neb: maximum number of neigboring element that share a element or a node.
%
% The number max_neb is chosen to be large enough 
% but still substantially smaller than total number of nodes or element,
% default is 50.
%
% These helper functions avoid computing the communication matrix, which
% can be a waste of memory. Instead, I use a [max_neb, NEL] array to describe
% the neighboring elements for each element, drastically more efficient in
% memory storage and searching.
%


function [C,ne,NumberOfColors]=GreedyColoringbyNeighbors(ELEM, NODE, threads, max_neb)
if nargin<4
    max_neb = 50;
end
ne=FindNeigbors(ELEM, max_neb);
[C,NumberOfColors]=ColorbyNeighbors(ne,threads,max_neb);
PlotColoring(ELEM,NODE,C,NumberOfColors,threads)
return

%% Function to plot the colored mesh
function []=PlotColoring(ELEM,NODE,C,NumberOfColors,threads)
RGBColor=hsv(NumberOfColors);
EdgeColor=[0 0 0];
% EdgeColor='none';
figure, hold on, axis off, axis equal
for i=1:size(ELEM,2)
    XY=NODE(:, ELEM(:,i));
    fill(XY(1,:),XY(2,:),RGBColor(C(i),:),'EdgeColor',EdgeColor)
    text(mean(XY(1,:)),mean(XY(2,:)),num2str(C(i)),...
        'HorizontalAlignment','center')
end
title(['Coloring for ' num2str(length(C)) ' elements using ' ...
    num2str(NumberOfColors) ' colors and ' num2str(threads) ' threads'])
return
    
%% Function that colors the elements using the Greedy Algorithm
function [C,NumberOfColors]=ColorbyNeighbors(ne,threads,max_neb)

% mark the first element as color 1
nel = size(ne,2);
C=zeros(1, nel);

C(1)=1;
NumberOfColors=1;
ColorCount     = zeros(1, nel); 
ColorCount(1) = 1;

% loop through the 2nd to the last element
for i=2:nel
    BlockedColors  = ones(1, max_neb)*-1;
    Neighbors=ne(:, i);% find the neighbors of element i
    for j=1:length(Neighbors) % loop through neigbors
        if (Neighbors(j)<=0) 
            continue
        end
        if C(Neighbors(j))~=0 % if a neighbor has been colored with C(Neighbors(j))
            BlockedColors(j)=C(Neighbors(j));
            % add this color to blocked color.
        end
    end
    for j=1:NumberOfColors % loop through existing colors
        if ColorCount(j)~=threads % the count for color j has not reached number of threads
            IsFree=isempty(find(BlockedColors==j, 1));% if color j is not blocked
            if IsFree
                C(i)=j; % color element i with j
                ColorCount(j)=ColorCount(j)+1; %number of colors increases
                break;
            end
        end
    end
    
    %C(i) has not been colored, either colors are blocked or reach number
    %of threads
    %color C(i) by a new color
    if C(i)==0
        NumberOfColors=NumberOfColors+1;
        C(i)=NumberOfColors;
        ColorCount(NumberOfColors)=1;
    end
end
return

%% Function that computes the neighboring elements
function ne=FindNeigbors(ELEM,max_neb)

% ne: 50 * nel array
nel = size(ELEM,2);
nnode = max(max(ELEM));

node2el = zeros(max_neb, nnode); % each node can be owned by maximum 4 elements
cnt         = zeros(1, nnode);

for i = 1: nel
    nodes = ELEM(:, i);
    for j = 1: length(nodes)
        i_node = nodes(j);
        if isempty(find(node2el(:,i_node)==i, 1))
            cnt(i_node) = cnt(i_node) + 1;
            node2el(cnt(i_node),i_node) = i;
        end
    end
end

ne  = ones(max_neb, nel)*-1; %initialize to be -1
cnt2 = zeros(1, nel); % number of neighbors

for i = 1:nnode 
    n2e  = node2el(:, i);
    cnt_i = cnt(i); % number of elements 
    for k = 1: cnt_i
        e_i = n2e(k);
        for j = k+1: cnt_i
            e_j = n2e(j);
            % check if e_i, e_j has already been included
            if isempty(find(ne(:, e_i)==e_j, 1))
            % e_i and e_j are neighbors!
            cnt2(e_i) = cnt2(e_i) + 1;
            cnt2(e_j) = cnt2(e_j) + 1;
            ne(cnt2(e_i), e_i) = e_j;
            ne(cnt2(e_j), e_j) = e_i;
            end
        end
    end
end

return
