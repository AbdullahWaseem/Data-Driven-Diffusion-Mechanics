%% Reading .msh file
%
% Author: Abdullah Waseem       
% Created: 19-August-2017       
% Contact: engineerabdullah@ymail.com


% This mesh reader is provided with Gmsh software.
msh = GmshReader([filename '.msh']);

% Extracting mesh data from the output. 
dim = 2;                        % Dimensionality of the Problem

%% Node Data
tnn   = msh.nbNod;                % Total number of nodes
nitrt = 1 : tnn;                % Node iterator
xy    = msh.POS(:, 1 : dim);       % Nodal coordinates
xy    = xy * sc;                   % Scale changing
V     = sc*sc;

%% Element Data
elementnumber = max(msh.ELE_INFOS(end,2));                              % Element number according to GMSH documentation.
eitrt         = find(msh.ELE_INFOS(:, 2) == elementnumber);             % Element Iterator
tne           = size(eitrt,1);                                          % Total number of elements
nne           = size(find(msh.ELE_NODES(eitrt(1),:)>0),2);              % Number of nodes in an elements.
nnbe          = size(find(msh.ELE_NODES(1,:)>0),2);                     % Number of nodes in an boundary elements.
egnn          = msh.ELE_NODES(eitrt, 1 : nne);                          % Element global node numbering
lnn           = 1 : nne;                                                % Local node numbering
ph_tag        = msh.ELE_TAGS(eitrt,1,1);                                % Provides the physical tags of the elements.

%% Boundary Data
% In gmsh when creating geometry file the physical tags are made for the material distribution and
% applying boundary conditions. 
% The physical tags for the boundary in the inclusion.geo file are BASE=11, RIGHT=12, TOP = 13, LEFT = 14
% Read gmshreader.m file to collect the nodes of these tags for applying boundary conditions.

% Boundary element iterator to access information from variable msh.ELE_NODES
bottom_e_itrt = find(msh.ELE_TAGS(:,1) == 11);
right_e_itrt  = find(msh.ELE_TAGS(:,1) == 12);
top_e_itrt    = find(msh.ELE_TAGS(:,1) == 13);
left_e_itrt   = find(msh.ELE_TAGS(:,1) == 14);
% Global node numbering for the boundary elements using boundary element iterators.
bottom_egnn   = msh.ELE_NODES(bottom_e_itrt, 1:nnbe);
right_egnn    = msh.ELE_NODES(right_e_itrt,  1:nnbe);
top_egnn      = msh.ELE_NODES(top_e_itrt,    1:nnbe);
left_egnn     = msh.ELE_NODES(left_e_itrt,   1:nnbe);

% We will not use structure msh again so clearing it 
clear msh

%% 
disp(' ')
disp(['Number of Elements : ' num2str(tne)]);
disp(['Number of Nodes    : ' num2str(tnn)]);
disp(' ')