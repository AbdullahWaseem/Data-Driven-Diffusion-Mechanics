function VtkWriter(filename, xyz, egnn, u, steps, elementtype)
%
% Author:   Abdullah WASEEM
% Email:    engineerabdullah@ymail.com
% Created : May 03, 2018
% Modified: May 30, 2018
% 
% VTK Writer for 2 Dimensional Elements. (can easily be extended to 3D).
% It will create a directory named "VTK_Output" and write *.vtk files in it. 
% 
% Only STRUCTURED meshes are allowed.
% Only scalar fields are allowed.
% 
% Nomenclature:
% 		tnn:			Total number of nodes in the mesh.
%		nne: 			Number of nodes in an element.
% 		tne:			Total number of elements.
%		tnts: 			Total number of time steps.
%		dim: 			Dimension of the problem.
%
% Inputs:    
%       filename: 		The name of the file/sequence of files.
%       xyz:    		The X, Y and Z coordinates (tnn x dim).
%       egnn:  	 		Element Global Node Numbering {Connectivity matrix} (tne x nne).
%       u:      		Solution vector (tnn x tnts).
%       steps:  		In transient problems it is the number of steps between two output files
%							i.e., 1:steps:tnts and steps=1 for steady-state problems.
%		elementtype: 	Type of the element, either TRI3, TRI6, QUAD4 or QUAD8 are allowed.
% 
% Reference: https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
% 
% ------------------------------------------------------------------------
% WARNING! IT WILL DELETE ALL THE EXISTING .vtk FILES COMMENT OUT LINE 59.
% ------------------------------------------------------------------------

% Checking for the elementtype
switch elementtype
    case 'TRI3'
        elnu = 5;
    case 'TRI6'
        elnu = 22;
    case 'QUAD4'
        elnu = 9;
    case 'QUAD8'
        elnu = 23;
end

nne = size(egnn,2);
tnn = size(u,1);
tne = size(egnn,1);
tnts = size(u,2);

% Check for the output folder, if its not present creates it.
thefolder = './PostProcessing/VTKOutput';
if exist(thefolder,'dir') == 0
	mkdir(thefolder);
end
cd(thefolder);
delete *.vtk

disp('  Writing VTK output files.');
for t = 1 : steps : tnts

    % Openning the file
    fileID = fopen([filename '_' num2str(t) '.vtk'], 'a');
    % Header
    fprintf(fileID, '# vtk DataFile Version 2.0\n');
    % Title
    fprintf(fileID, 'Cool Simulation Data\n');
    % Data Type
    fprintf(fileID, 'ASCII\n');
    % Geometry/Topology
    fprintf(fileID, 'DATASET UNSTRUCTURED_GRID\n');

    % ----- Starting Attributes -----

    % xyz coordinates
    fprintf(fileID, ['POINTS ' num2str(tnn) ' float\n']);
    fprintf(fileID, '%f %f %f\n', xyz');
    fprintf(fileID, '\n');

    % Connectivity Array
    fprintf(fileID, ['CELLS ' num2str(tne) ' ' num2str(size(egnn,1)*size(egnn,2)+size(egnn,1)) '\n']);
    vec = [nne*ones(tne,1) egnn-1];
    if elnu == 5
        fprintf(fileID, '%d %d %d %d\n', vec');
    elseif elnu == 9
        fprintf(fileID, '%d %d %d %d %d\n', vec');
    elseif elnu == 22
        fprintf(fileID, '%d %d %d %d %d %d %d\n', vec');
    elseif elnu == 23
        fprintf(fileID, '%d %d %d %d %d %d %d %d %d\n', vec');
    end
    fprintf(fileID, '\n');

    % Cell types
    fprintf(fileID, ['CELL_TYPES ' num2str(tne) '\n']);
    vec = elnu*ones(tne,1);
    fprintf(fileID, '%d\n', vec);
    fprintf(fileID, '\n');

    % Scalar Field
    fprintf(fileID, ['POINT_DATA ' num2str(tnn) '\n']);
    fprintf(fileID, 'SCALARS scalar1 float 1\n');
    fprintf(fileID, 'LOOKUP_TABLE default\n');
    vec = u(:,t);
    fprintf(fileID, '%f\n', vec);

    % Closing the file.
    fclose(fileID);

end

disp(' ');
disp('  Done writing files.');

cd ../..

end
