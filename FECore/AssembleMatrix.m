function [ barK ] = AssembleMatrix( egnn, nne, tne, eK)
% This function assembles the FEM matrix
% 
% Author:           Abdullah Waseem       
% Created:          20-August-2017       
% Last Modified:    17-December-2017
% Contact: engineerabdullah@ymail.com
%
% References: 
% https://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/#view_comments
% RCV methos is used here.
%
% KEEPING THE ASSEMBLY SEPARATE BECAUSE SPARSE INDEXING IS SLOW.


% Length of the C,R and V columns   len = (#nodes in X) x (#nodes in Y) x tne
len = nne * nne * tne;
% Initializing the location and storage vectors.
Ridx = zeros(len,1);       
Cidx = zeros(len,1);      
VK = zeros(len, 1);
% Indexer
idx = 1;
% Loop over all the elements
for en = 1 : tne
    % Calling the global nodes in that element
    gn = egnn(en, :);
    % Loop over the element nodes
    for I = 1 : nne
        for J = 1 : nne
            % Storing
            Ridx(idx,1)   = gn(I);
            Cidx(idx,1)   = gn(J);
            VK(idx,1)  = eK(I,J,en);
            idx = idx + 1;
        end
    end
end
% Sparse K.
barK = sparse(Ridx, Cidx, VK);
    
end
