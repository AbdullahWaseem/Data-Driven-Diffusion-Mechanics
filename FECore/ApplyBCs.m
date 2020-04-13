%% Applying boundary conditions
%
% Author: Abdullah Waseem       
% Created: 19-August-2017
% Modified:17-January-2020
% Contact: engineerabdullah@ymail.com

% Global node iterator of the boundary segments
bottom_n_itrt = unique(bottom_egnn, 'stable');
right_n_itrt  = unique(right_egnn,  'stable');
top_n_itrt    = unique(top_egnn,    'stable');
left_n_itrt   = unique(left_egnn,   'stable');
   
% Prescribed nodes
p = unique([left_n_itrt; 
            right_n_itrt; 
            top_n_itrt; 
            bottom_n_itrt]);

% Free nodes
f = setdiff(1:tnn, p); 
%Converting it to column.
f = f';
