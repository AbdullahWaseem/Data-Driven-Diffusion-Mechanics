%% Applying boundary conditions
%
% Author: Abdullah Waseem       
% Created: 19-August-2017
% Modified:17-January-2020
% Contact: engineerabdullah@ymail.com

% Global node iterator of the boundary segments
bottom_n_itrt   = unique(bottom_egnn, 'stable');
right_n_itrt    = unique(right_egnn,  'stable');
top_n_itrt      = unique(top_egnn,    'stable');
left_n_itrt     = unique(left_egnn,   'stable');
   
% Prescribed nodes
p = unique([left_n_itrt; 
            right_n_itrt; 
            top_n_itrt; 
            bottom_n_itrt]);

% Free nodes
f = setdiff(1:tnn, p); 
f = f'; %Converting it to column.

% Matrix and Vector partitions according to boundary conditions

Kpp = barK(p,p); 
Kpf = barK(p,f); 
Kfp = barK(f,p); 
Kff = barK(f,f);

Fg_p_str = barFg_str(p,1); 
Fg_f_str = barFg_str(f,1);

Fj_p_str = barFj_str(p,1); 
Fj_f_str = barFj_str(f,1);