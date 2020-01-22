%% Precalculation and storage of element data. 
%
% Author: Abdullah Waseem       
% Created: 19-August-2017       
% Contact: engineerabdullah@ymail.com
% 
% In a normal FE calculations it should be performed at each gauss point
% here it is done to gain speeds during post-processing.

% Initializing
Jcbn    = zeros(dim,dim,tne,ngp);
Be      = zeros(nne,dim,tne,ngp);
detJcbn = zeros(tne,ngp);

% Element loop.
for en = 1 : tne
    % Integration loop
	for gs = 1 : ngp

        % Jacobian Matrix
		Jcbn(:,:,en,gs) = xy(egnn(en,:),:)' * DPHI(:,:,gs);
		
        % Its determinant "the differential area"
		detJcbn(en,gs)  = det(Jcbn(:,:,en,gs))/Ad;
		
        % Strain Gradient Be				  % This is the inv(Jcbn) 
		Be(:,:,en,gs)   = DPHI(:,:,gs) * (Jcbn(:,:,en,gs)\eye(dim,dim));
        
	end
end
