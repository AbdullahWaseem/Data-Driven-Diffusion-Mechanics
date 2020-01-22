%% Force Vectors for DD
%
% Author: Abdullah Waseem       
% Created: 16-January-2020
% Contact: engineerabdullah@ymail.com

% Initialization
barFg_str = zeros(tnn,1);
barFj_str = zeros(tnn,1);

Ke = zeros(nne, nne, tne);

% Element loop.
for en = 1 : tne
    
    % Calling the global node numbering
    gnn = egnn(en,:);
    
    % Integration loop
    for gs = 1 : ngp
        
		% Force Vector for Mu problem
        barFg_str(gnn,1) = barFg_str(gnn,1) + Be(:,:,en,gs)*Dref*gradMu_str(:,en,gs)*glw(gs)*detJcbn(en,gs);
        
        % Force Vector for Lam problem
        barFj_str(gnn,1) = barFj_str(gnn,1) + Be(:,:,en,gs)*       tenJ_str(:,en,gs)*glw(gs)*detJcbn(en,gs); 
		
        % The Stiffness Matrix
        Ke(:,:,en) = Ke(:,:,en) + Be(:,:,en,gs)*Dref*Be(:,:,en,gs)'* glw(gs)*detJcbn(en,gs);
    end
end

% Assembly
barK = AssembleMatrix(egnn, nne, tne, Ke);
