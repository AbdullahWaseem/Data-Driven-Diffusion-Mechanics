%% STEADY-STATE DATA DRIVEN SOLVER FOR DIFFUSION TYPE PROBLEMS.
%  ------------
% Author: Abdullah Waseem       
% Created: 16-Jan-2020
% Contact: engineerabdullah@ymail.com
% 
% -----------------------
% I N T R O D U C T I O N
% -----------------------
%
% Equation: $\nabla\cdot j =0$ where $j=-M\cdot\nabla\mu$ is assumed.

%% Initializing Code

clear; 
clc; 
close all;
format long;

path(pathdef);
addpath FECore/ PreProcessing/ PreProcessing/MeshFiles/

%% Pre-Processing

% .msh file in PreProcessing/MeshFiles/
filename = 'Square';

% Scale .geo file lengths
sc = 1;
run('MeshData');

%% Initializing Finite Element

% QUAD: 1,4,9 -- TRI: 3,6,12
ngp = 3;
run('GaussQuadrature');

run('ShapeFunctions');

%% Initializing DATA-DRIVEN.

% The reference mobility tensor (isotropy assumed).
Mref = 10 * eye(dim);

% Number of data points
nDP  = 10;

% Maximum Iteration.
maxIter = 50;

% Tolerance and distance value. 
tol     = eps;
prvDis  = 1;

% The set of conjugate quantities: s in S
% -----------------------------------------------
% * * *   I S O T R O P I C    M O D E L    * * *
% -----------------------------------------------

% Data-set S (this data should comply with the boundary conditions)
gradMu_prm = linspace(-0.1,0,nDP+1);
% gradMu_prm(find(gradMu_prm==0))=gradMu_prm(find(gradMu_prm==0)-1);
  tenJ_prm = linspace(Mref(1,1)*abs(min(gradMu_prm)),-Mref(1,1)*abs(max(gradMu_prm)),nDP+1);

% Noise in the data
tenJ_prm = tenJ_prm + (rand(1,length(tenJ_prm)).*tenJ_prm*1e-1);

% Plot the data
figure(1); plot(gradMu_prm, tenJ_prm, '.', 'MarkerSize', 17);
xlabel('g'), ylabel('j')

% s \in S at Gauss points
% -----------------------
% Randomly assiging s\inS
    % Matrix of random integers.
    RanMat = randi(nDP+1, [tne,ngp]);
    % Initializing
    gradMu_str = zeros(tne,ngp);
    tenJ_str  = gradMu_str;
    % Assigning data randomnly
    gradMu_str(:,:) = gradMu_prm(RanMat(:,:));
      tenJ_str(:,:) =   tenJ_prm(RanMat(:,:));

%% Finite Element Calculations

run('FECore');
run('FE_Mat_Vec');

%% Apply the Boundary Conditions

% Global node iterator of the boundary segments
bottom_n_itrt = unique(bottom_egnn, 'stable');
right_n_itrt  = unique(right_egnn,  'stable');
top_n_itrt    = unique(top_egnn,    'stable');
left_n_itrt   = unique(left_egnn,   'stable');

%    Mu   
% -------
% Prescribed nodes
pm = unique([left_n_itrt; 
            right_n_itrt]);

% Free nodes
fm = setdiff(1:tnn, pm); 
%Converting it to column.
fm = fm';

%    Lam    
% ---------
% Prescribed nodes
pl = pm;
% pl = unique([left_n_itrt; 
%             right_n_itrt; 
%             top_n_itrt; 
%             bottom_n_itrt]);

% Free nodes
fl = setdiff(1:tnn, pl); 
%Converting it to column.
fl = fl';

% P A R T I T I O N N I N G 
% ----->

% Stiffness Matrix for Mu
Kmmpp = barK(pm,pm); Kmmpf = barK(pm,fm); Kmmfp = barK(fm,pm); Kmmff = barK(fm,fm);
% Stiffness Matrix for Lam
Kllpp = barK(pl,pl); Kllpf = barK(pl,fl); Kllfp = barK(fl,pl); Kllff = barK(fl,fl);
% Internal Flux Vemtor for Mu
Fg_p_str = barFg_str(pm,1); Fg_f_str = barFg_str(fm,1);
% Internal Flux Vemtor for Lam
Fj_p_str = barFj_str(pl,1); Fj_f_str = barFj_str(fl,1);
% External Flux Vemtor (appear in Lam problem only).
rF = zeros(tnn,1); rFp = rF(pl,1); rFf = rF(fl,1);
% <-----

Kmmff = decomposition(Kmmff);
Kllff = decomposition(Kllff);

% Initialize Mu and Lam
Mu  = zeros(tnn,1);
Lam = Mu;

% Prescribing the Drichlet part of the boundary
% (be consistant with the data)
Mu(left_n_itrt,1)  = 1;
Mu(right_n_itrt,1) = 0;

%% Steady State Data-Driven Solver 

k = 1;
% Main data driven loop
while k < maxIter
    
    % Solve for Mu
    Mu (fm,1) = Kmmff \ ( Fg_f_str - Kmmfp*Mu(pm,1) );
 
    % Solve for Lam
    Lam(fl,1) = Kllff \ ( rFf - Fj_f_str - Kllfp*Lam(pl,1) );
    
    % Initializing
    Dis = 0;
    barFg_str = zeros(tnn,1);
    barFj_str = zeros(tnn,1);
    
    % Element loop.
    for en = 1 : tne
        
        % Calling the global node numbering
        gnn = egnn(en,:);
        
        % Integration loop
        for gs = 1 : ngp
            
            % Calculating local states at (k) 
            % -------------------------------
            gradMu =                              Be(:,:,en,gs)' * Mu (gnn,1);
            tenJ   = [tenJ_str(en,gs) 0]' + Mref*(Be(:,:,en,gs)' * Lam(gnn,1));
            
            % FIND ROTATIONS
            % --------------
            % Numerical
            [rot,krot] = findrotationNR(gradMu,[gradMu_str(en,gs) 0]',tenJ,[tenJ_str(en,gs) 0]',Mref);
%             % Analytical
%             rot = -atan( (gradMu_str(1)*gradMu(2)*Mref(1,1) + tenJ_str(1)*tenJ(2)*1/Mref(1,1)) / ...
%                          (gradMu_str(1)*gradMu(1)*Mref(1,1) + tenJ_str(1)*tenJ(1)*1/Mref(1,1)) );
            % Rotation Matrix
            Q = [ cos(rot) -sin(rot); sin(rot)  cos(rot)];
            
            % DISTANCE
            Pi = dot((gradMu-Q'*[gradMu_prm; zeros(1,length(gradMu_prm))])'*Mref,(gradMu-Q'*[gradMu_prm; zeros(1,length(gradMu_prm))])',2) + ...
                 dot((  tenJ-Q'*[  tenJ_prm; zeros(1,length(  tenJ_prm))])'*Mref,(  tenJ-Q'*[  tenJ_prm; zeros(1,length(  tenJ_prm))])',2);
            
            % Finding s\inS which minimizes distance.
            [gpDis, indx] = min(Pi);
            
            % Integrating global Distances
            Dis = Dis + gpDis*glw(gs)*detJcbn(en,gs);
            
            % Assiging new gradMu_str and tenJ_str
            gradMu_str(en,gs) = gradMu_prm(indx);
              tenJ_str(en,gs) =   tenJ_prm(indx);
            % <-----
            
            % Force Vector for Mu  problem
            barFg_str(gnn,1) = barFg_str(gnn,1) + Be(:,:,en,gs)*Mref*[gradMu_str(en,gs) 0]'*glw(gs)*detJcbn(en,gs);
            
            % Force Vector for Lam problem
            barFj_str(gnn,1) = barFj_str(gnn,1) + Be(:,:,en,gs)*      [tenJ_str(en,gs) 0]'*glw(gs)*detJcbn(en,gs);
            
        end
    end
    
    % Partitioning the Internal Flux Vectors Again
    Fg_f_str = barFg_str(fm,1);
    Fj_f_str = barFj_str(fl,1);
    
    % Post-processing
    disp(['Iteration: ' num2str(k) '    Tolerance: ' num2str(abs(Dis-prvDis))]);
    
    % Check convergence
    if abs(Dis-prvDis) < tol
        break;
    else
        k = k + 1;
        prvDis = Dis;
    end

end

%% Post-Processing

% Reference Solution Finite Element
% ---------------------------------
Mu_fe = zeros(tnn,1);
Mu_fe(pm,1) = Mu(pm,1);
Mu_fe(fm,1) = barK(fm,fm)\(-barK(fm,pm)*Mu(pm,1));

% The Gradient Calculations
% to check the consistency between the (made up) data and applied boundary
% conditions.
gradMu_fe = zeros(dim,tne,ngp);
for en = 1 : tne
    gnn = egnn(en,:);
    for gs = 1 : ngp 
        gradMu_fe(:,en,gs) = Be(:,:,en,gs)' * Mu_fe(gnn,1);
    end
end
%

% Plotting Lam at each iteration.
figure(2);
subplot(4,1,1);
pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Lam,'Colormap','jet','Mesh','on');
axis equal; axis tight;
title('\lambda');
drawnow;

% Plotting Mu at convergence
figure(2)
subplot(4,1,2);
pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Mu,'Colormap','jet','Mesh','off');
axis equal; axis tight;
title('\mu_{dd}');
drawnow;

% Plotting FE Solution
figure(2)
subplot(4,1,3);
pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Mu_fe,'Colormap','jet','Mesh','off');
axis equal; axis tight;
title('\mu_{fe}');
drawnow;

% Plotting the difference
figure(2)
subplot(4,1,4);
pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',(Mu_fe-Mu),'Colormap','jet','Mesh','off');
axis equal; axis tight;
title('\mu_{fe}-\mu_{dd}');
drawnow;
 
