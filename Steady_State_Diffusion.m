%% STEADY-STATE DATA DRIVEN SOLVER FOR DIFFUSION TYPE PROBLEMS.
% 
% Author: Abdullah Waseem       
% Created: 16-Jan-2020
% Contact: engineerabdullah@ymail.com
% 
% -----------------------
% I N T R O D U C T I O N
% -----------------------
% something about tht code.

%% Initializing Code

clear; 
clc; 
close all;
format short;

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

% The reference diffusivity.
Dref    = 10 * eye(dim);

% Number of data points
nDP     = 50;

% Maximum Iteration.
maxIter = 200;
% Tolerance and displacement value. 
tol     = 1e-9;
Dis_prv = 1e3;

% The set of conjugate quantities s in S
% -----------------------------------------------
% * * *   I S O T R O P I C    M O D E L    * * *
% -----------------------------------------------

% Data-set S
gMu = linspace(-10,10,nDP+1);
tnJ = linspace(-10,10,nDP+1);

% Noisy data
tnJ = tnJ + (rand(1,length(tnJ)).*tnJ*1e0);

% Plot the Data
figure(1); plot(gMu, tnJ, 'o');

gradMu_prm = [gMu; gMu];
  tenJ_prm = [tnJ; tnJ];

% s \in S at Gauss points
% -----------------------

% Matrix of random integers.
RanMat = randi(nDP+1, [tne,ngp]); 

% Randomly assiging s\inS .
gradMu_str = zeros(dim,tne,ngp);
tenJ_str   = gradMu_str;

% Assigning data randomnly. 
gradMu_str(1,:,:) = gMu(RanMat(:,:));
gradMu_str(2,:,:) = gMu(RanMat(:,:));
  tenJ_str(1,:,:) = Dref(1,1)*tnJ(RanMat(:,:));
  tenJ_str(2,:,:) = Dref(2,2)*tnJ(RanMat(:,:));

%% FINITE ELEMENT CALCULATIONS

run('FECore');
run('FE_Mat_Vec');

%% Apply the Boundary Conditions

run('ApplyBCs');

% Initialize Mu
Mu = zeros(tnn,1);
% Prescribing the Drichlet part of the boundary.
Mu(left_n_itrt,1)  = 100;
Mu(right_n_itrt,1) = 40;

% Initialize reaction fluxes
Fr = zeros(tnn,1);

%% DATA-DRIVEN SOLVER

k = 1;
% Main data driven loop
while k < maxIter
    
    % Solve for Mu
    Mu(f,1) = Kff \ ( Fg_f_str - Kfp*Mu(p,1) );
    
    % Post-process reaction fluxes on the prescribed domain.
    % q.n   = -q_n
    Fr(p,1) = -(Kpp*Mu(p,1) + Kpf*Mu(f,1) - Fg_p_str);
 
    % Solve for Lam
    Lam = barK \ (barFj_str - Fr);
    
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
            gradMu =                           Be(:,:,en,gs)' * Mu (gnn,1);
            tenJ   = tenJ_str(:,en,gs) - Dref*(Be(:,:,en,gs)' * Lam(gnn,1));
            
            % Seaching for gradMu_str, tenJ_str
            % ---------------------------------
            
            % Calculating ditances from each s in S
            Pi = 0.5*dot(gradMu-gradMu_prm(:,:),Dref*(gradMu-gradMu_prm(:,:))) + ...
                    0.5*dot(tenJ-tenJ_prm(:,:),(Dref\eye(dim))*(tenJ-tenJ_prm(:,:)));
            
            % minimization problem
            [Dis_gp, indx] = min(Pi);
            
            % Integrating global Distance
            Dis = Dis + Dis_gp * glw(gs) * detJcbn(en,gs);
            
            % Assiging new gradMu_str and tenJ_str
            gradMu_str(:,en,gs) = gradMu_prm(:,indx);
              tenJ_str(:,en,gs) =   tenJ_prm(:,indx);
            
            % ---------------------------------
            
            % Force Vector for Mu  problem
            barFg_str(gnn,1) = barFg_str(gnn,1) + Be(:,:,en,gs)*Dref*gradMu_str(:,en,gs)*glw(gs)*detJcbn(en,gs);

            % Force Vector for Lam problem
            barFj_str(gnn,1) = barFj_str(gnn,1) + Be(:,:,en,gs)*       tenJ_str(:,en,gs)*glw(gs)*detJcbn(en,gs);

        end
    end
    
    % Post-processing
    disp(['Iteration: ' num2str(k)]);
    disp(['Tolerance: ' num2str(abs(Dis-Dis_prv)/abs(Dis_prv))]);
    disp(' ');
    figure(2);
    pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Mu,'Colormap','jet','Mesh','on');
    drawnow
    
    % Check convergence
    if abs(Dis-Dis_prv)/abs(Dis_prv) < tol
        break;
    else
        k = k + 1;
        Dis_prv = Dis;
    end

end