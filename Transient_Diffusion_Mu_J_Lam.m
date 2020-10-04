%% STEADY-STATE DATA DRIVEN SOLVER FOR DIFFUSION TYPE PROBLEMS.
%  ------------
% Author: Abdullah Waseem       
% Created: 29-Jan-2020
% Contact: engineerabdullah@ymail.com
% 
% -----------------------
% I N T R O D U C T I O N
% -----------------------
%
% Equation: $\dot{\mu} + \nabla\cdot j =0$ 
% $j=-M\cdot\nabla\mu$ and $\mu = \Lambda c$
% $\mathcal{D}=\Lambda M$
% j: flux vector, M: mobility tensor, mu: chemical potential, \Lambda: chemical modulus, c: concentration

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

% TEMPORAL DISCRETIZATION
%
% Total Time
T    = 10;
% Number of time steps
nTs  = 100;
% Time Step Size
dT   = T/nTs;
% Time iterator
tItrt= 0:dT:T;

%% Initializing Finite Element

% QUAD: 1,4,9 -- TRI: 3,6,12
ngp = 3;
run('GaussQuadrature');

run('ShapeFunctions');

%% Initializing DATA-DRIVEN.

% The reference diffusivity.
Dref = 10 * eye(dim);

% Number of data points
nDP  = 101;

% Maximum Iteration.
maxIter = 100;
% Tolerance and distance value. 
tol     = 1e-10;
prvDis  = 1;

% The set of conjugate quantities s in S
% -----------------------------------------------
% * * *   I S O T R O P I C    M O D E L    * * *
% -----------------------------------------------

% Data-set S (this data should comply with the boundary conditions)
gradMu_prm = linspace(-1,0.0001,nDP+1);
% gradMu_prm(find(gradMu_prm==0))=gradMu_prm(find(gradMu_prm==0)-1);
  tenJ_prm = linspace(Dref(1,1)*abs(min(gradMu_prm)),-Dref(1,1)*abs(max(gradMu_prm)),nDP+1);

% Noisy in the data
tenJ_prm = tenJ_prm + (rand(1,length(tenJ_prm)).*tenJ_prm*1e-1);

% Plot the Data
figure(1); plot(gradMu_prm, tenJ_prm, '.', 'MarkerSize', 17);
xlabel('g'), ylabel('j')
drawnow

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

%% Boundary Conditions

% Global node iterator of the boundary segments
bottom_n_itrt = unique(bottom_egnn, 'stable');
right_n_itrt  = unique(right_egnn,  'stable');
top_n_itrt    = unique(top_egnn,    'stable');
left_n_itrt   = unique(left_egnn,   'stable');
  
%    Mu   
% -------
% Prescribed nodes
pm = unique(left_n_itrt);

% Free nodes
fm = setdiff(1:tnn, pm); 
%Converting it to column.
fm = fm';

%     Lam     
% ------------
% Prescribed nodes
pl = pm;
% Free nodes
fl = setdiff(1:tnn, pl);
% Converting to column
fl = fl';

% P A T I T I O N I N G 
% ----->
% Stiffness Matrix for C
Kmmpp = barK(pm,pm); Kmmpf = barK(pm,fm); Kmmfp = barK(fm,pm); Kmmff = barK(fm,fm);
% Stiffness Matrix for Lam
Kllpp = barK(pl,pl); Kllpf = barK(pl,fl); Kllfp = barK(fl,pl); Kllff = barK(fl,fl);
% Mass "Coupling" for C,Lam
Mmlpp = barM(pm,pl); Mmlfp = barM(fm,pl); Mmlpf = barM(pm,fl); Mmlff = barM(fm,fl);
% Mamm "Coupling" for Lam,C
Mlmpp = barM(pl,pm); Mlmfp = barM(fl,pm); Mlmpf = barM(pl,fm); Mlmff = barM(fl,fm);
% Internal Flux Vector: C
Fg_p_str = barFg_str(pm,1); Fg_f_str = barFg_str(fm,1);
% Internal Flux Vector: Lam
Fj_p_str = barFj_str(pl,1); Fj_f_str = barFj_str(fl,1);
% Flux Vector
Fr = zeros(tnn,1); Frp = Fr(pl,1); Frf = Fr(fl,1);
% <-----

% Decomposition of the Coupled Matrix (it won't change during the solution).
CoupMat = decomposition([   Mlmff   -dT*Kllff;
                         dT*Kmmff       Mmlff]); 

% Initialize Mu and Lam
Mu  = zeros(tnn,length(tItrt));
Lam = Mu;
dMu = zeros(tnn,1);

% Prescribing Mu in time.
% ----->

% % Sinusoidal Loading
% osc  = 1;
% a    = 2*pi*osc*tItrt/T;
% Muosc = 1*sin(a);
% for n = 1 : length(tItrt)
%     Mu(left_n_itrt,n) = Muosc(n);
% end

% Constant Loading
Mu(left_n_itrt,:) = 1;

% <-----

% For the reference solution.
Mu_fe= Mu;

%% Transient Data Driven Solver

% Time Stepping (Euler Backward).
for n = 2 : length(tItrt)
    
    disp(['Time Step: ' num2str(n)])
    
    k = 1;
    % Main data driven loop
    while k < maxIter

        % Solve for Mu and Lam
        % ----->
        % RHS
        Fmf = dT*Fj_f_str - dT*Frf              + Mlmfp*(Mu(pm,n)-Mu(pm,n-1)) - dT*Kllfp*Lam(pl,n);
        Flf = dT*Fg_f_str - dT*Kmmfp*Mu(pm,n-1) - dT*Kmmff*Mu(fm,n-1)         - dT*Kmmfp*(Mu(pm,n)-Mu(pm,n-1)) - Mmlfp*Lam(pl,n);
        % SOVLE
        U = CoupMat \ [Fmf; Flf];
        % Assiging the solution
        dMu(fm,1) = U(1:length(fm));
         Mu(fm,n) = Mu(fm ,n-1) + dMu(fm,1);
        Lam(fl,n) = U(length(fm)+1:end);
        % <-----

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
                gradMu =                                 Be(:,:,en,gs)' * Mu (gnn,n) ;
                  tenJ = [tenJ_str(en,gs) 0]' + dT*Dref*(Be(:,:,en,gs)' * Lam(gnn,n));

                % FIND ROTATIONS
                % --------------
%                 % Numerical
%                 [rot,~] = findrotationNR(gradMu,[gradMu_str(en,gs) 0]',tenJ,[tenJ_str(en,gs) 0]',Dref);
                % Analytical (CAREFUL ABOUT THE ZERO IN THE DENOMINATOR)
                rot = -atan( (gradMu_str(en,gs)*gradMu(2)*Dref(1,1) + tenJ_str(en,gs)*tenJ(2)*1/Dref(1,1)) / ...
                             (gradMu_str(en,gs)*gradMu(1)*Dref(1,1) + tenJ_str(en,gs)*tenJ(1)*1/Dref(1,1)) );

                % Rotation Matrix
                Q = [ cos(rot) -sin(rot); sin(rot)  cos(rot)];

                % DISTANCE
                Pi = dot((gradMu-Q'*[gradMu_prm; zeros(1,length(gradMu_prm))])'* Dref        ,(gradMu-Q'*[gradMu_prm; zeros(1,length(gradMu_prm))])',2) + ...
                     dot((  tenJ-Q'*[  tenJ_prm; zeros(1,length(  tenJ_prm))])'*(Dref\eye(2)),(  tenJ-Q'*[  tenJ_prm; zeros(1,length(  tenJ_prm))])',2);

                % Finding s\inS which minimizes distance.
                [gpDis, indx] = min(Pi);

                % Integrating global Distances
                Dis = Dis + gpDis*glw(gs)*detJcbn(en,gs);

                % Assiging new gradMu_str and tenJ_str
                gradMu_str(en,gs) = gradMu_prm(indx);
                  tenJ_str(en,gs) =   tenJ_prm(indx);
                % <-----

                % Force Vector for Mu  problem
                barFg_str(gnn,1) = barFg_str(gnn,1) + Be(:,:,en,gs)*Dref*[gradMu_str(en,gs) 0]'*glw(gs)*detJcbn(en,gs);

                % Force Vector for Lam problem
                barFj_str(gnn,1) = barFj_str(gnn,1) + Be(:,:,en,gs)*     [  tenJ_str(en,gs) 0]'*glw(gs)*detJcbn(en,gs);

            end
        end
        
        % Partitioning the Internal Flux Vectors Again
        Fg_f_str = barFg_str(fm,1);
        Fj_f_str = barFj_str(fl,1);

        % Post-processing
        disp(['    Iteration: ' num2str(k) '    Tolerance: ' num2str(abs(Dis-prvDis)) ]);

        % Check convergence
        if abs(Dis-prvDis) < tol
            break;
        else
            k = k + 1;
            prvDis = Dis;
        end

    end
    
    disp(' ');
    
    % Plotting Lam at each k-th iteration
    figure(2);
    subplot(4,1,1);
    pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Lam(:,n),'Colormap','jet','Mesh','off');
    axis equal; axis off;
    title(['Iteration ' num2str(k) '   (\lambda)']);
    drawnow;
    
    % Plotting converged Mu at each time step n
    figure(2);
    subplot(4,1,2);
    pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Mu(:,n),'Colormap','jet','Mesh','off');
    axis equal; axis off;
    title('(\mu_{dd})');
    caxis([-1 1]);
    drawnow;
    
    % Finite Element Solution
    Mu_fe(fm,n) = (dT*Kmmff+barM(fm,fm)) \ (-dT*Kmmfp*Mu_fe(pm,n)-barM(fm,pm)*(Mu_fe(pm,n)-Mu_fe(pm,n-1)) + barM(fm,fm)*Mu_fe(fm,n-1));
    
    % Plotting Finite Element Solution 
    figure(2);
    subplot(4,1,3);
    pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',Mu_fe(:,n),'Colormap','jet','Mesh','off');
    axis equal; axis off;
    title('(\mu_{fe})');
    caxis([-1 1]);
    drawnow;
    
    % Difference in FE and DD Solvers
    figure(2);
    subplot(4,1,4);
    pdeplot(xy(:,1:2)',[],[egnn ones(size(egnn,1),1)]','XYData',(Mu_fe(:,n)-Mu(:,n)),'Colormap','jet','Mesh','off');
    axis equal; axis off;
    title('(\mu_{fe} - \mu_{dd})');
%     caxis([-1 1]);
    drawnow;
    
    % Plotting the solution on the bottom nodes (Its like 1D).
    figure(3);
    clf;
    plot(Mu(bottom_n_itrt,n));
    hold on;
    plot(Mu_fe(bottom_n_itrt,n));
    ylim([-1 1]);
    drawnow;
    
end
