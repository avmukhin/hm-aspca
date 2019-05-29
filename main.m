%% Clear variables, figures
clear variables
clc;
close all;
clf
set(0, 'defaultfigurewindowstyle', 'docked');
%% Adding path
addpath(strcat(pwd, '/NonLinDiff1D'));  % PATH to 1D non-linear diffusion solver
addpath(strcat(pwd, '/swap_pca'));      % PATH to adaptive strategies arXiv:1903.07220

swap_or_seq = 'swap';                   % Special parameter for testing extension_PCA vs. sequential basis extension                  

load IGEMS_1FMMUP1.mat                  % Dataset of prior model realizations
dataset = 0.01 * dataset;

D = dataset(:, 1);                      % Choosing real model
Dinit = dataset(:,20);                  % Choosing initial model
dataset(:,1) = [];                      % Deleting real model from dataset

L = pi;                                 % Spatial size of 1D problem
T = 10;                                 % Time
Nx = length(D)+1;                       % Number of cells in x
NxD = Nx-1;                             % Number of cells for diffusion coefficients (values at center => -1)
Nt = 51;                                % Number of timesteps
x = linspace(0.1, L+0.1, Nx)';          % Linspace for x
xD = linspace(0.1, L+0.1, NxD)';        % Linspace for x for diffusion coefficients

uinit = 2*ones(Nx, 1) + cos(0.03*(1:Nx))';  % Initial condition
uL = 0.01;                                  % Boundary condition of 2nd order on LEFT
uR = -0.02;                                 % Boundary condition of 2nd order on RIGHT
%%
par = PCA_par(dataset, 'svd', 0.7);         % Special class for parametrization (data, algorithm of decomposition, percent of energy to save)
%% Covariance matrix and its eigenvalues of dataset
figure('Name', 'Covariance matrix and its eigenvalues')
subplot(2,1,1)
imagesc(par.eigvector_full);
colorbar;
subplot(2,1,2);
semilogy(abs(diag(par.sigma_full)));
drawnow limitrate
%}
%%

eps = 1e-8;     % Convergence criterion for Newton-Raphson method
Kmin = 4;       % Lower coefficient for timestep adaptation for Newton-Raphson method
Kmax = 12;      % Higher coefficient for timestep adaptation for Newton-Raphson method
ml = 1.0;       % Multiplier for timestep adaptation for Newton-Raphson method

maxNewton = 100;    % Maxim number of iterations for Newton-Raphson method

solver_params = struct( ...
    'L', L, ...
    'T', T, ...
    'Nx', Nx, ...
    'NxD', NxD, ...
    'Nt', Nt, ...
    'uinit', uinit, ...
    'uL', uL, ...
    'uR', uR, ...
    'D', D, ...
    'eps', eps, ...
    'Kmin', Kmin, ...
    'Kmax', Kmax, ...
    'multiplier', ml, ...
    'maxNewton', maxNewton ...
    );
%% Solving Forward Problem for real model


fprintf('Building u dustribution...');
[u_true, timesteps3] = solve(solver_params, D);     % Saving resulting observations (u_true) and timesteps
disp('done');
fprintf('--q----\n');



figure('Name', 'Solution of non-linear diffusion problem');
hold on;
for i = 1:Nt
    plot(u_true(:,i), 'linewidth', 1 + (i==1));
end
drawnow limitrate;
%%

max_iter = 100;     % Maximum number of iterations for optimization algorithm
eps = 1e-6;         % Epsilon value for convergence criteria
useplot = false;    % To plot intermediate results or not

%% Optimization algorithms
% In this code, two algorithms for gradient optimization are presented
% - conjugate gradients method  \optimize_CG
% - quasi-newton method (SR1)   \optimize_qN

%% History Matching without parametrization
[a_no_cg, b_no_cg, ~] = optimize_qN(solver_params, u_true, Dinit, max_iter, eps, par, 'No', useplot, false, false); 
fprintf('NO is done\n')
fprintf('-----------------\n')

%% History Matching with simle PCA parametrization
[a_pca_cg, b_pca_cg, ~] = optimize_qN(solver_params, u_true, Dinit, max_iter, eps, par, 'PCA', useplot, false, false);  
fprintf('PCA is done\n')
fprintf('-----------------\n')

%% History Matching with Rotaton Strategy PCA (see preprint arXiv:1903.07220) 
[a_rot_cg, b_rot_cg, ex_rot_cg] = optimize_qN(solver_params, u_true, Dinit, max_iter, eps, par, 'PCA', useplot, 5, 'rotation');
fprintf('ROTATION is done\n')
fprintf('-----------------\n')

%% History Matching with Swap Strategy PCA (see preprint arXiv:1903.07220)
[a_swap_cg, b_swap_cg, ex_swap_cg] = optimize_qN(solver_params, u_true, Dinit, max_iter, eps, par, 'PCA', useplot, 5, 'swap');
fprintf('SWAP is done\n')
fprintf('-----------------\n')

%% History Matching with Extension Strategy PCA (see preprint arXiv:1903.07220)
[a_ext_cg, b_ext_cg, ex_ext_cg] = optimize_qN(solver_params, u_true, Dinit, max_iter, eps, par, 'PCA', useplot, 5, 'extend');
fprintf('EXT is done\n')
fprintf('-----------------\n')

results;