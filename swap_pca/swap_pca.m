function [x_prev, par] = swap_pca(x_prev, gradF, par, mode, n)
%% Andrei Mukhin, Alexey Khlyupin, 2019, CETX, CET MIPT
% Swap PCA is an algorithm that allows adaptive modification of PCA basis
% during parametrization process. Paper can be found at: arXiv:1903.07220




%% Input parameters:
%   x_prev;                   % previous value of parametrized vector. Changes only in "swap" mode
%   gradF;                    % gradient of objective function w.r.t. x_prev
%   par;                      % parametrization class that contains all information about previously used parameters
%   mode;                     % = 'rotate' - rotation of initial basis is performed
%                             % = 'swap' - swap is performed
%   n;                        % number of vectors that should be swaped by algorithm. Only for "swap" mode
% 

% % Class PAR parameters:
%   N;                        % number of state vectors in parametrized ensemble
%   method;                   % eig or svd
%   n_threshold;              % number of valuable eigenvectors
%   eigvector;                % matrix of n_threshold eigenvectors
%   sigma;                    % diagonal matrix of n_threshold eigenvalues
%   eigvector_full;           % full eigenvectors matrix
%   sigma_full;               % full eigenvalues matrix
%   ES;                       % eigvector*sigma
%   invES;                    % inv(sigma)*eigvector.'
%   mean_data;                % mean value of data
%   data;                     % input dataset
 
eigvector = par.eigvector;
eigvector_full = par.eigvector_full; % full eigenvectors matrix
sigma_full = par.sigma_full; % full eigenvalues matrix

%% Firstly, we want to get indexes of old basis and missed components in full basis
pca_index = [];
non_index = [];

if strcmp(mode,'rotate')
    for j=1:size(eigvector_full,2)
        if j <= size(eigvector,2)
            pca_index(end+1) = j;
        else
            non_index(end+1) = j;
        end
       
    end
elseif strcmp(mode,'swap') || strcmp(mode,'extend')
    for j=1:size(eigvector_full,2)
        non_index(end+1) = j;
        for i=1:size(eigvector,2)
            if isequal(eigvector(:,i), eigvector_full(:,j))
                pca_index(end+1) = j;
                non_index(non_index == j) = [];
            end   
        end
    end
end




%% Now start calculations of coefficients required by S-PCA
coef_pca = zeros(1, length(pca_index));
coef_non = zeros(1, length(non_index));
coef_c1 = zeros(length(pca_index), length(non_index));
phi_upd = zeros(size(eigvector));

for ii=1:length(pca_index)
   
    pca_id = pca_index(ii);
    
    if (size(gradF,1) < size(gradF,2))
        coef_pca(ii) = gradF * eigvector_full(:, pca_id);
    else
        coef_pca(ii) = eigvector_full(:, pca_id)' * gradF;
    end
     
    for jj=1:length(non_index)
       
        non_id = non_index(jj);
        if (size(gradF,1) < size(gradF,2))
            coef_non(jj) = gradF * eigvector_full(:, non_id);
        else
            coef_non(jj) = eigvector_full(:, non_id)' * gradF;
        end
        
        coef_c1(ii,jj) = (coef_pca(ii)*coef_non(jj)) * sigma_full(non_id,non_id) / (sigma_full(pca_id,pca_id) - sigma_full(non_id,non_id));
        
        phi_upd(:, ii) = phi_upd(:, ii) + coef_c1(ii,jj) * eigvector_full(:, non_id);
        
    end
    
end

%% Basis update accordingly to chosen mode
if strcmp(mode,'rotate')
    
    gamma = min(0.6, 0.6 / max(sum(coef_c1, 2)));
    par.eigvector = par.eigvector + gamma * phi_upd;
    par.ES = par.eigvector * par.sigma;
    par.invES = inv(par.sigma) * (par.eigvector)'; % inv(sigma)*eigvector.'
    fprintf('Basis has been rotated \n');
elseif strcmp(mode,'swap')
    
    norm_upd = sqrt(sum(phi_upd.^2, 1));
    [~, pca_sort_idxs] = sort(norm_upd);
    norm_cnon = sqrt(sum(coef_c1.^2, 1));
    [~, non_sort_idxs] = sort(norm_cnon, 'descend');
    
    par.eigvector = [eigvector(:,pca_sort_idxs(1:end-n)), eigvector_full(:, non_index(non_sort_idxs(1:n)))];
    sigma_idxs = [pca_index(pca_sort_idxs(1:end-n)), non_index(non_sort_idxs(1:n))];
    sigma_full = diag(sigma_full); sigma_full = sigma_full(sigma_idxs);
    par.sigma = diag(sigma_full);
    par.ES = par.eigvector * par.sigma; %+ 1e-5 * ones(size(par.ES));
    par.invES = inv(par.sigma) * (par.eigvector)'; % inv(sigma)*eigvector.'
    x_prev = x_prev(pca_sort_idxs);
    x_prev(end-n+1, end) = 0;
    fprintf('Basis has been swaped \n');
    
    disp(['sigma_idx ', num2str(sigma_idxs)]);
elseif strcmp(mode,'extend')
    norm_cnon = sqrt(sum(coef_c1.^2, 1));
    [~, non_sort_idxs] = sort(norm_cnon, 'descend');
    par.eigvector = [eigvector, eigvector_full(:, non_index(non_sort_idxs(1:n)))];
    sigma_idxs = [pca_index, non_index(non_sort_idxs(1:n))];
    sigma_full = diag(sigma_full); sigma_full = sigma_full(sigma_idxs);
    par.sigma = diag(sigma_full);
    par.ES = par.eigvector * par.sigma; %+ 1e-5 * ones(size(par.ES));
    par.invES = inv(par.sigma) * (par.eigvector)'; % inv(sigma)*eigvector.'
    x_prev(end+1) = 0;
    fprintf('Basis has been extended by \n');
    disp(['sigma_idx ', num2str(sigma_idxs)]);
end
    
end

