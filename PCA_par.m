classdef PCA_par %< handle
    properties
        N; % number of state vectors in parametrized ensemble
        method; % eig or svd
        n_threshold; % number of valuable eigenvectors
        eigvector; % matrix of n_threshold eigenvectors
        sigma; % diagonal matrix of n_threshold eigenvalues
        eigvector_full; % full eigenvectors matrix
        sigma_full; % full eigenvalues matrix
        ES; % eigvector*sigma
        invES; % inv(sigma)*eigvector.'
        mean_data; % mean value of data
        data; % input dataset
    end
    methods
        function self = PCA_par(data, method, t)
            disp('Initializing PCA');
            
            self.data = data;
            self.N = size(data,2);
            if strcmp(method, 'eig') || strcmp(method, 'svd')
                self.method = method;
            else
                error(['Wrong decomposition method: ', method]);
            end
            
            self.mean_data = mean(self.data, 2);
            p = self.data - repmat(self.mean_data, 1, self.N);
            P =  p * p.' / (self.N - 1);
            
            if (t > 1 && t <= self.N)
                [self.eigvector_full, self.sigma_full, ~] = self.es(P, 1, self.method);
                self.n_threshold = t;
            elseif (t > 0 && t <= 1)
                [self.eigvector_full, self.sigma_full, self.n_threshold] = self.es(P, t, self.method);
            else
                warning('s must be > 0 but less than number of vectors in data');
            end
            
            self.eigvector = self.eigvector_full(:, 1:self.n_threshold);
            self.sigma = self.sigma_full(1:self.n_threshold, 1:self.n_threshold);
            self.ES = self.eigvector * self.sigma;
            self.invES = diag(1./diag(self.sigma)) * self.eigvector.';
            
            fprintf('number of valuable eigenvectors = %d\n', self.n_threshold);
        end
        
        function [P2,D2]=sortem(self, P,D)
            % this function takes in two matrices P and D, presumably the output
            % from Matlab's eig function, and then sorts the columns of P to
            % match the sorted columns of D (going from largest to smallest)
            
            [D_tmp, ind] = sort(diag(D),'descend');
            D2 = diag(D_tmp);
            P2 = P(:,ind);
        end
        
        function [eigvector_full, sigma_full, n_threshold] = es(self, P, t, method)
            % INPUT:
            % P - covariation matrix
            % t - threshold - ratio of valuable singular or eigenvalues to the total sum valuable singular or eigenvalues
            % if t <= 1 then number of eigenvectors is calculated
            % automatically and t stands for energy threshold
            % if t is an integer >= 2 then it represends number of
            % valuable eigenvectors
            % method - defines decomposition method
            % OUTPUT:
            % eigvector_full - matrix of all eigenvectors
            % sigma_full - diagonal matrix of all eigenvalues
            
            if (strcmp(method, 'svd'))
                [eigvector_full, sigma_full] = svd(P);
            elseif (strcmp(method, 'eig'))
                [eigvector_full_unsorted, sigma_full_unsorted] = eig(P);
                [eigvector_full, sigma_full] = self.sortem(eigvector_full_unsorted, sigma_full_unsorted);
            end
            sigma_full_sum = sum(abs(sigma_full(:)));
            n_threshold = 0;
            tmp_sum = 0;
            if (t <= 1)
                for i = 1:size(sigma_full, 1)
                    n_threshold = n_threshold + 1;
                    if (tmp_sum/sigma_full_sum >= t)
                        break;
                    end
                    tmp_sum = tmp_sum + abs(sigma_full(i, i));
                end
            else
                n_threshold = t;
            end
        end
        
        function Ds = par_to_solver(self, Dp)
            %INPUT:
            %Ds - vectors suitable for solver
            %OUTPUT:
            %Dp - parametrized vectors
            Ds = self.ES * Dp + repmat(self.mean_data, 1, size(Dp, 2));
        end
        function Dp = solver_to_par(self, Ds)
            Dp = self.invES * (Ds - repmat(self.mean_data, 1, size(Ds, 2)));
        end
        function Dp = solver_to_par_jaco(self, Ds) % parametrizing jacobian for ajoint
            if size(Ds,1) > size(Ds,2)
                Dp =  self.invES * Ds;
            else
                Dp =  self.invES * Ds';
            end
        end
    end
end