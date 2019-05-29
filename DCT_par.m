classdef DCT_par %< handle
    properties
%         N; % number of state vectors in parametrized ensemble
        dct_id; 
        cosma_full; % full cosine matrix
        cosma % truncated cosin matrix
%         mean_data; % mean value of data
        data; % input dataset
    end
    methods
        function self = DCT_par(data, t)
            disp('Initializing DCT');
            
            self.cosma_full = dct1d(length(data));
            data_cos = self.cosma_full * data;
            [~, idxs] = sort(data_cos, 'descend');
            idxs_use = idxs(1:t);
            self.cosma = self.cosma_full(idxs_use,:);
            self.dct_id = idxs_use;
        end
       
        
        function Ds = par_to_solver(self, Dp)
            %INPUT:
            %Ds - vectors suitable for solver
            %OUTPUT:
            %Dp - parametrized vectors
            Ds = (self.cosma)' * Dp;
        end
        function Dp = solver_to_par(self, Ds)
            Dp = self.cosma * Ds;
        end
        function Dp = solver_to_par_jaco(self, Ds) % parametrizing jacobian for ajoint
            if size(Ds,1) > size(Ds,2)
                Dp = self.cosma * Ds;
            else
                Dp = self.cosma * Ds';
            end
        end
    end
end