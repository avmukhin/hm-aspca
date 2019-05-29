function [C, tau, Rx, G] = solve(init_struct, D)
Nx = init_struct.Nx;
Nt = init_struct.Nt;
C = init_struct.uinit;
C_prev = init_struct.uinit;
% Time = 0;
tau = init_struct.T / (Nt - 1);
% timeit = [];
% timesteps = tau;

G = {};

for tt=1:Nt
    
    [C_new, ~, g] = newton_iteration(init_struct, C_prev, D);
%     G = cat(3, G, g);
    G{end+1} = g;
    
    C = cat(2, C, C_new);
    C_prev = C_new;
    
    tt = tt + 1;
%     Time = Time + tau;
    
%     if Time > init_struct.T
%         timeit = cat(2, timeit, init_struct.T);
%     else
%         timeit = cat(2, timeit, Time);
%         timesteps = cat(2, timesteps, tau);
%     end
    
end
h = init_struct.L / (init_struct.Nx - 1);
tau = init_struct.T / (init_struct.Nt - 1);
Rx = - 2*h^2/tau * eye(init_struct.Nx);
Rx(1,1) = 0.5*Rx(1,1); Rx(Nx,Nx) = 0.5*Rx(Nx,Nx);
C = C(:,2:end);
end