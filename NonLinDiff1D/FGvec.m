function [F, G] = FGvec(init_struct, u, un, D)
Nx = init_struct.Nx;
NxD = init_struct.NxD;
uL = init_struct.uL;
uR = init_struct.uR;
h = init_struct.L / (init_struct.Nx - 1);
tau = init_struct.T / (init_struct.Nt - 1);

F = zeros(Nx, 1);
G = zeros(Nx, NxD);

%% u{L/R} = du{left/right} / dh
F(1) = ( u(1)-un(1) ) * h^2/tau ...
    - D(1) * (u(1)^2 + u(2)^2) * (u(2) - u(1)) ...
    + D(1) * (u(1)^2) * (2*uL*h);

F(Nx) = ( u(Nx)-un(Nx) ) * h^2/tau ...
    + D(Nx-1) * (u(Nx)^2 + u(Nx-1)^2) * (u(Nx) - u(Nx-1)) ...
    - D(Nx-1) * (u(Nx)^2) * (2*uR*h);

G(1,1) = - (u(1)^2 + u(2)^2) * (u(2) - u(1)) ;%+ (u(1)^2) * (2*uL*h);
G(Nx,NxD) = (u(Nx)^2 + u(Nx-1)^2) * (u(Nx) - u(Nx-1)) ;%- (u(Nx)^2) * (2*uR*h);
% G(Nx-1,Nx-2) = - (u(Nx)^2) * (2*uR);

%%
for i=2:Nx-1
    F(i) = (u(i) - un(i)) * 2*h^2/tau - ...
        D(i)*( u(i+1)^2 + u(i)^2 ) * ( u(i+1) - u(i) ) + ...
        D(i-1)*( u(i)^2 + u(i-1)^2 ) * (u(i) - u(i-1)) ;

end

for i=2:Nx-1
    G(i,i) = - ( u(i+1)^2 + u(i)^2 ) * ( u(i+1) - u(i) );
    G(i, i-1) = ( u(i)^2 + u(i-1)^2 ) * (u(i) - u(i-1));   
end

end