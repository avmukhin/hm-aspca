function jaco = jacobian(init_struct, u, D)
Nx = init_struct.Nx; 
h = init_struct.L / (init_struct.Nx - 1);
tau = init_struct.T / (init_struct.Nt - 1);
uL = init_struct.uL;
uR = init_struct.uR;

jaco = zeros(Nx);

%% u{L/R} = du{left/right} / dh
jaco(1,1) = h^2/tau ...
    - D(1) * ( 2*u(1) * (u(2) - u(1)) - (u(1)^2 + u(2)^2) ) ...
    + D(1) * ( 2*u(1)) * ( 2*uL*h);
jaco(1,2) = - D(1) * ( 2*u(2) * (u(2) - u(1)) + (u(1)^2 + u(2)^2) );

jaco(Nx, Nx) =  h^2/tau ...
    + D(Nx-1) * ( 2*u(Nx) * (u(Nx) - u(Nx-1)) + (u(Nx)^2 + u(Nx-1)^2) ) ...
    - D(Nx-1) * ( 2*u(Nx) ) * ( 2*uR*h );
jaco(Nx, Nx-1) = D(Nx-1) * ( 2*u(Nx-1) * (u(Nx) - u(Nx-1)) - (u(Nx)^2 + u(Nx-1)^2) );

%%
for i = 2:Nx-1
    jaco(i,i-1) = D(i-1) * ( 2*u(i-1)*(u(i) - u(i-1)) - (u(i)^2 + u(i-1)^2) );
    jaco(i,i) = 2*h^2/tau ...
        - D(i) * ( 2*u(i)*(u(i+1) - u(i)) - (u(i+1)^2 + u(i)^2) ) ...
        + D(i-1) * ( 2*u(i)*(u(i) - u(i-1)) + (u(i)^2 + u(i-1)^2) ) ;
    jaco(i,i+1) = - D(i) * ( 2*u(i+1)*(u(i+1) - u(i)) + (u(i+1)^2 + u(i)^2) ) ;
end

end