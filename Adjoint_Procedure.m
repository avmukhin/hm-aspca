function gradient_return = Adjoint_Procedure(init_struct, HistoryData, D)

%%
% Main purpose of this function is to calculate gradient of objective
% function using adjoint procedure
%%
[resM, ~, Rx, G] = solve(init_struct, D);

resS = HistoryData;  % y_data

L = @(t) (resM(:, t) - resS(:, t))'; % derivative of objective function

% [Nx, Nt] = size(resS);
Nx = init_struct.Nx; Nt = init_struct.Nt;
lamb = zeros(Nx, Nt); %adjoint vectors

J = @(t) jacobian(init_struct, resM(:, t), D);

lamb(:, Nt) = - L(Nt) * J(Nt)^(-1);
% lamb(:, Nt) = zeros(Nx, 1); % _refactoring_ ????

for t = Nt-1:-1:1
    lT = (lamb(:, t+1))';
    jI = J(t)^(-1);
    dL = L(t);
    lamb(:, t) = (-(dL + lT * Rx) * jI)';
end

gradient = zeros(1, init_struct.Nx-1);
for t = 1:Nt
    Gt = G{t};
    lamb_t = lamb(:, t);
    add = lamb_t' * Gt;
    gradient = gradient + add;
end

gradient_return = - gradient;

% gradient_return = D' .* gradient(1:init_struct.NxD);
end