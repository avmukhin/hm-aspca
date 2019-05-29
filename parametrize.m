% function D_output = parametrization_method(method, ps, D_input, par)
function D_output = parametrize(D_input, method, ps, par)
if (strcmp(method, 'No'))
    if size(D_input,1) > size(D_input,2)
        D_output = D_input;
    else
        D_output = D_input';
    end
elseif (strcmp(method, 'PCA'))
    if (strcmp(ps, 'p_jaco'))
        D_output = par.solver_to_par_jaco(D_input);
    elseif (strcmp(ps, 'p'))
        D_output = par.solver_to_par(D_input);
    elseif (strcmp(ps, 's'))
        D_output = par.par_to_solver(D_input);
    end
    
elseif (strcmp(method, 'S-PCA'))
    if (strcmp(ps, 'p_jaco'))
        D_output = par.solver_to_par_jaco(D_input);
    elseif (strcmp(ps, 'p'))
        D_output = par.solver_to_par(D_input);
    elseif (strcmp(ps, 's'))
        D_output = par.par_to_solver(D_input);
    end
    
elseif (strcmp(method, 'DCT'))
    
    if (strcmp(ps, 'p_jaco'))
        D_output = par.solver_to_par_jaco(D_input);
    elseif (strcmp(ps, 'p'))
        D_output = par.solver_to_par(D_input);
        cosma = dct1d(length(D_input)) * D_input;
    elseif (strcmp(ps, 's'))
        D_output = par.par_to_solver(D_input);
    end
end

end