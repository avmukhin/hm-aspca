function [un1n, iteration, g] = newton_iteration(init_struct, un, D)
iteration = 0;
un1n = un;
min_delta = Inf;
while iteration <= init_struct.maxNewton
    j = jacobian(init_struct, un1n, D);
    
    [f, g] = FGvec(init_struct, un1n, un, D);
    delta_u = j \ f;
    
    u = un1n - delta_u;
    if (norm(f)==0)
        un1n = u;
        break
    end
    norm_delta = norm( (u'-un1n')/(un1n' + 1e-8*ones(size(u'))) );
    if norm_delta < init_struct.eps
        un1n = u;
        break
    end
    if norm_delta < min_delta
        un1n = u;
    end
    iteration = iteration + 1;
end
end