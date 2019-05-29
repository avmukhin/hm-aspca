function [reslist, convlist, exit_iter] = optimize_CG(init_struct, HistoryData, x_prev, max_iter, eps, par, method, useplot, spca, spca_method)
alpha0 = 2;
alpha = alpha0;
global s1;
global s2;
%%
% spca_method = 'rotate';
% spca_method = 'swap';
% spca_method = 'extend';

%%
color_palette = ['b', 'k', 'r', 'g', 'y', 'c', 'm', 'w'];

%%
x_prev_p = parametrize(x_prev, method, 'p', par);
x_prev = parametrize(x_prev_p, method, 's', par);
x_min = 0.001*ones(size(x_prev));
x_max = 2*ones(size(x_prev));

alpha_loop = 0;
%%
f_best = Objective_Function(init_struct, HistoryData, x_prev);
g_prev = Adjoint_Procedure(init_struct, HistoryData, x_prev); 
g_prev = g_prev / (norm(g_prev') + 1e-12);
g_prev_p = parametrize(g_prev, method, 'p_jaco', par);
x_prev_p = parametrize(x_prev, method, 'p', par);
%%
convlist = f_best;
reslist = [];
real_itr = 0;
itr = 0;
exit_iter = [];
if useplot
    curr_title = get(gca, 'Title');
    curr_title = curr_title.String;
end

while itr < max_iter
    
    x_new_p = x_prev_p + alpha * g_prev_p;
    x_new = parametrize(x_new_p, method, 's', par);
    
    if any(x_new > x_max) || any(x_new < x_min)
        x_new = x_new - min(x_new(:));
        x_new = (x_new/(max(x_new(:))-min(x_new(:))))*(max(x_max) - min(x_min));
        x_new = x_new + min(x_min);
    end
    
    if useplot
        if exist('sp', 'var')
            delete(sp) 
        end
        sp = plot(s1, x_new, 'color', 'k');
        drawnow limitrate;
    end
        
    f_new = Objective_Function(init_struct, HistoryData, x_new);
    
    if (abs(f_new - f_best) < 1e-12) || (alpha_loop == 2)
        fprintf('Objective convergence \n');
        fprintf('iteration num is %d\n', itr);
        exit_iter(end+1) = itr;
        if spca
            alpha_loop = 0;
            spca = spca - 1;
%             figure(40)
%             plot(x_new);
%             hold on;
%             plot(par.ES(:, 1:end-1) * x_prev_p(1:end-1))
            [x_prev_p, par] = swap_pca(x_prev_p, g_new, par, spca_method, 1);
            
%             figure(40)
%             plot(par.ES * x_prev_p);
            
            alpha = alpha0;
            x_prev = parametrize(x_prev_p, method, 's', par);
            if any(x_prev > x_max) || any(x_prev < x_min)
                x_prev = x_prev - min(x_prev(:));
                x_prev = (x_prev/(max(x_prev(:))-min(x_prev(:))))*(max(x_max) - min(x_min));
                x_prev = x_prev + min(x_min);
            end
            f_best = Objective_Function(init_struct, HistoryData, x_prev);
            f_new = f_best;
            convlist = cat(2, convlist, f_new);
            reslist = cat(2, reslist, x_prev);
            if strcmp(spca_method, 'extend')
                x_new_p = x_prev_p;
            end
        else
            break
        end
     
    end
    
    if (f_new <= f_best)
        
        convlist = cat(2, convlist, f_new);
        reslist = cat(2, reslist, x_new);
        f_best = f_new;
        x_prev_p = x_new_p;
        
        itr = itr + 1;
        
        g_new = Adjoint_Procedure(init_struct, HistoryData, x_new);
        g_new = g_new / (norm(g_new') + 1e-12);
        g_new_p = parametrize(g_new, method, 'p_jaco', par);
        
        if strcmp(spca_method, 'extend')
            g_prev_p = g_new_p;
        else
            pr = (g_new_p' * (g_new_p - g_prev_p)) / (g_prev_p' * g_prev_p);
            w = max(0, pr);
            g_prev_p = g_new_p + g_prev_p * w;
        end
       
        
        if useplot
            
%             if exist('g', 'var')
%                 delete(g)
%             end
%             g = plot(s2, g_prev, 'color', 'b');
%             drawnow limitrate;    
%             if exist('gn', 'var')
%                 delete(gn)
%             end
%             gn = plot(s2, g_new, 'color', 'r');
            try
                semilogy(s2, itr, f_new, 'o-', 'color', color_palette(spca - length(color_palette)*floor(spca/length(color_palette))));
            catch
                semilogy(s2, itr, f_new, 'o-', 'color', 'r');
            end
            xlim([1,itr+1]);
            drawnow limitrate;

%             title(s2, ['red - real, blue - CG   ', 'alpha = ', num2str(alpha), ]);
            title(s2, ['Convergence']);


%             
%             if exist('p', 'var')
%                 delete(p)
%             end
%             p = plot(s1, x_new, 'color', 'b');
            drawnow limitrate;
            title(s1, [curr_title, 'f = ', num2str(f_new), ', Iteration ', num2str(itr)]);
            pause(0.5);
        end
        alpha = min(alpha0, alpha*(2^5));
        
    else
        alpha = alpha / 2;
        
        
        if alpha < 1e-8
            if alpha_loop == 0
                alpha_loop = 1;
                alpha = alpha0;
%                 alpha = min(alpha0, alpha*(2^5));

                g_prev = Adjoint_Procedure(init_struct, HistoryData, x_prev);
                g_prev = g_prev / (norm(g_prev') + 1e-12);
                g_prev_p = parametrize(g_prev, method, 'p_jaco', par);
            else
                alpha_loop = 2;
            end

        end
    end
       
    real_itr = real_itr + 1;
   

end

end