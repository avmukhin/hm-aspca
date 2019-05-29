function costval = Objective_Function(init_struct, HistoryData, D)


[resM, ~, ~, ~] = solve(init_struct, D);
resS = HistoryData;  
result = 0;

L = @(t) (resM(:, t) - resS(:, t))';


for j = 1:size(resS, 2)
    dif = L(j);
    result = result + dif * dif.';
end
costval = result / 2;
end

