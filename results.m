% hFig = figure('Name', 'Results');
hFig_c = figure('Name', 'Convergence');
set(hFig_c, 'Position',[0.13 0.11 0.295 0.815])
% hFig = figure(1);
% set(hFig, 'Position', [x y width height])

% set(0, 'CurrentFigure', hFig)
sNo = figure(10);
set(sNo, 'Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483])
% sNo = figure('Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483]);
xlim([0,length(Dinit)]);
hold on;
sPCA = figure(11);
set(sPCA, 'Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483])
% sPCA = figure('Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483]);
xlim([0,length(Dinit)]);
hold on;
sRot = figure(12);
set(sRot, 'Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483])
% sRot = figure('Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483]);
xlim([0,length(Dinit)]);
hold on;
sSwap = figure(13);
set(sSwap, 'Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483])
% sSwap = figure('Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483]);
xlim([0,length(Dinit)]);
hold on;
sExt = figure(14);
set(sExt, 'Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483])
% sExt = figure('Position',[0.13 0.225982532751092 0.334659090909091 0.526347975723483]);
xlim([0,length(Dinit)]);
hold on;
    
%     D = D_123(:,iii);
    
set(0, 'CurrentFigure', sNo)
plot(D, 'linewidth', 2, 'color', 'r'); 
plot(Dinit, 'linewidth', 2, 'color', 'g');
xlabel('x','FontSize',20)
ylabel('D(x)','FontSize',20)

set(0, 'CurrentFigure', sPCA)
plot(D, 'linewidth', 2, 'color', 'r'); 
plot(Dinit, 'linewidth', 2, 'color', 'g');
xlabel('x','FontSize',20)
ylabel('D(x)','FontSize',20)

set(0, 'CurrentFigure', sRot)
plot(D, 'linewidth', 2, 'color', 'r'); 
plot(Dinit, 'linewidth', 2, 'color', 'g');
xlabel('x','FontSize',20)
ylabel('D(x)','FontSize',20)

set(0, 'CurrentFigure', sSwap)
plot(D, 'linewidth', 2, 'color', 'r'); 
plot(Dinit, 'linewidth', 2, 'color', 'g');
xlabel('x','FontSize',20)
ylabel('D(x)','FontSize',20)

set(0, 'CurrentFigure', sExt)
plot(D, 'linewidth', 2, 'color', 'r'); 
plot(Dinit, 'linewidth', 2, 'color', 'g');
xlabel('x','FontSize',20)
ylabel('D(x)','FontSize',20)

%     start_comp;

set(0, 'CurrentFigure', sNo)

plot(a_no_cg(:,end), 'linewidth', 1, 'color', 'k'); 
legend('Real', 'Initial', 'Result')
title('Без параметризации')
saveas(sNo, 'figures/test/sNo_2.png', 'png')
%     plot(sNo, Dinit, 'linewidth', 2, 'color', 'g');
set(0, 'CurrentFigure', sPCA)
plot(a_pca_cg(:,end), 'linewidth', 1, 'color', 'k');
legend('Real', 'Initial', 'Result')
title('PCA параметризация, 6 компонент')
saveas(sPCA, 'figures/test/sPCA_2.png', 'png')
%     plot(sPCA, Dinit, 'linewidth', 2, 'color', 'g');
set(0, 'CurrentFigure', sRot)
plot(a_rot_cg(:,ex_rot_cg), 'linewidth', 1, 'color', 'k'); 
legend('Real', 'Initial', 'Result')
title('Rotation - PCA параметризация, 6 компонент')
saveas(sRot, 'figures/test/sRot_2.png', 'png')
%     plot(sRot, Dinit, 'linewidth', 2, 'color', 'g');
set(0, 'CurrentFigure', sSwap)
plot(a_swap_cg(:,ex_swap_cg), 'linewidth', 1, 'color', 'k'); 
legend('Real', 'Initial', 'Result')
title('Swap - PCA параметризация, 6 компонент')
saveas(sSwap, 'figures/test/sSwap_2.png', 'png')
%     plot(sSwap, Dinit, 'linewidth', 2, 'color', 'g');
set(0, 'CurrentFigure', sExt)
plot(a_ext_cg(:, ex_ext_cg), 'linewidth', 1, 'color', 'k');
legend('Real', 'Initial', 'Result')
title('Extension - PCA параметризация, 6-11 компонент')
saveas(sExt, 'figures/test/sExt_2.png', 'png')
%     plot(sExt, Dinit, 'linewidth', 2, 'color', 'g');

set(0, 'CurrentFigure', hFig_c)
set(hFig_c, 'Position',[0.13 0.11 0.295 0.815])
semilogy(b_no_cg, 'b');
hold on;
semilogy(b_pca_cg, 'k', 'linewidth',4);
semilogy(b_rot_cg, 'r', 'linewidth',2);
semilogy(b_swap_cg, 'm');
semilogy(b_ext_cg, 'g');
% Create xlabel
xlabel('№ итерации');

% Create ylabel
ylabel('log Значение целевой функции');


legend('Полная размерность', 'PCA','ROTATION', 'SWAP', 'EXTENSION');
title('Сходимость методов для модели 2')
saveas(hFig_c, 'figures/test/conv_2.png', 'png')



