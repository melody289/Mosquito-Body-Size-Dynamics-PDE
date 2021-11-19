
% Get Error, must have all files
% lowM, HighM from data - dataFormat
% Anewg13, Anewgl13
load('VarSepFit.mat')
ch = 1;
ParError
% base nn2, ms, as, nn, rs, g1
ord = [1, 4, 10, 5, 11,  3 9, 2 8 7 13, 6 12];
s = {'Base', '0.8n_2', '0.8\alpha_2', '0.8\alpha_1', '0.8n_1', '0.8r_*', '0.8\gamma_j', '1.2n_2', '1.2\alpha_2', '1.2\alpha_1',  '1.2n_1', '1.2r_*', '1.2\gamma_j'};
b = bar(err(ord));
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xticklabels(s(ord))
axis([0.5, 13.5, 9, 18])
title('High Density Samples')
xlabel('Parameter Varied')
ylabel('Error')
hold on
plot([0, 14], [10.0962, 10.0962], '--k')
set(gca, 'FontSize', 16)
set(gca, 'Position', [0.0800    0.1283    0.8500    0.7970])
set(gcf, 'Position', [457   357   924   420])


figure

ch = 0;
ParError


b = bar(err(ord));
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xticklabels(s(ord))
axis([0.5, 13.5, 5,6])
title('Low Density Samples')
xlabel('Parameter Varied')
ylabel('Error')
hold on
plot([0, 14], [5.5997, 5.5997], '--k')
set(gca, 'FontSize', 16)
set(gca, 'Position', [0.0800    0.1283    0.8500    0.7970])
set(gcf, 'Position', [457   357   924   420])
