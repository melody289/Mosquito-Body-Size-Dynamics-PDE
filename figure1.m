% The creates figure 1
% First we need to load data, both high (HighM) and low (LowM).
% HighM is a matrix of days in increments of two, so that the first box is
% those from days 0 to 2, and the second box is those just after 2 to day
% 4,  LowM days are incremented by one day 
% mass is incremented from 0.3 to 5.25 by 0.3
% See file dataFormat for all necessary files
load('VarSepFit.mat')
ch = 1;
ParError
ch = 0;
ParError


ac(1) = subplot(2,3,1);
imagesc(HighM)
colormap parula
colorbar
title('High Density Data')
ylabel('Development Time')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16)
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
set(gca, 'Position',[0.1 0.58 0.22 0.34])
yticks(1.5:2:28)
yticklabels({'2', '6',  '10',  '14', '18', '22', '26'})




ac(2) = subplot(2,3,2);
imagesc(HighS(:,:,1))
colorbar
title('High Density Model')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16)
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
set(gca, 'Position',[0.4 0.58 0.22 0.34]);
yticks(1.5:2:28)
yticklabels({'2', '6',  '10',  '14', '18', '22', '26'})



ac(3) = subplot(2,3,3);
imagesc(HighS(:,:,1)- HighM )
caxis([-1,1])
colorbar
title('Error')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16)
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
set(gca,'Position',[0.71 0.58 0.22 0.34])
yticks(1.5:2:28)
yticklabels({'2', '6',  '10',  '14', '18', '22', '26'})



ac(4) = subplot(2,3,4);
imagesc(LowM)
colorbar
title('Low Density Data')
ylabel('Development Time')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16, 'Position',[0.1 0.1 0.22 0.34])
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
yticks(2.5:2:16.5)
yticklabels({'2', '4',  '6',  '8', '10', '12',  '14', '16'})


ac(5) = subplot(2,3,5);
imagesc(LowS(:,:,1))
colorbar
title('Low Density Model')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16,'Position',[0.41 0.1 0.22 0.34])
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
yticks(2.5:2:16.5)
yticklabels({'2', '4',  '6',  '8', '10', '12',  '14', '16'})


ac(6) = subplot(2,3,6);
imagesc(LowS(:,:,1)-LowM)
caxis([-1,1])
colorbar
title('Error')
xlabel('Mass')
set(gca,'YDir','normal', 'FontSize', 16, 'Position',[0.71 0.1 0.22 0.34])
xticks(0.5:4:17)
xticklabels({'0.45',  '1.65',  '2.85',   '4.05',  '5.25'})
yticks(2.5:2:16.5)
yticklabels({'2', '4',  '6',  '8', '10', '12',  '14', '16'})

colormap(ac(3), pink)
colormap(ac(6), pink)
set(gcf, 'Position', [34   245   966   552])