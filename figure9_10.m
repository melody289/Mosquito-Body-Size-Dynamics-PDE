load('cyclicgroup.mat')
boxplot(mass25')
set(gca, 'FontSize', 16, 'Position', [0.1000    0.1303    0.8500    0.7947])
set(gcf, 'Position', [ 116   436   884   355])
ylabel('Mass > 2.5 mg')


figure 
boxplot(mass15')
set(gca, 'FontSize', 16, 'Position', [0.1000    0.1303    0.8500    0.7947])
set(gcf,'Position', [ 116   436   884   355])
ylabel('Mass > 1.5 mg')
xlabel('Period of Cyclic Resources (Days)')


%%


boxplot(age10')
set(gca, 'FontSize', 16, 'Position', [0.1000    0.1303    0.8500    0.7947])
set(gcf, 'Position', [ 116   436   884   355])
ylabel('Age > 10 days')


figure 
boxplot(age7')
set(gca, 'FontSize', 16, 'Position', [0.1000    0.1303    0.8500    0.7947])
set(gcf,'Position', [ 116   436   884   355])
ylabel('Age > 7 days')
xlabel('Period of Cyclic Resources (Days)')
