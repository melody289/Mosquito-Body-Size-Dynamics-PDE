%% Load data

data = readtable('Larval_counts_emergence_time_and_adult_body_mass.xlsx', 'sheet', 'VarE');



% Remove males
Ds = data.sex;
Index = find(contains(Ds,'F'));
data2 = data(Index,:);

% Remove na
Ds = data2.b_mass_em;
Index = (contains(Ds,'n'));
data2 = data2(~Index,:);
data2.b_mass_em = str2double(data2.b_mass_em); % change to numerical value
data2.dev_time = data2.dev_time./24 + 3; % change to days and three days for estimated egg days

% Separate by low and high
Ds = data2.lar_density;
Index = (contains(Ds,'low'));
datal = data2(Index,:);
datah = data2(~Index,:);

% group by mass and development time
[N,Xedgel,Yedgel] = histcounts2(datal.dev_time, datal.b_mass_em);

% 20.6 is the average survival of the data for low density IV = 26
N = 20.6.*N./sum(sum(N));
LowM = zeros(16,17);
LowM(10:15,6:11) = N;


% Group for high density 
[N,Xedges,Yedges] = histcounts2(datah.dev_time, datah.b_mass_em);
N = 50.5.*N./sum(sum(N));


HighM = zeros(13,17);
HighM(5:12,2:7) = N;

clearvars N data data2 data2  Ds Index Xedges Yedges  datah datal