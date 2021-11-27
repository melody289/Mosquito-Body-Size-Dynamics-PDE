dt = 0.01; da= 2*dt; 
dm = 0.02; 

% Maximum time
tmax = 175;


% smallest mass, juvenile growth
m0 = 0.01; 


% Parameters
% period of cycle
omg = 1:2:31; 

% Max values
maxm = 5; % max mass larvae
maxAa = 31; % max age
maxAj = maxAa;
% Add one because iterate 1 is age zero

%
a2 = 0:da:maxAa ; %a2 =a2';

astepsA = length(a2); 

mass = m0:dm:(maxm);
msteps = length(mass);
m25 = mass > 2.4999;
m15 = mass > 1.4999;
a10 = a2 > 9.999;
a7 = a2 > 6.999;
mass25 = zeros(16,176);
mass15 = zeros(16,176);
age7 = zeros(16,176);
age10 = zeros(16,176);
for kkj = 1:16
ss = sprintf('C1par_cylicmass_%d.mat', omg(kkj));
load(ss)


massT = squeeze(sum(Ad));
mass25(kkj,:) = sum(massT(m25,:))./sum(massT);
mass15(kkj,:) = sum(massT(m15,:))./sum(massT);

ageT = squeeze(sum(Ad,2));
age10(kkj,:) = sum(ageT(a10,:))./sum(ageT);
age7(kkj,:) = sum(ageT(a7,:))./sum(ageT);

end
%%
