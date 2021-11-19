dt = 0.01; da= 2*dt; 
dm = 0.02; ga = 0;% Adults will not grow
%
tmax = 150;

rs = [ 0.02, 0.05, 0.2];

Bs = 10*dt;
b1 = 1*dt;
% smallest mass, juvenile growth
m0 = 0.01; 


% Parameters
% For mass maturation exponent and coef
nn2 = 10.8;
ms = 0.6;
% For age maturation exponent and coef
nn = 3;
as = 6e-4;
% growth parameter
g1 = 0.33; 


% this is the rounded value of mass that is kept
dmk = dm;

% Max values
maxm = 5; % max mass larvae
maxAa = 31; % max age
maxAj = maxAa;
% Add one because iterate 1 is age zero

%
a2 = 0:da:maxAa ; %a2 =a2';

astepsA = length(a2); astepsJ = length(a2);
%agem = find(a2> 4.99 , 1);

mass = m0:dm:(maxm);
msteps = length(mass);  % for adults
%
% 1 is the time steps I am keeping
tm1= 1;
tm = round(1/dt);
stept = astepsA;

%
% Parameters
% Age can reproduce, birth rate
ys = 8;   

% death rates
d0 = 0.05; dj = 0.05;
d1 = 0.025;

% Here is the simpson weights
if(mod(astepsA,2) ==0)
    disp('A Not even')
    return
elseif(mod(astepsJ,2) ==0)
    disp('J Not even')
    return
end
w = 2* ones ( astepsA  ,1); % column vector , starts all 2 ’s
w (2:2:end)=4; 
w(1) = 1; w(end) = 1; % set ends to 1 ’s


w2 = w(1:(msteps)); w2(end) = 1;
W = w*w2';
W = W(find(a2 >= ys - 0.001,1):end,:);
 w = w(1:astepsJ); w(end) = 1;

%
%
%


Adg = zeros(astepsA, msteps, (tmax + 1),3);
%Juvg = zeros(astepsJ, msteps, (tmax + 1),3);
%Anewg = zeros((tmax + 1),msteps,3);

 [m22, a22] = meshgrid( mass, a2);
% This is the birth function
beta = @(A22) sum(sum((W.*(b1.*A22(find(a2 >= ys - 0.001,1):end, :).*m22(find(a2 >= ys - 0.001,1):end, :) + Bs.*A22(find(a2 >= ys - 0.001,1):end, :)))))/9;  %(a2>(10-da));
%beta = @(A22) Bs*sum((w(find(a2 >= ys,1):end).*sum(A22(find(a2 >= ys,1):end, :),2)))*da/3;  %(a2>(10-da));
ggam = (1- exp(-m22.^nn2.*ms)).*(1- exp(-a22.^nn.*as));
%

for kkj = 1:3
    tic
A1 = Aiv; %zeros(astepsA, msteps);
A2 = zeros(astepsA, msteps);
J2 = zeros(astepsJ, msteps); 
J1 = Jiv; %J2;

gj = @(m, tt) g1.*m.^(2/3);

rt = @(tt,mm) 0.25.*(maxm-mm).*(1- exp(-tt*rs(kkj))).^5./(maxm.*((1- exp(-tt*rs(kkj))).^5 + 0.75^5)) + 0.02;

%
% Matrix of (age, mass,time)
Ad = zeros(astepsA, msteps, ceil(tmax*tm1));
Juv = zeros(astepsJ, msteps, ceil(tmax*tm1));
Anew = zeros((tmax + 1),msteps);
%




Ad(:,:, 1)  = A1;
Juv(:,:, 1)  = J1;
        
%   %  
%     %
for i = dt:dt:(tmax)
     
    for k = 1
        j =1;
         A2(j,k) = A1(j,k) - dt*(A1(j,k)/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
         J2(j,k) = J1(j,k) - dt*(J1(j,k)/da + gj(((k-1).*dm + m0), i)*J1(j,k)/dm + (dj*rt(i,((k-1).*dm + m0)) +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
       
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + gj(((k-1).*dm + m0),i)*J1(j,k)/dm + (dj*rt(i, ((k-1).*dm + m0)) +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        end
    end
    
    for k = 2:(msteps-1)
        j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)- gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj*rt(i, ((k-1).*dm + m0)) +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)-gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj*rt(i, ((k-1).*dm + m0)) +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
          
        end

    end
    
    
    % This is only so that I do not loose those extra large individuals. As
    % individuals should not be this large, if this for large, then the
    % parameters are unrealistic.
   k = msteps;
   
         j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da - gj(((k-2).*dm + m0),i)*J1(j,k-1)/dm + (dj*rt(i, ((k-1).*dm + m0)) +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da - (gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj*rt(i,((k-1).*dm + m0)) +   (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        end
    
 
 cc = w'*(J1.*ggam).*dt./3;
 A2(1,:) = A2(1,:)+ cc;  
 
 Anew(ceil(i),:) = Anew(ceil(i),:)+ (cc);
    A1 = A2;
    J1 = J2;
    if( sum(sum(A1 <0)) >0)
        i
        break
    elseif(dt/da + gj(maxm,i+2)*dt/dm > 0.994)
        error('CFL issue')
    end


born = beta(A1);
%Anewd2(round(i/dt)) = J1(1,1);
J1(1,1) = J1(1,1) + born;
%Anewd(round(i/dt),:) = (cc);




%Jnewd(round(i/dt)) = born;

    if(mod(round((i + 0.001),2),tm1) ==0)
        Ad(:,:, round(i)+1)  = A1;
        Juv(:,:, round(i)+1)  = J1;
    end

end
        Adg(:,:, :,kkj)  = Ad;
       
        toc
        save('C1par_decayDeath.mat',  'Adg')
end
            
 %%