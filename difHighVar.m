%
dt = 0.01; da= 2*dt; 
dm = 0.05; ga = 0;% Adults will not grow
%

tmax = 26;

Bs = 0;
b1 = 0;
% smallest mass, juvenile growth
m0 = 0.01; 
dmk = dm;

% Max values
maxm = 5.5; % max mass larvae
maxAa = 31; % max age
maxAj = maxAa;
% Add one because iterate 1 is age zero

%
a2 = 0:da:maxAa ; %a2 =a2';

astepsA = length(a2); astepsJ = length(a2);


mass = m0:dm:(maxm);
msteps = length(mass);  % for adults

% 1 is the time steps I am keeping
tm1= 1;
tm = round(1/dt);
stept = astepsA;

% Age can reproduce, birth rate
 
%Adult
d0 = 0; 
d1 = 0;

Anewg = zeros(tmax + 1,msteps,1);

% This is to vary each parameter by 20\% below and above
II = ones(13,6) + [0 0 0 0 0 0 ; -0.2.*eye(6); eye(6).*0.2];

%
% In order to run the final run, dm needs to increase for CFL reasons. 
for kk = 1:12
    tic
IV1 = 78;

% All initial parameters 
%netpar = [ 22 1 5e-7 5 0.06 0.3].*II(kk,:); % Original 
%netpar = [ 20 1.2 6.2e-7 6 0.044 0.3].*II(kk,:); % 2 10.1094
%netpar = [ 20 1.2 6.2e-7 6 0.035 0.3].*II(kk,:); % 3 10.0962
netpar = [ 16 1.2 6.2e-7 6 0.035 0.3].*II(kk,:); % 4 10.0962
nn2 = netpar(1);
ms = netpar(2); 
as = netpar(3);
nn= netpar(4);
rs = netpar(5);
g1 = netpar(6);
dj = 0.05;
ys = 8;


gj = @(m, tt) g1.*m.^(2/3).*exp(-tt*rs);



% This is the current and future time step
% age,mass
A1 = zeros(astepsA, msteps);
%A1 = Aiv;
A2 = zeros(astepsA, msteps);
J2 = zeros(astepsJ, msteps); 

%J1 = Jiv; 
J1 = J2;




%
% Matrix of (age, mass,time)
Juv = zeros(astepsJ, msteps, ceil(tmax*tm1));
Anew = zeros((tmax + 1),msteps);



% Initial values 
 J1(1,1) = IV1/(dm*da);


% Here is the simpson weights
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

 [m22, a22] = meshgrid( mass, a2);
% This is the birth function
beta = @(A22) sum(sum((W.*(b1.*A22(find(a2 >= ys - 0.001,1):end, :).*m22(find(a2 >= ys - 0.001,1):end, :) + Bs.*A22(find(a2 >= ys - 0.001,1):end, :)))))/9;  %(a2>(10-da));
%beta = @(A22) Bs*sum((w(find(a2 >= ys,1):end).*sum(A22(find(a2 >= ys,1):end, :),2)))*da/3;  %(a2>(10-da));


m22a = m22;  a22a = a22; 
% Maturation function
ggam = (1- exp(-m22a.^nn2.*ms)).*(1- exp(-a22a.^nn.*as));
%
Juv(:,:, 1)  = J1;

    
%   %  
%     %


for i = dt:dt:tmax
     
    for k = 1
        j =1;
         A2(j,k) = A1(j,k) - dt*(A1(j,k)/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
         J2(j,k) = J1(j,k) - dt*(J1(j,k)/da + gj(((k-1).*dm + m0), i)*J1(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
       
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + gj(((k-1).*dm + m0),i)*J1(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        end
    end
    
    for k = 2:(msteps-1)
        j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)- gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)-gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
          
        end

    end
    
    
    % This is only so that I do not loose those extra large individuals. As
    % individuals should not be this large, if this for large, then the
    % parameters are unrealistic.
   k = msteps;
   
         j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da - gj(((k-2).*dm + m0),i)*J1(j,k-1)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da - (gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)))*J1(j,k));
        end
    
 
 cc = w'*(J1.*ggam).*dt./3;
 A2(1,:) = A2(1,:)+ cc;  
 
 Anew(ceil(i),:) = Anew(ceil(i),:)+ (cc);
    A1 = A2;
    J1 = J2;
    if(dt/da + gj(maxm,i+2)*dt/dm > 0.994)
        error('CFL issue')
    end


born = beta(A1);
%Anewd2(round(i/dt)) = J1(1,1);
J1(1,1) = J1(1,1) + born;
%Anewd(round(i/dt),:) = (cc);

   

end
     
        Anewg(:,:,1) = Anew;
       
        
        
        save('ParHighfit.mat',  'Anewg')
  toc
    
end




%%