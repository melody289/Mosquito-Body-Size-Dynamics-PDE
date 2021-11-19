kkj =1;  %  10.6952

%
dt = 0.01; da= 2*dt; 
dm = 0.05; ga = 0;% Adults will not grow
%

tmax2 = 20;

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

% No Adults
d0 = 0; 
d1 = 0;

Juvgl = zeros(astepsJ, msteps, 21,5);
Anewgl = zeros(22,msteps,5);
% 2nd smaller 3rd larger  10.741
%
II = ones(13,6) + [0 0 0 0 0 0 ; -0.2.*eye(6); eye(6).*0.2];
%
for kk = 1:13 
    tic
IV2 = 26;

%netpar = [ 8    0.4250    0.0006    2.5    0.0510    0.28 ].*II(kk,:);
%netpar = [ 14    0.62    9e-4   2.5    0.05   0.28 ].*II(kk,:);
%
%netpar = [  2.5    0.1    7e-4    2.5    0.03    0.4 ].*II(kk,:);  % a  6.0880
%netpar = [  3    0.12    8e-4    3.5    0.03    0.4 ].*II(kk,:);  % b 5.7393
%netpar = [  3    0.12    8e-4    2.8    0.03    0.4 ].*II(kk,:);  % c  5.724
%netpar = [  3    0.096    9e-4    2.8    0.03    0.4 ].*II(kk,:);  % d 5.7028
%netpar = [  3    0.06   11e-4    2.8    0.03    0.4 ].*II(kk,:);  % f 5.6519
%netpar = [  3.6    0.06   11e-4    2.8    0.03    0.4 ].*II(kk,:);  % g 5.6341 3, 6, 10
netpar = [  3.6    0.038   2e-3    2.8    0.028    0.4 ].*II(kk,:);  % h 5.6011 0.028 back to 0.03


nn2 = netpar(1); % mass exponent
ms = netpar(2); % mass coefficient
as = netpar(3); % age coef
nn= netpar(4); % age exponent
rs = netpar(5);
g1 = netpar(6); % max growth parameter
dj = 0.025;
ys = 8;


gj = @(m, tt) g1.*m.^(2/3).*exp(-tt*rs);


% This is the current and future time step
% age,mass
A1b = zeros(astepsA, msteps);
A2b = zeros(astepsA, msteps);
J2b = zeros(astepsJ, msteps); 
J1b = J2b;


%
% Matrix of (age, mass,time)
Juvl = zeros(astepsJ, msteps, ceil(tmax2*tm1));
Anewl = zeros(22,msteps);

%
% Initial values 
 J1b(1,1) = IV2/(dm*da);

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


m22a = m22;  a22a = a22; 
ggam = (1- exp(-m22a.^nn2.*ms)).*(1- exp(-a22a.^nn.*as));
%

Juvl(:,:, 1)  = J1b;
        

for i = dt:dt:(tmax2)
     
    for k = 1
        j =1;
         A2b(j,k) = A1b(j,k) - dt*(A1b(j,k)/da + ga*A1b(j,k)/dm + (d0+ da*d1*(j-1))*A1b(j,k));
         J2b(j,k) = J1b(j,k) - dt*(J1b(j,k)/da + gj(((k-1).*dm + m0), i)*J1b(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
       
        for j = 2:astepsJ
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da + ga*A1b(j,k)/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da + gj(((k-1).*dm + m0),i)*J1b(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        end
    end
    
    for k = 2:(msteps-1)
        j = 1;
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k))/da + ga*(A1b(j,k)-A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k))/da + (gj(((k-1).*dm + m0),i)*J1b(j,k)- gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        
        for j = 2:astepsJ
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da + ga*(A1b(j,k)-A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da + (gj(((k-1).*dm + m0),i)*J1b(j,k)-gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
          
        end

    end
    
    
    % This is only so that I do not loose those extra large individuals. As
    % individuals should not be this large, if this for large, then the
    % parameters are unrealistic.
   k = msteps;
   
         j = 1;
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k))/da - ga*(A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k))/da - gj(((k-2).*dm + m0),i)*J1b(j,k-1)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        
        for j = 2:astepsJ
           
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da - ga*(A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da - (gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)))*J1b(j,k));
        end
    
   
 ccb = w'*(J1b.*ggam).*dt./3;
 A2b(1,:) = A2b(1,:)+ ccb;  
 

 Anewl(ceil(i),:) = Anewl(ceil(i),:)+ (ccb);
    A1b = A2b;
    J1b = J2b;
    if(dt/da + gj(maxm,i+2)*dt/dm > 0.994)
        error('CFL issue')
    end


    % no births
% born = beta(A1);
% J1(1,1) = J1(1,1) + born;
% born = beta(A1b);
% J1b(1,1) = J1b(1,1) + born;




%Jnewd(round(i/dt)) = born;

    if(mod(round((i + 0.001),2),tm1) ==0)
        
        Juvl(:,:, round(i)+1)  = J1b;
    end
    
end




       
        Anewgl(:,:,kk) = Anewl;
       
        
     

LowS = zeros(16,17);
    AAl = Anewl(:,:).*dm.*da;
    for  jj = 1:17
        for ii = 1:16
            LowS(ii,jj) = sum(sum(AAl(ii,(jj*6+1):(jj*6+6)) ));
        end
        
    end

        save('ParLowFith.mat',    'Anewgl')

         er= norm(LowM-LowS, 'fro');
       
      
  toc

    
end



%%


