
% Step sizes
dt = 0.01; da= 2*dt; 
dm = 0.02; ga = 0;% Adults will not grow

% Maximum time for high and low
tmax = 26;
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

Juvg = zeros(astepsJ, msteps, 27,5);
Anewg = zeros(28,msteps,5);
Juvgl = zeros(astepsJ, msteps, 21,5);
Anewgl = zeros(22,msteps,5);
%
II = ones(13,6) + [0 0 0 0 0 0 ; -0.2.*eye(6); eye(6).*0.2];

% In order to run the final run, dm needs to increase for CFL reasons.
for kk = 1:12  
    tic
IV1 = 78;
IV2 = 26;
%netpar = [ 8.5 0.65 9e-4 2.5 0.05 0.28].*II(kk,:); % original
%netpar = [  6.8    0.78    0.0011    3   0.04    0.336 ].*II(kk,:);
% netpar = [  6.8    0.24   0.0011    3   0.04    0.336 ].*II(kk,:); b
%netpar = [  6.8    0.28   0.009    3   0.048    0.332 ].*II(kk,:);  c 4-7
%netpar = [8.2  0.28   0.0009    3   0.048    0.332 ].*II(kk,:); c 8-11
%
%netpar = [  8.5000    0.5200    0.0009    3.0000    0.0480    0.3320 ].*II(kk,:);  % d
%netpar = [  8.5000    0.5    0.0006    3.0000    0.05    0.33 ].*II(kk,:);  % e 
netpar = [  10.8    0.6    0.0006    3.0000    0.05    0.33 ].*II(kk,:);  % f 24.8411
%netpar = [  10.8    0.6    0.0006    3.0000    0.05    0.396 ];

nn2 = netpar(1); % mass exponent alpha2
ms = netpar(2); % mass coefficient n2
as = netpar(3); % age coef alpha1
nn= netpar(4); % age exponent n1
rs = netpar(5); % decay rate of resource
g1 = netpar(6); % max growth parameter gamma_j
dj = 0.025;
ys = 8;


gj = @(m, tt) g1.*m.^(2/3).*exp(-tt*rs);


% This is the current and future time step
% age,mass
A1 = zeros(astepsA, msteps);
A2 = zeros(astepsA, msteps);
J2 = zeros(astepsJ, msteps); 
J1 = J2;
A1b = zeros(astepsA, msteps);
A2b = zeros(astepsA, msteps);
J2b = zeros(astepsJ, msteps); 
J1b = J2;


%
% Matrix of (age, mass,time)
Juv = zeros(astepsJ, msteps, ceil(tmax*tm1));
Anew = zeros(28,msteps);
Juvl = zeros(astepsJ, msteps, ceil(tmax2*tm1));
Anewl = zeros(22,msteps);

%
% Initial values 
 J1(1,1) = IV1/(dm*da);
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

Juv(:,:, 1)  = J1;
Juvl(:,:, 1)  = J1b;
        

for i = dt:dt:(tmax2)
     
    for k = 1
        j =1;
         A2(j,k) = A1(j,k) - dt*(A1(j,k)/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
         J2(j,k) = J1(j,k) - dt*(J1(j,k)/da + gj(((k-1).*dm + m0), i)*J1(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
         A2b(j,k) = A1b(j,k) - dt*(A1b(j,k)/da + ga*A1b(j,k)/dm + (d0+ da*d1*(j-1))*A1b(j,k));
         J2b(j,k) = J1b(j,k) - dt*(J1b(j,k)/da + gj(((k-1).*dm + m0), i)*J1b(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
       
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*A1(j,k)/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + gj(((k-1).*dm + m0),i)*J1(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da + ga*A1b(j,k)/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da + gj(((k-1).*dm + m0),i)*J1b(j,k)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        end
    end
    
    for k = 2:(msteps-1)
        j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)- gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
          A2b(j,k) = A1b(j,k) - dt*((A1b(j,k))/da + ga*(A1b(j,k)-A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k))/da + (gj(((k-1).*dm + m0),i)*J1b(j,k)- gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da + ga*(A1(j,k)-A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da + (gj(((k-1).*dm + m0),i)*J1(j,k)-gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da + ga*(A1b(j,k)-A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da + (gj(((k-1).*dm + m0),i)*J1b(j,k)-gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
          
        end

    end
    
    
    % This is only so that I do not loose those extra large individuals. As
    % individuals should not be this large, if this for large, then the
    % parameters are unrealistic.
   k = msteps;
   
         j = 1;
           A2(j,k) = A1(j,k) - dt*((A1(j,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k))/da - gj(((k-2).*dm + m0),i)*J1(j,k-1)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1(j,k));
          A2b(j,k) = A1b(j,k) - dt*((A1b(j,k))/da - ga*(A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k))/da - gj(((k-2).*dm + m0),i)*J1b(j,k-1)/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)) )*J1b(j,k));
        
        for j = 2:astepsJ
           A2(j,k) = A1(j,k) - dt*((A1(j,k)-A1(j-1,k))/da - ga*(A1(j,k-1))/dm + (d0+ da*d1*(j-1))*A1(j,k));
           J2(j,k) = J1(j,k) - dt*((J1(j,k)-J1(j-1,k))/da - (gj(((k-2).*dm + m0),i)*J1(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)))*J1(j,k));
   
           A2b(j,k) = A1b(j,k) - dt*((A1b(j,k)-A1b(j-1,k))/da - ga*(A1b(j,k-1))/dm + (d0+ da*d1*(j-1))*A1b(j,k));
           J2b(j,k) = J1b(j,k) - dt*((J1b(j,k)-J1b(j-1,k))/da - (gj(((k-2).*dm + m0),i)*J1b(j,k-1))/dm + (dj +  (1- exp(-((k-1)*dm).^nn2.*ms)).*(1- exp(-((j-1)*da).^nn.*as)))*J1b(j,k));
        end
    
 
 cc = w'*(J1.*ggam).*dt./3;
 A2(1,:) = A2(1,:)+ cc;  
 ccb = w'*(J1b.*ggam).*dt./3;
 A2b(1,:) = A2b(1,:)+ ccb;  
 
 Anew(ceil(i),:) = Anew(ceil(i),:)+ (cc);
    A1 = A2;
    J1 = J2;
 Anewl(ceil(i),:) = Anewl(ceil(i),:)+ (ccb);
    A1b = A2b;
    J1b = J2b;
    if( sum(sum(A1 <0)) >0)
        i
        break
    elseif(dt/da + gj(maxm,i+2)*dt/dm > 0.994)
        error('CFL issue')
    end


    % no births
% born = beta(A1);
% J1(1,1) = J1(1,1) + born;
% born = beta(A1b);
% J1b(1,1) = J1b(1,1) + born;




%Jnewd(round(i/dt)) = born;

    if(mod(round((i + 0.001),2),tm1) ==0)
        Juv(:,:, round(i)+1)  = J1;
        Juvl(:,:, round(i)+1)  = J1b;
    end
    
end


for i = (tmax2 + dt):dt:tmax
     
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
        %Ad(:,:, round(i)+1)  = A1;
        Juv(:,:, round(i)+1)  = J1;
   end
%         Ad(:,:, round(i/dt)+1)  = A1;
%        Juv(:,:, round(i/dt)+1)  = J1;

end
       
        Anewg(:,:,kk) = Anew;
        Anewgl(:,:,kk) = Anewl;
       
        
     
HighS = zeros(13,17);
LowS = zeros(16,17);
    AA = Anew(:,:).*dm.*da;
    AAl = Anewl(:,:).*dm.*da;
    for  jj = 1:17
        for ii = 1:13
            HighS(ii,jj) = sum(sum(AA((2*ii-1):(2*ii),((jj-1)*15+16):(jj*15+15)) ));
            LowS(ii,jj) = sum(sum(AAl(ii,((jj-1)*15+16):(jj*15+15)) ));
        end
        for ii = 14:16
            LowS(ii,jj) = sum(sum(AAl(ii,((jj-1)*15+16):(jj*15+15)) ));
        end
        
    end

        save('ParCombFitf13.mat',  'Anewg',  'Anewgl')

       er= norm(HighM-HighS, 'fro') + 2.*norm(LowM-LowS, 'fro');
       
      
  toc

    
end



%%


