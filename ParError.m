
% This will find the error between the data and simulations shown
% please choose the following value 
% ch= 0 is low density run only
% ch= 1 is high density run only
% ch= 2 is both low and high density run
% Note that if the step sizes are changed then this needs to be changed as
% it is hard coded

%load('ParLowFitg.mat')
% comment out or assign
%ch = 1;
runit = 13;

% get Data
dataFormat

% Make simulation into similar groups as data




if ch == 2
    HighS = zeros(13,17,runit);
    LowS = zeros(16,17,runit);
    for kk = 1:(runit-1)
        AA = Anewg(:,:,kk).*0.02.*0.02;
        AAl = Anewgl(:,:,kk).*0.02.*0.02;
        for  jj = 1:17
            for ii = 1:13
                HighS(ii,jj,kk) = sum(sum(AA((2*ii-1):(2*ii),((jj-1)*15+16):(jj*15+15)) ));
                LowS(ii,jj,kk) = sum(sum(AAl(ii,((jj-1)*15+16):(jj*15+15)) ));
            end
            for ii = 14:16
                LowS(ii,jj,kk) = sum(sum(AAl(ii,((jj-1)*15+16):(jj*15+15)) ));
            end
        
        end
    end



% Make simulation into similar groups as data, special run needed for change in step for cfl issues


    kk = runit;
        AA = Anewg13.*0.02.*0.05;
        AAl = Anewgl13.*0.02.*0.05;
        for  jj = 1:17
            for ii = 1:13
                HighS(ii,jj,kk) = sum(sum(AA((2*ii-1):(2*ii),(jj*6+1):(jj*6+6)) ));
                LowS(ii,jj,kk) = sum(sum(AAl(ii,(jj*6+1):(jj*6+6)) ));
            end
            for ii = 14:16
                LowS(ii,jj,kk) = sum(sum(AAl(ii,(jj*6+1):(jj*6+6)) ));
            end
        
        end
        
    err = zeros(runit,1);
%

        for ii = 1:runit
            err(ii) =  norm(HighM-HighS(:,:,ii), 'fro') + 2.*norm(LowM-LowS(:,:,ii), 'fro');
        end

else
    




    err = zeros(runit,1);
%
    if  ch == 1
        HighS = zeros(13,17,runit);

        for kk = 1:(runit-1)
        AA = Anewg(:,:,kk).*0.02.*0.02;
        for  jj = 1:17
            for ii = 1:13
                HighS(ii,jj,kk) = sum(sum(AA((2*ii-1):(2*ii),((jj-1)*15+16):(jj*15+15)) ));
            end
        
        end
        end
    
        kk = 13;
        AA = Anewg13.*0.05.*0.02;
        for  jj = 1:17
            for ii = 1:13
                HighS(ii,jj,kk) = sum(sum(AA((2*ii-1):(2*ii), (jj*6+1):(jj*6+6) ) ));
            end
        
        end

        for ii = 1:runit
            err(ii) =  norm(HighM-HighS(:,:,ii), 'fro') ;
        end
    
    elseif ch == 0
        LowS = zeros(16,17,runit);
    
        for kk = 1:(runit)
            AAl = Anewgl(:,:,kk).*0.02.*0.05;
            for  jj = 1:17
                for ii = 1:16
                    LowS(ii,jj,kk) = sum(sum(AAl(ii,(jj*6+1):(jj*6+6)) ));
                end
        
            end
        end
    
        for ii = 1:runit
            err(ii) =  norm(LowM-LowS(:,:,ii), 'fro');
        end  
    end
end
  
clearvars ii jj kk AA ch runit
%%


