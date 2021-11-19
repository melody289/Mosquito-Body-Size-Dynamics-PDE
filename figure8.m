
maxAa = 31; da = 0.02;
%
a2 = 0:da:maxAa ;  astepsA = length(a2);
Tpoint = 1:15:151;
scf = 100;

load('C1par_cyclicMass.mat')
%
sp = zeros(11,1);
ii = 1;
subplot(3,2,ii)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end
xticks(-1)
ylabel('\omega = 7')
set(gca, 'FontSize', 16, 'Position', [0.1    0.7090    0.3700    0.2500])
title('Cyclic Resource in Growth')
plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
axis([ 1    12     0     30])
%
ii = 2;
subplot(3,2,3)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end

plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
xticks(-1)
set(gca, 'FontSize', 16, 'Position', [0.1    0.4    0.3700    0.2500])
axis([ 1    12     0     30])
ylabel('\omega = 19')

ii = 3;
subplot(3,2,5)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end

plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
xticks(1:11)

set(gca, 'FontSize', 16, 'Position', [0.1    0.1    0.3700    0.2500])
set(gcf,  'Position', [ 55    94   945   697])
axis([ 1    12     0     30])
xticklabels(0:15:150)
xlabel('Time')
ylabel('\omega = 27')

%
load('C1par_cyclicDeath.mat')
sp = zeros(11,1);
ii = 1;
subplot(3,2,ii.*2)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end
xticks(-1)

set(gca, 'FontSize', 16, 'Position', [0.5650    0.7090    0.3700    0.2500])
title('Cyclic Resource in Death')
plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
axis([ 1    12     0     30])

ii = 2;
subplot(3,2,ii.*2)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end

plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
xticks(-1)
set(gca, 'FontSize', 16, 'Position', [0.5650    0.4    0.3700    0.2500])
axis([ 1    12     0     30])

ii = 3;
subplot(3,2,ii.*2)
Adm = Adg(1:(end-1),:,:,ii);
Adm = squeeze(sum(Adm,2))./squeeze(sum(sum(Adm)))';
for ij = 1:11
    Adm2 = Adm(:,Tpoint(ij));
 fill( [0 ;Adm2(1); Adm2; 0].*scf + ij , [0, 0, a2(1:end-1),  a2(end-1)], [ 0.4660    0.6740    0.1880])
    hold on
    sp(ij) = find(cumsum(Adm2) > 0.5, 1);
    
end

plot(1:12, [a2((sp )),  a2((sp(end) ))], '--k')
xticks(1:11)

set(gca, 'FontSize', 16, 'Position', [0.5650    0.1    0.3700    0.2500])
set(gcf,  'Position', [ 55    94   945   697])
axis([ 1    12     0     30])
xticklabels(0:15:150)
xlabel('Time')
%%