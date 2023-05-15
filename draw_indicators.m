function [] = draw_indicators(object,index,N)
% Function draws one row given in 'index' from 'object'.
A = readtable(append(object,'/aRMSE.csv'));
B = readtable(append(object,'/Jx.csv'));
C = readtable(append(object,'/Jy.csv'));
D = readtable(append(object,'/eps.csv'));

mean_aRMSE = table2array(A(index,1:7));
aRMSE_odch = table2array(A(index,8:14));
mean_Jx = table2array(B(index,1:7));
Jx_odch = table2array(B(index,8:14));
mean_Jy = table2array(C(index,1:7));
Jy_odch = table2array(C(index,8:14));
mean_eps = table2array(D(index,1:7));
eps_odch = table2array(D(index,8:14));
cycles = table2array(A(index,15));


aRMSE_err(1) = 2.2281 * (aRMSE_odch(1))/sqrt(cycles);
Jx_err(1) = 2.2281 * (Jx_odch(1))/sqrt(cycles);
Jy_err(1) = 2.2281 * (Jy_odch(1))/sqrt(cycles);
eps_err(1) = 2.2281 * (eps_odch(1))/sqrt(cycles);
for i=2:size(mean_aRMSE,2)
    aRMSE_err(i) = 1.96 * (aRMSE_odch(i))/sqrt(cycles);
    Jx_err(i) = 1.96 * (Jx_odch(i))/sqrt(cycles);
    Jy_err(i) = 1.96 * (Jy_odch(i))/sqrt(cycles);
    eps_err(i) = 1.96 * (eps_odch(i))/sqrt(cycles);
end

subplot(2,2,1)
semilogx(N,mean_aRMSE, '-*')
hold on
errorbar(N,mean_aRMSE,aRMSE_err, 'Marker','none','LineStyle', 'none')
title('aRMSE')


subplot(2,2,2)
semilogx(N,mean_Jx, '-*')
hold on
errorbar(N,mean_Jx,Jx_err, 'Marker','none','LineStyle', 'none')
title('Jx')

subplot(2,2,3)
semilogx(N,mean_Jy, '-*')
hold on
errorbar(N,mean_Jy,Jy_err, 'Marker','none','LineStyle', 'none')
title('Jy')

subplot(2,2,4)
semilogx(N,mean_eps, '-*')
hold on
errorbar(N,mean_eps,eps_err, 'Marker','none','LineStyle', 'none')
title('eps')
end