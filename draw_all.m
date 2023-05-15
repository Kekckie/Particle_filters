function [] = draw_all(object,N)
% Function combines and draws all data
A = readtable(append(object,'/aRMSE.csv'));
B = readtable(append(object,'/Jx.csv'));
C = readtable(append(object,'/Jy.csv'));
D = readtable(append(object,'/eps.csv'));

mean_aRMSE = table2array(A(:,1:7));
aRMSE_odch = table2array(A(:,8:14));
mean_Jx = table2array(B(:,1:7));
Jx_odch = table2array(B(:,8:14));
mean_Jy = table2array(C(:,1:7));
Jy_odch = table2array(C(:,8:14));
mean_eps = table2array(D(:,1:7));
eps_odch = table2array(D(:,8:14));
cycles_aRMSE = table2array(A(:,15));
cycles_Jx = table2array(B(:,15));
cycles_Jy = table2array(C(:,15));
cycles_eps = table2array(D(:,15));

while (size(mean_aRMSE,1)>1)
    [temp_std, temp_cycles, temp_mean] = combine_2_stds(aRMSE_odch(1,:),aRMSE_odch(2,:)...
        ,mean_aRMSE(1,:),mean_aRMSE(2,:),cycles_aRMSE(1),cycles_aRMSE(2));
    mean_aRMSE(2,:) = temp_mean;
    cycles_aRMSE(2) = temp_cycles;
    mean_aRMSE = mean_aRMSE(2:end,:);
    cycles_aRMSE = cycles_aRMSE(2:end);
    aRMSE_odch(2,:) = temp_std;
    aRMSE_odch = aRMSE_odch(2:end,:);
end
while (size(mean_Jx,1)>1)
    [temp_std, temp_cycles, temp_mean] = combine_2_stds(Jx_odch(1,:),Jx_odch(2,:)...
        ,mean_Jx(1,:),mean_Jx(2,:),cycles_Jx(1),cycles_Jx(2));
    mean_Jx(2,:) = temp_mean;
    cycles_Jx(2) = temp_cycles;
    mean_Jx = mean_Jx(2:end,:);
    cycles_Jx = cycles_Jx(2:end);
    Jx_odch(2,:) = temp_std;
    Jx_odch = Jx_odch(2:end,:);
end
while (size(mean_Jy,1)>1)
    [temp_std, temp_cycles, temp_mean] = combine_2_stds(Jy_odch(1,:),Jy_odch(2,:)...
        ,mean_Jy(1,:),mean_Jy(2,:),cycles_Jy(1),cycles_Jy(2));
    mean_Jy(2,:) = temp_mean;
    cycles_Jy(2) = temp_cycles;
    mean_Jy = mean_Jy(2:end,:);
    cycles_Jy = cycles_Jy(2:end);
    Jy_odch(2,:) = temp_std;
    Jy_odch = Jy_odch(2:end,:);
end
while (size(mean_eps,1)>1)
    [temp_std, temp_cycles, temp_mean] = combine_2_stds(eps_odch(1,:),eps_odch(2,:)...
        ,mean_eps(1,:),mean_eps(2,:),cycles_eps(1),cycles_eps(2));
    mean_eps(2,:) = temp_mean;
    cycles_eps(2) = temp_cycles;
    mean_eps = mean_eps(2:end,:);
    cycles_eps = cycles_eps(2:end);
    eps_odch(2,:) = temp_std;
    eps_odch = eps_odch(2:end,:);
end

aRMSE_err(1) = 2.2281 * (aRMSE_odch(1))/sqrt(cycles_aRMSE);
Jx_err(1) = 2.2281 * (Jx_odch(1))/sqrt(cycles_Jx);
Jy_err(1) = 2.2281 * (Jy_odch(1))/sqrt(cycles_Jy);
eps_err(1) = 2.2281 * (eps_odch(1))/sqrt(cycles_eps);
for i=2:size(mean_aRMSE,2)
    aRMSE_err(i) = 1.96 * (aRMSE_odch(i))/sqrt(cycles_aRMSE);
    Jx_err(i) = 1.96 * (Jx_odch(i))/sqrt(cycles_Jx);
    Jy_err(i) = 1.96 * (Jy_odch(i))/sqrt(cycles_Jy);
    eps_err(i) = 1.96 * (eps_odch(i))/sqrt(cycles_eps);
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