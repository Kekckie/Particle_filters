% Vd = [10]; %standard deviation of the state variables noise
% Vn = [1]; %Standard deviation of the output system noise
% f=@(x)[x(1),x(2)]; %state variables function
% h=@(x)[x(1)+x(2),x(1)-x(2),x(1)*x(2)]; %output function
% f=@(x)[0.5*x(1)+(25*x(1))/(1+x(1)^2)+8*cos(1.2*x(1))];
% h=@(x)[x(1)^2/20];
% x_pocz = [0.1]; %initial values of state variables

%%
clear, clc, close all
load('ob1.mat'); % load object
M = 40; %number of iterations
N = [10,30,50,100,300,500,1000]; %number of particles
resampling_type = 'stratified';
location = append('ob1/',resampling_type); %dont forget to change when chaning object
cycles = 50;

for i=1:cycles
    [aRMSE(i,:), Jx(i,:), Jy(i,:), eps(i,:)] = run_particle_filter(M,Vd,Vn,N,f,h,x_pocz,resampling_type);   
end

for i=1:size(aRMSE,2)
    mean_aRMSE(i) = mean(aRMSE(:,i));
    aRMSE_odch(i) = std(aRMSE(:,i));
    mean_Jx(i) = mean(Jx(:,i));
    Jx_odch(i) = std(Jx(:,i));
    mean_Jy(i) = mean(Jy(:,i));
    Jy_odch(i) = std(Jy(:,i));
    mean_eps(i) = mean(eps(:,i));
    eps_odch(i) = std(eps(:,i));
end

aRMSE_err(1) = 2.2281 * (aRMSE_odch(1))/sqrt(cycles);
Jx_err(1) = 2.2281 * (Jx_odch(1))/sqrt(cycles);
Jy_err(1) = 2.2281 * (Jy_odch(1))/sqrt(cycles);
eps_err(1) = 2.2281 * (eps_odch(1))/sqrt(cycles);
for i=2:size(aRMSE,2)
    aRMSE_err(i) = 1.96 * (aRMSE_odch(i))/sqrt(cycles);
    Jx_err(i) = 1.96 * (Jx_odch(i))/sqrt(cycles);
    Jy_err(i) = 1.96 * (Jy_odch(i))/sqrt(cycles);
    eps_err(i) = 1.96 * (eps_odch(i))/sqrt(cycles);
end

save_to_file(location,mean_aRMSE,aRMSE_odch,mean_Jx,Jx_odch,mean_Jy,Jy_odch,...
    mean_eps,eps_odch,cycles);
% draw_all('ob1/residual',N);