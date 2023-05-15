clear, clc, close all
A = readtable('test.csv');
z = table2array(A(1,1:300));
n = 100;
m = 200;
x1 = table2array(A(1,1:100)); %n
x2 = table2array(A(1,101:300)); %m

mean_z = mean(z);
std_z = std(z);
mean_x1 = mean(x1);
mean_x2 = mean(x2);
std_x1 = std(x1);
std_x2 = std(x2);
var_x1 = std_x1^2;
var_x2 = std_x2^2;

combined_mean = (n*mean_x1+m*mean_x2)/(n+m);
combined_var_part1 = ((n-1)*var_x1+(m-1)*var_x2)/(n+m-1);
combined_var_part2 = (n*m*(mean_x1-mean_x2)^2)/((n+m)*(n+m-1));
combined_std = sqrt(combined_var_part1 + combined_var_part2);

% 0.111160840541184 - prawdziwe
% 0.111160840541184 - wyliczone