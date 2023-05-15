function [new_std,new_cycles,new_mean] = combine_2_stds(std1,std2,mean1,mean2,n,m)
    var1 = std1.^2;
    var2 = std2.^2;
    new_mean = (n*mean1+m*mean2)/(n+m);
    new_cycles = n+m;
    combined_var_part1 = ((n-1)*var1+(m-1)*var2)/(n+m-1);
    combined_var_part2 = (n*m*(mean1-mean2).^2)/((n+m)*(n+m-1));
    new_std = sqrt(combined_var_part1 + combined_var_part2);
end