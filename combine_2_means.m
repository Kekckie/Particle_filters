function [new_mean,new_cycles] = combine_2_means(mean1,mean2,n,m)
    new_mean = (n*mean1+m*mean2)/(n+m);
    new_cycles = n+m;
end