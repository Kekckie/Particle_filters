function [] = draw_state_variable(i,M)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
plot(1:M,x(i,:))
hold on
plot(1:M,x_hat(i,:))
end

