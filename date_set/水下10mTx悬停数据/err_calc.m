% Code to calculate the MSE of a given vector x with respect to a target y
function [mse, me] = err_calc(x,y)
    
    n = length(x);
    
    sum_mse = 0;
    sum_me = 0;
    
    for i = 1:n
        sum_mse = sum_mse + (x(i) - y(i))^2;
        sum_me = sum_me + abs(x(i) - y(i));
    end
    
    mse = sqrt(sum_mse / n);
    me = sum_me / n;

end