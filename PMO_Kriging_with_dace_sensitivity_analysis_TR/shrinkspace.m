function [ output_args ] = shrinkspace(Y_star_k_1,Y_star_k,Y_star_predict_k,...
    X_star_k,X_star_k_1,lower_bound,upper_bound,lower_bound_current,...
    upper_bound_current)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
c1 = 0.75;
c2 = 1.25;
r1 = 0.1;
r2 = 0.75;
lambda = 0.05;
Delta = upper_bound - lower_bound;
if(Y_star_k_1 > Y_star_k)
    X_c = X_star_k;
else
    X_c = X_star_k_1;
end
r = GetTRFactor(Y_star_k_1,Y_star_k,Y_star_predict_k);
delta_k = UpdateDelta(r,X_star_k,X_star_k_1,c1,c2,r1,r2,Delta);
delta_k = max(delta_k,lambda*Delta);
lower_bound_current = X_c - delta_k;
upper_bound_current = X_c + delta_k;
lower_bound_current = max(lower_bound_current,lower_bound);
upper_bound_current = min(upper_bound_current,upper_bound);
output_args = [lower_bound_current;upper_bound_current];
end

