function [ output_args ] = GetTRFactor( Y_star_k_1,Y_star_k,Y_star_predict_k )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
output_args = (Y_star_k_1 - Y_star_k)/(Y_star_k_1 - Y_star_predict_k);

end

