function [ output_args ] = Get_MSE( predictor, dmodel, x )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
[~,~,c,~] = predictor(x,dmodel);
output_args = c;

end

