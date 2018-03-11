function [ output_args ] = phy( x1,x2,c )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
output_args = ((x1 - x2)'*(x1 - x2)+c^2)^(-0.5);

end

