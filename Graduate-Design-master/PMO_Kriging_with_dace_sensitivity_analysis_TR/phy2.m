function [ output_args ] = phy2( x,X,c)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
[m n] = size(X);%有多少组数
[p q] = size(x);
for j = 1:p
    for i = 1:m
        output_args(j,i) = phy(x(j,:)',X(i,:)',c); 
    end
end
output_args = output_args';
end

