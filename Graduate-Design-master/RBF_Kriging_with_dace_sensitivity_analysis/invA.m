function [ output_args ] = invA(x,c)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
[m n] = size(x);
A = ones(m,m);
for i = 1:m
    for j = 1:m
        A(i,j) = phy(x(i,:)',x(j,:)',c);
    end
end
output_args = inv(A);
end

