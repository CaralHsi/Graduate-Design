function [ output_args ] = constructR( theta, D_square, m)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
index = 1;
R = zeros(m);
for i = 1:m
    for j = 1:m
        if i == j
            R(i,j) = 1 + 1*eps;
         elseif i < j
             D_square_ij = D_square(index,:);
             index = index + 1;
            R(i,j) = exp(-theta*D_square_ij');
        else
            R(i,j) = R(j,i);
        end
    end
end
output_args = R;

end

