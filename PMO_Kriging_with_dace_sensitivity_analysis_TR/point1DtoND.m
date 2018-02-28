function [ output_args ] = point1DtoND(sample_point_i_j,x_0,choose_i_j)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
[m n] = size(sample_point_i_j);
dimension_concerned = choose_i_j;
size_d = size(dimension_concerned,2);
output_args = repmat(x_0,m,1);
for i = 1:size_d
    output_args(:,dimension_concerned(i)) = sample_point_i_j(:,i);
end

end

