function [ output_args] = HDMR_function(f,f_0,x,number_of_constructed)
%UNTITLED HDMR_function is a function that sumerize all member function in
%the first order and the const term, dimension_array is the dimension order
%when construct the member functions for the first order terms. 
%   此处显示详细说明
[m p] = size(x);
n = number_of_constructed; % number of total dimension
sum = repmat(f_0,m,1);
for j = 1:m
    for i = 1:n
        sum(j) = sum(j) + f{1,i}.func(x(j,i));
    end
end
output_args = sum;

end

