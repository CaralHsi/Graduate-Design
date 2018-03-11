function [ output_args ] = detect( new_point, old_point )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if abs((new_point - old_point)/(old_point+eps)) <= 0.001
    output_args = 1;
else
    output_args = 0;
end
end

