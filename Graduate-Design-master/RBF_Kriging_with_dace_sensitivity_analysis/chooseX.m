function [ S ] = chooseX( LB,UB,N,D )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

S = lhsdesign(N,D,'criterion','maximin');
S = S.*repmat(UB - LB,N,1) + repmat(LB,N,1);

end

