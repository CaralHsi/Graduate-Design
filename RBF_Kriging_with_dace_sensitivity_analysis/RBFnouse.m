function [output_args] = RBFnouse( sample_point, sample_value)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
[sample_point,ia] = unique(sample_point,'rows');
sample_value = sample_value(ia);
[m,n] = size(sample_point);
A = zeros(m,m);
for i = 1:m
    for j = 1:m
        A(i,j) = (sample_point(i,:) - sample_point(j,:))*(sample_point(i,:) - sample_point(j,:))'...
            * log(norm(sample_point(i,:) - sample_point(j,:)));
         if norm(sample_point(i,:) - sample_point(j,:)) == 0
            A(i,j) = 0;
         end
    end
end
P = cell(n+1,1);
P{1} = @(x)1;
for i = 1:n
    P{i+1} = @(x)x(i);
end
P_bar = zeros(m,n+1);
for i = 1:m
    for j = 1:n+1
        P_bar(i,j) = P{j}(sample_point(i,:));
    end
end
invA = ([A P_bar;P_bar' zeros(n+1,n+1)]);
beta = (invA)\[sample_value;zeros(n+1,1)];
alpha = beta(m+1:m+n+1);
beta = beta(1:m);
output_args.P = P;
output_args.sample_point = sample_point;
output_args. beta = beta;
output_args.alpha = alpha;
output_args.m = m;
output_args.n = n;
end

