function [ output_args ] = RBF_predictor(RBF_model, x)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
output_args = 0;
for i = 1:RBF_model.m
    temp = norm(x - RBF_model.sample_point(i,:));
    if temp <= 0.0001
        continue;
    else
        output_args = output_args + ...
        RBF_model.beta(i)*(temp).^2*log(norm(x - RBF_model.sample_point(i,:)));

    end
end
for i = 1:RBF_model.n+1
    output_args = output_args + RBF_model.alpha(i)*RBF_model.P{i}(x);
end

end

