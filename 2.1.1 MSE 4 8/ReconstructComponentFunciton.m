function [ output_args, sample_point, sample_value] = ReconstructComponentFunciton( sample_point, sample_value, new_sample, new_value)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
sample_point = [sample_point; new_sample];
sample_value = [sample_value; new_value];
sample_value = sample_value - repmat(new_value,length(sample_value),1);
[sample_point,ia,ic] = unique(sample_point,'rows');
sample_value = sample_value(ia);
theta = [10]; lob = [1e-1]; upb = [20];
[dmodel, perf] = dacefit(sample_point, sample_value, @regpoly0, @corrgauss, theta, lob, upb);
pre.dmodel = dmodel;
pre.perf = perf;
pre.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
pre.func = @(x)predictor(x,pre.dmodel);
output_args = pre;


end

