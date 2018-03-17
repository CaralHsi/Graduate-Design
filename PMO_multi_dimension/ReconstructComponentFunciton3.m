function [ output_args, sample_point, sample_value] = ReconstructComponentFunciton3(sample_point, sample_value, new_sample, new_value, f_0_old, f_0_new, second_order_count)
sample_point{2,second_order_count} = [sample_point{2,second_order_count}; new_sample];
sample_value{2,second_order_count} = [sample_value{2,second_order_count}; new_value];
sample_value{2,second_order_count} = sample_value{2,second_order_count} - repmat(f_0_new - f_0_old,length(sample_value{2,second_order_count}),1);
[sample_point{2,second_order_count},ia,~] = unique(sample_point{2,second_order_count},'rows');
sample_value{2,second_order_count} = sample_value{2,second_order_count}(ia);
ll = size(sample_value{2, second_order_count}, 2);
theta = repmat(10, 1, ll); lob = repmat(1e-1, 1, ll); upb = repmat(20, 1, ll);
[dmodel, perf] = dacefit(sample_point{2,second_order_count}, sample_value{2,second_order_count}, @regpoly0, @corrgauss, theta, lob, upb);
pre.dmodel = dmodel;
pre.perf = perf;
pre.func = @(x)predictor(x, pre.dmodel);
pre.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
output_args = pre;
 
 
end
