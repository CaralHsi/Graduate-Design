function [NOE,output_args, sample_point_temp, sample_value_temp] = construct_second_order( HDMR, f_append, sample_point_whole,sample_value,x_0,f_0,objective_function,NOE,K,lower_bound, upper_bound, second_order_count)
%   construct_second_order is a function that calculate the correlation
%   function  between the ath and bth first-order member functions.
%   function to be biult will be saved in number second_order_count
%   functions in f{2}
%   此处显示详细说明
%% set ga options
ga_opts = gaoptimset('Display','off');

%% sample point for correlation function
%   j = second_order_count; % the sequence number of where output_args will be saved
sample = sample_point_whole{1,1};
for i = 2:(second_order_count + 1)  
    sample = [sample; sample_point_whole{1,i}];
end
sample_point_temp = sample(:,K(1:second_order_count + 1));

%% test point is ramdomly selected
t = rand(1, second_order_count + 1);
test_point = t.*(upper_bound(K(1:second_order_count + 1)) - lower_bound(K(1:second_order_count + 1))) + lower_bound(K(1:second_order_count + 1));

%% construct RBF function using sample_point_temp and sample_value_temp
l = size(sample_point_temp,1);
sample_value_temp = zeros(l,1);
for kk = 1:l
    sample_value_temp(kk) = objective_function(sample(kk,:)) - HDMR(sample_point_temp(kk,1:end - 1)) - f_append.func(sample_point_temp(kk,end));
end
% sample_point_temp = sample_point_temp(:,choose);
for loop = 1:8*(second_order_count)
    [sample_point_temp,ia,~] = unique(sample_point_temp,'rows');
    sample_value_temp = sample_value_temp(ia);
    [dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, ones(1, second_order_count + 1),repmat(1e-3, 1, second_order_count + 1), repmat(100, 1, second_order_count + 1));
    pre.dmodel = dmodel;
    pre.perf = perf;
    pre.func = @(x)predictor(x, pre.dmodel);
    pre.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
    test_value = objective_function(point1DtoND(test_point,x_0,K(1:second_order_count + 1))) - HDMR(test_point(1:end - 1)) - f_append.func(test_point(end));
    NOE = NOE + 1;
    detector = detect(pre.func(test_point),test_value);
    if loop == 8*(second_order_count) || detector == 1 
        break;
    else
        sample_point_temp = [sample_point_temp; test_point];
        sample_value_temp = [sample_value_temp;test_value];
        test_point = ga(pre.mse,second_order_count + 1,[],[],[],[],lower_bound(K(1:second_order_count + 1)),upper_bound(K(1:second_order_count + 1)), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
    end     
end
output_args = pre;
end