function [NOE,output_args, sample_point_temp, sample_value_temp] = construct_second_order_substitute2( HDMR, f_append, sample_point_whole,sample_value,x_0,f_0,objective_function,NOE,K,lower_bound, upper_bound, second_order_count, sample_point, x_0_old)
%   construct_second_order is a function that calculate the correlation
%   function  between the ath and bth first-order member functions.
%   function to be biult will be saved in number second_order_count
%   functions in f{2}
%   此处显示详细说明
%% set ga options
ga_opts = gaoptimset('Display','off');

%% construct cut planes
d = second_order_count + 1;
sample_point_temp = [];
sample_point_value = [];
prefunc = cell(1, d);
prefunctotal = @(x)0;
if second_order_count > 1 % for 2 dimensional problem, this term is not necessary to do
    for i = 1:second_order_count % for 3 dimensional problem, second_order_count = 2, ...
    %and we need to construct 13 and 23 plane, i.e. 2 planes in total. for 4
    %dimensional problem, second_order_count = 3, and we need to construct 14 24
    %34 plane, i.e. 3 planes in total. For n dimensional problem,
    %second_order_count = d-1, and we need to construct 1 d, 2 d, ..., d - 1 d
    %plane, i.e. d - 1 plane in total.
        sam_t = point1DtoND(sample_point{1, d}, x_0(K([i, d])), 2);
        l = size(sam_t, 1);
        val_t = zeros(l,1);
        for kk = 1:l
            val_t(kk) = objective_function(point1DtoND(sam_t(kk, :), x_0, K([i, d]))) - HDMR(point1DtoND(sam_t(kk, :), x_0(K(1:d - 1)), K(i))) - f_append.func(sam_t(kk, 2));
        end
    %     val_t = zeros(l, 1); % after verification, use this line instead of
    %     the above 3 line.
        t = rand(1, 2);
        test_t = t.*(upper_bound(K([i, d])) - lower_bound(K([i, d]))) + lower_bound(K([i, d]));
        detector = 0;
        count = 0;
        while (count < 8)
            count = count + 1;
            [sam_t,ia,~] = unique(sam_t,'rows');
            val_t = val_t(ia);
            [dmodel, perf] = dacefit(sam_t, val_t, @regpoly0, @corrgauss, ones(1, 2),repmat(1e-3, 1, 2), repmat(100, 1, 2));
            prefunc{i}.dmodel = dmodel;
            prefunc{i}.perf = perf;
            prefunc{i}.func = @(x)predictor(x, prefunc{i}.dmodel);
            prefunc{i}.mse = @(x)-Get_MSE(@predictor, prefunc{i}.dmodel, x);
            test_v_t = objective_function(point1DtoND(test_t,x_0,K([i, d]))) - HDMR(point1DtoND(test_t,x_0(K(1: d - 1)),K(i))) - f_append.func(test_t(2));
            NOE = NOE + 1;
            detector = detect(prefunc{i}.func(test_t),test_v_t);
%             detector = detect2(prefunc{i}.func(test_t),test_v_t, objective_function(point1DtoND(test_t,x_0,K([i, d]))));
%             if detector == 1
%                 count = count + 1;
%             else
%                 count = 0;
%             end
            sam_t = [sam_t; test_t];
            val_t = [val_t; test_v_t];
            test_t = ga(prefunc{i}.mse, 2, [], [], [], [], lower_bound(K([i, d])),upper_bound(K([i, d])), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
        end
        prefunctotal  = @(x)prefunctotal(x) + prefunc{i}.func(x(K([i, d])));
        sample_point_temp = [sample_point_temp; point1DtoND(sam_t, x_0, K([i, d]))];
        sample_point_value = [sample_point_value; val_t];    
    end
end

%% sample point for correlation function
%   j = second_order_count; % the sequence number of where output_args will be saved
sample = sample_point_whole{1,1};
for i = 2:(second_order_count + 1)  
    sample = [sample; sample_point_whole{1,i}];
end
for i = 1:second_order_count - 1
   sample =  [sample; sample_point_whole{2,i}];
end
sample = [sample_point_temp; sample];
sample_point_temp = sample(:,K(1:second_order_count + 1));

%% test point is ramdomly selected
t = rand(1, second_order_count + 1);
test_point = t.*(upper_bound(K(1:second_order_count + 1)) - lower_bound(K(1:second_order_count + 1))) + lower_bound(K(1:second_order_count + 1));

%% construct RBF function using sample_point_temp and sample_value_temp
l = size(sample_point_temp,1);
sample_value_temp = zeros(l,1);
for kk = 1:l
    sample_value_temp(kk) = objective_function(sample(kk,:)) - prefunctotal(sample(kk,1:end)) - HDMR(sample_point_temp(kk,1:end - 1)) - f_append.func(sample_point_temp(kk,end));
end
count = 0;
% sample_point_temp = sample_point_temp(:,choose);
for loop = 1:8%8*(second_order_count)
    [sample_point_temp,ia,~] = unique(sample_point_temp,'rows');
    sample_value_temp = sample_value_temp(ia);
    [dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, ones(1, second_order_count + 1),repmat(1e-3, 1, second_order_count + 1), repmat(100, 1, second_order_count + 1));
    pre.dmodel = dmodel;
    pre.perf = perf;
    pre.func = @(x)predictor(x, pre.dmodel);
    pre.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
    test_value = objective_function(point1DtoND(test_point,x_0,K(1:second_order_count + 1))) - prefunctotal(point1DtoND(test_point,x_0,K(1:second_order_count + 1))) - HDMR(test_point(1:end - 1)) - f_append.func(test_point(end));
    NOE = NOE + 1;
    detector = detect(pre.func(test_point),test_value);
%     detector = detect2(pre.func(test_point),test_value, objective_function(point1DtoND(test_point,x_0,K(1:second_order_count + 1))) );
    if detector == 1 %|| loop ==  % satisfy the convergent condition
        count = count + 1;
    else
        count = 0;
    end
    if count > 0
        break;
    end
    sample_point_temp = [sample_point_temp; test_point];
    sample_value_temp = [sample_value_temp;test_value];
    test_point = ga(pre.mse,second_order_count + 1,[],[],[],[],lower_bound(K(1:second_order_count + 1)),upper_bound(K(1:second_order_count + 1)), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
end
%% the true sample_point_temp and its correspond value
l = size(sample_point_temp);
for i = 1:l
    sample_value_temp(i) = sample_value_temp(i) + prefunctotal(point1DtoND(sample_point_temp(i,:),x_0,K(1:second_order_count + 1)));
end
[sample_point_temp,ia,~] = unique(sample_point_temp,'rows');
sample_value_temp = sample_value_temp(ia);
[dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, ones(1, second_order_count + 1),repmat(1e-3, 1, second_order_count + 1), repmat(100, 1, second_order_count + 1));
pre.dmodel = dmodel;
pre.perf = perf;
pre.func = @(x)predictor(x, pre.dmodel);

%% return term
loop

output_args = pre;
end