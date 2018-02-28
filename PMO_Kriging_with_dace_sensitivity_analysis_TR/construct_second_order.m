function [NOE,output_args] = construct_second_order( f,choose,sample_point,sample_value,x_0,f_0,objective_function,NOE,K,lower_bound, upper_bound)
%   construct_second_order is a function that calculate the correlation
%   function  between the ath and bth first-order member functions.
%   function to be biult will be saved in number second_order_count
%   functions in f{2}
%   此处显示详细说明
%% set ga options
ga_opts = gaoptimset('Display','off');

%% sample point for correlation function
%   j = second_order_count; % the sequence number of where output_args will be saved
dimension_concerned = choose; % two dimensions to be conbined
sample_point_temp = point1DtoND(sample_point{dimension_concerned(1)},x_0,K(dimension_concerned(1)));
sample_point_temp = [sample_point_temp;point1DtoND(sample_point{dimension_concerned(2)},x_0,K(dimension_concerned(2)))];

%% biuld test library
x1 = sample_point{dimension_concerned(1)};
x2 = sample_point{dimension_concerned(2)};
test_point_lib = zeros(length(x1)*length(x2),2);
for k1 = 1:length(x1)
    for k2 = 1:length(x2)
        test_point_lib((k1-1)*length(x2) + k2,:) = [x1(k1),x2(k2)];
    end
end
randIndex = randperm(size(test_point_lib,1));
test_point_lib = test_point_lib(randIndex,:);
num = 1; % first test point is the first row in test library
test_point = test_point_lib(num,:);

%% construct RBF function using sample_point_temp and sample_value_temp
sample = sample_point_temp;sample_point_temp = sample(:,[K(dimension_concerned(1)) K(dimension_concerned(2))]);
sample_value_temp =  objective_function(sample) - repmat(f_0,size(sample,1),1);
for k = 1:2
    temp = sample_point_temp(:,k);
    l = length(temp);
    temp_value = zeros(l,1);
    for kk = 1:l
        temp_value(kk) = f{1,dimension_concerned(k)}.func(temp(kk));
    end
    sample_value_temp = sample_value_temp - temp_value;
end
for loop = 1:8
    [sample_point_temp,ia,ic] = unique(sample_point_temp,'rows');
    sample_value_temp = sample_value_temp(ia);
    [dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, [1 1], [1e-3 1e-3], [100 100]);
    pre.dmodel = dmodel;
    pre.perf = perf;
    pre.func = @(x)predictor(x, pre.dmodel);
    pre.mse = @(x)Get_MSE(@predictor, pre.dmodel, x);
    test_value = objective_function(point1DtoND(test_point,x_0,[K(dimension_concerned(1)) K(dimension_concerned(2))])) - f_0;
    NOE = NOE + 1;
    for k = 1:2
        temp = test_point(k);
        test_value = test_value - f{1,dimension_concerned(k)}.func(temp);
    end
    detector = detect(pre.func(test_point),test_value);
    if detector == 1
        break;
    else
        sample_point_temp = [sample_point_temp; test_point];
        sample_value_temp = [sample_value_temp;test_value];
        num = num + 1;
        if num == length(test_point_lib)
            break;
        else
            test_point = ga(pre.mse,2,[],[],[],[],lower_bound(K(choose)),upper_bound(K(choose)), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
            %test_point = test_point_lib(num,:);
            continue;
        end
    end     
end

output_args = pre;
end