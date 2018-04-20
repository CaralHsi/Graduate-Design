function [ NOE,output_args,sample_point,sample_value] = construct_first_order(j,dth,sample_point,sample_value,x_0,f_0,lower_bound,upper_bound,objective_function,NOE)
%   construct_first_order is a function for constructing the first order
%member function for jth dimension. j is the dimension to be constructed,
%sample_point and sample_value are the database saving the points for
%constructing member function and true values of those points, x_0 and f_0
%are the current const terms
%   此处显示详细说明
%% set ga options
ga_opts = gaoptimset('Display','off');

%% sample point for correlation function
sample_point{1,j} = [lower_bound(dth) + eps; upper_bound(dth) - eps]; % choose first two points at two terminals
sample_value{1,j}= ...
    objective_function(point1DtoND(sample_point{1,j},x_0,dth)) - repmat(f_0,size(sample_point{1,j},1),1);
% values at two chosen points
NOE = NOE + 1; % two times function call
x_test = x_0(dth); % x_test is the jth dimension value of x_0
test_value = objective_function(point1DtoND(x_test,x_0,dth)); % value at test point
count = 0;
min_temp = zeros(1,1000);
for loop = 1:10000000%6 % for-loop to find the member function in dth dimension and saved in the jth order of database f...
%   which satisfies the convergent condition
    [sample_point{1,j},ia] = unique(sample_point{1,j},'rows');
    sample_value{1,j} = sample_value{1,j}(ia);
    theta = 10; lob = 1e-1; upb = 20;
    [dmodel, perf] = dacefit(sample_point{1,j}, sample_value{1,j}, @regpoly0, @corrgauss, theta, lob, upb);
    pre.dmodel = dmodel;
    pre.perf = perf;
    pre.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
    pre.func = @(x)predictor(x,pre.dmodel);
%     count = count + 1;
%     [~, min_temp(count)] = ga(pre.func, 1, [], [], [], [], lower_bound(dth), upper_bound(dth), [], ga_opts);
%     if count > 3
%         if abs((min_temp(count) - min_temp(count - 1)) / (min_temp(count) - min_temp(1))) < 0.02 && abs((min_temp(count - 1) - min_temp(count - 2)) / (min_temp(count) - min_temp(1))) < 0.02 %&& Max_value_EI(count - 1) / Max_value_EI(1) < 0.1%std(Max_value_EI(end - 2:end)) / std(Max_value_EI(1:6)) < 0.1
%             break;
%         end
%     end
        detector = detect(pre.func(x_test),test_value - f_0); % convergent condition
        if detector == 1 %|| loop ==  6 
           % break;% satisfy the convergent condition
            count = count + 1;
        else
            count = 0;
        end
 % not satisfy the convergent condition
        if count > 1
            break;
        end
%         else
            sample_point{1,j} = [sample_point{1,j}; x_test];
            sample_value{1,j} = [sample_value{1,j};test_value - repmat(f_0,size(x_test,1),1)];
            x_test = (lower_bound(dth) + rand()*(upper_bound(dth) - lower_bound(dth)));
            x_test = ga(pre.mse,1,[],[],[],[],lower_bound(dth),upper_bound(dth), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
            test_value = objective_function(point1DtoND(x_test,x_0,dth));
            NOE = NOE + 1;
%         end
end  
output_args = pre;
end

