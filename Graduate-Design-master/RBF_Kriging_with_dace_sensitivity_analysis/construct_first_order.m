function [ NOE,output_args,sample_point,sample_value] = construct_first_order(j,dth,sample_point,sample_value,x_0,f_0,lower_bound,upper_bound,objective_function,NOE)
%   construct_first_order is a function for constructing the first order
%member function for jth dimension. j is the dimension to be constructed,
%sample_point and sample_value are the database saving the points for
%constructing member function and true values of those points, x_0 and f_0
%are the current const terms
%   此处显示详细说明
sample_point{j} = [lower_bound(dth); upper_bound(dth)]; % choose first two points at two terminals
sample_value{j}= ...
    objective_function(point1DtoND(sample_point{j},x_0,dth)) - repmat(f_0,size(sample_point{j},1),1);
% values at two chosen points
NOE = NOE + 2; % two times function call
x_test = x_0(dth); % x_test is the jth dimension value of x_0
test_value = f_0; % value at test point
for loop = 1:4 % for-loop to find the member function in dth dimension and saved in the jth order of database f...
%   which satisfies the convergent condition
    RBF.model = RBFnouse( sample_point{j}, sample_value{j});
    RBF.predic_func = @(x)RBF_predictor(RBF.model, x);
    detector = detect(RBF.predic_func(x_test),test_value - f_0); % convergent condition
    if detector == 1 % satisfy the convergent condition
        break;
    else % not satisfy the convergent condition
        sample_point{j} = [sample_point{j}; x_test];
        sample_value{j} = [sample_value{j};test_value - repmat(f_0,size(x_test,1),1)];
        x_test = (lower_bound(dth) + rand()*(upper_bound(dth) - lower_bound(dth)));
        test_value = objective_function(point1DtoND(x_test,x_0,dth));
        NOE = NOE + 1;
        continue;
    end     
end  
output_args = RBF;
end

