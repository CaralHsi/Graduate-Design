function [NOE,output_args] = construct_second_order( f,choose,sample_point,sample_value,x_0,f_0,objective_function,NOE,K)
%   construct_second_order is a function that calculate the correlation
%   function  between the ath and bth first-order member functions.
%   function to be biult will be saved in number second_order_count
%   functions in f{2}
%   此处显示详细说明
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
        temp_value(kk) = f{1,dimension_concerned(k)}.predic_func(temp(kk));
    end
    sample_value_temp = sample_value_temp - temp_value;
end
for loop = 1:8
    RBF.model = RBFnouse( sample_point_temp, sample_value_temp);
    RBF.predic_func = @(x)RBF_predictor(RBF.model, x);
    test_value = objective_function(point1DtoND(test_point,x_0,[K(dimension_concerned(1)) K(dimension_concerned(2))])) - f_0;
    NOE = NOE + 1;
    for k = 1:2
        temp = test_point(k);
        test_value = test_value - f{1,dimension_concerned(k)}.predic_func(temp);
    end
    detector = detect(RBF.predic_func(test_point),test_value);
    if detector == 1
        break;
    else
        sample_point_temp = [sample_point_temp; test_point];
        sample_value_temp = [sample_value_temp;test_value];
        num = num + 1;
        if num == length(test_point_lib)
            break;
        else
            test_point = test_point_lib(num,:);
            continue;
        end
    end     
end

output_args = RBF;
end