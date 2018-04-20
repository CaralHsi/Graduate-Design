%% clean and clear
clc;
clear;  
%% 
syms x;
NOE = 0; % initialize NOE
dimension = 20; % the number of concerned dimension
dimension_array = 1:dimension;
lower_bound = repmat(-4,1,dimension);
upper_bound = repmat(5,1,dimension); % the feasible zone of variables
lower_bound_current = lower_bound;
upper_bound_current = upper_bound;
objective_function = @(x)loadfunc(15,x); % load the objective function
% x_0 = repmat(0,1,dimension); % initialize const term x_0
% set 5 initial points in total and find the best one.
temp1 = 0.18;temp2 = 5;
x_0 = lhsdesign(temp2,dimension);
x_0 = repmat(0.5 - temp1/2, temp2, dimension) + x_0.*temp1;
x_0 = (x_0.*repmat(upper_bound_current(1:dimension) - lower_bound_current(1:dimension),temp2,1) + repmat(lower_bound_current(1:dimension),temp2,1));
f_0 = objective_function(x_0);% calculate f_0 at x_0
[f_0, ia] = min(f_0)
x_0 = x_0(ia,:);
NOE = NOE + 1;
sample_point = cell(1,dimension); % initialize sample_point
sample_point_whole = cell(1, dimension); % initialize sample_point_whole(for all dimension)
sample_value = cell(1,dimension); % initialize sample_value
f = cell(2,1); % initialize member functions
ga_opts = gaoptimset('Display','off');
for j = 1:dimension
    [NOE,f{1,j},sample_point,sample_value] = construct_first_order(j,j,sample_point,sample_value,x_0,f_0,lower_bound_current,upper_bound_current,objective_function,NOE);
end
%     if i == 1
%         test_high_order_point = zeros(dimension,1)';
%         for j = 1:nchoose(i)
%             Xtemp = sample_point{i}{j};
%             Xtemp = Xtemp(randperm(numel(Xtemp)));
%             test_high_order_point(j) = Xtemp(1);
%         end
%         RBFtemp = f_0;
%         for j = 1:nchoose(i)
%             RBFtemp = RBFtemp + f{i,j}(test_high_order_point(j));
%         end
%         detector = detect(RBFtemp,objective_function(test_high_order_point));
%         if detector == 1
%             break;
%         else
%             continue;
%         end
%     end

%% function constructed in first-order optimazation
HDMR = @(x)HDMR_function(f,f_0,x,dimension);
K = 1:dimension;
dimension_number = dimension;
x_0_new = x_0;
TEST_ACCURACY
x_0_old = x_0; f_0_old = f_0;
[x_0_new, f_0_predict] = ga(HDMR,dimension,[],[],[],[],lower_bound',upper_bound',[], ga_opts); % optimization
f_0_new = objective_function(x_0_new); % f_0_new
NOE = NOE + 1;
for i = 1:dimension
    [f_0_new, ia] = min([sample_value{1,i} + repmat(f_0, length(sample_value{1,i}), 1); f_0_new]);
    if(ia <= length(sample_value{1,i}))
        x_0_new = sample_point{1,i}(ia,:);
        x_0_new = point1DtoND(x_0_new, x_0, i);
    end
end
if (f_0_new > f_0_old)
    f_0_new = f_0_old;
    x_0_new = x_0_old;
end
f_0_new


%% roulette wheel selection
K =  zeros(1,dimension);
% selected_dimension = select_dimension_func(); % select the first dimension to construct member function
[S0,V] = sensitivityAnalyse(dimension, HDMR, lower_bound, upper_bound, K);
selected_dimension = roullete(S0);
% selected_dimension = 3;
K(1) = selected_dimension
S0(K(1)) = 0;
newbound = shrinkspace(f_0_old,f_0_new,f_0_predict,x_0_new,x_0_old,...
    lower_bound,upper_bound,lower_bound_current,upper_bound_current);
lower_bound_current = newbound(1,:);
upper_bound_current = newbound(2,:);

%% reset
sample_point = cell(1,dimension); sample_value = cell(1,dimension); f = cell(2,1); 
% reset sample point database and member functions
% x_0_old = x_0; f_0_old = f_0; % update old x_0 and f_0

%% first-loop roulette wheel member function construction
i = 1;j = 1; % first order && first dimension
x_0 = x_0_new; f_0 = f_0_new; % renew the zeroth order point and value
[NOE,f{i,j},sample_point,sample_value] = construct_first_order(j,K(1),sample_point,sample_value,x_0,f_0,lower_bound_current,upper_bound_current,objective_function,NOE);
% construct first-order member functions of first dimension
sample_point_whole{i, j} = point1DtoND(sample_point{i,j},x_0,K(1));

%% function constructed in first-dimension optimazation
HDMR =@(x)HDMR_function(f,f_0,x,1); % HDMR function
dimension_number = 1;
TEST_ACCURACY
x_0_old = x_0_new; f_0_old = f_0_new; % update old x_0 and f_0
[x_0_dimension, f_0_predict] = ga(HDMR,1,[],[],[],[],lower_bound(K(1)),upper_bound(K(1)), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
x_0_new(selected_dimension) = x_0_dimension; % update x_0_new
f_0_new = objective_function(x_0_new); % f_0_new
for i = 1:1
    [f_0_new, ia] = min([sample_value{1,i} + repmat(f_0, length(sample_value{1,i}), 1); f_0_new]);
    if(ia <= length(sample_value{1,i}))
        x_0_new = sample_point{1,i}(ia,:);
        x_0_new = point1DtoND(x_0_new, x_0, K(i));
    end
end
if (f_0_new > f_0_old)
    f_0_new = f_0_old;
    x_0_new = x_0_old;
end
[ f{1,1}, sample_point{1,1}, sample_value{1,1}] = ReconstructComponentFunciton( sample_point{1,1}, sample_value{1,1}, x_0_dimension, f_0_new - f_0_old);
NOE = NOE + 1;
% what if f_0_new is bigger than f_0
newbound = shrinkspace(f_0_old,f_0_new,f_0_predict,x_0_new,x_0_old,...
    lower_bound,upper_bound,lower_bound_current,upper_bound_current);
lower_bound_current(K(1)) = newbound(1,K(1));
upper_bound_current(K(1)) = newbound(2,K(1));
f_0_new

%% construct remained member functions
dimension_number = 1;
second_order_count = 0;
choose = zeros(nchoosek(dimension,2),2); % initialize choose array, each row represents the current two dimension need to be chosen
while dimension_number < dimension && dimension_number < 16 && NOE < 2900 %detect(f_0_new,f_0_old) == 1
    %% 
    dimension_number = dimension_number + 1; % consider which number of dimension
    
    %% roulette wheel selection
    selected_dimension = roullete(S0);
    K(dimension_number) = selected_dimension
    S0(K(dimension_number)) = 0;
    
    %% first-order member function construction
    i = 1;j = dimension_number; % first order && first dimension
    x_0 = x_0_new; f_0 = f_0_new; % renew the zeroth order point and value
    [NOE,f{i,j},sample_point,sample_value] = construct_first_order(j,K(dimension_number),sample_point,sample_value,x_0,f_0,lower_bound_current,upper_bound_current,objective_function,NOE); % construct first-order member functions of first dimension
    sample_point_whole{i, j} = point1DtoND(sample_point{i,j},x_0,K(dimension_number));
    if dimension_number == 2
        HDMR =@(x)HDMR_function(f,f_0,x,1);
    else
        HDMR = @(x)HDMR_function_2(f,f_0,x,dimension_number - 1,second_order_count,K(1:end - 1));
    end
    
    %% dimension_concerned-order member function construction
    i = 2; % second-order function to be biuld
    second_order_count = second_order_count + 1;
     [NOE,f{i,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}] = construct_second_order(HDMR,f{1,dimension_number}, sample_point_whole,sample_value,x_0,f_0,objective_function,NOE,K,lower_bound_current, upper_bound_current, second_order_count, sample_point, x_0_old);
%     count = 0;
%     err_count = 0;
%     while count == err_count
%         try
%             [NOE,f{i,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}] = construct_second_order_substitute(HDMR,f{1,dimension_number}, sample_point_whole,sample_value,x_0,f_0,objective_function,NOE,K,lower_bound_current, upper_bound_current, second_order_count, sample_point, x_0_old);
%         catch MyErr
%             err_count = err_count + 1;
%         end
%         count = count + 1;
%     end
    sample_point_whole{i, second_order_count} = point1DtoND(sample_point{i,second_order_count}, x_0, K(1:dimension_number));
    
    %% function constructed in first-dimension optimazation
    HDMR = @(x)HDMR_function_2(f,f_0,x,dimension_number,second_order_count,K); % HDMR function
    %% test fidelity
    x_0_old = x_0_new; f_0_old = f_0_new; % update old x_0 and f_0
    [x_0_dimension, f_0_predict] = ga(HDMR,dimension_number,[],[],[],[],lower_bound(K(1:dimension_number)),upper_bound(K(1:dimension_number)), [], ga_opts);  % optimize HDMR functions and get new x_0 and f_0
    f_0_predict
    for i = 1:dimension_number
        x_0_new(K(i)) = x_0_dimension(i); % update x_0_new
    end
    f_0_new = objective_function(x_0_new); % f_0_new
    %% biuld the sample and value library
    value = objective_function(sample_point_whole{1,1});
    sample = sample_point_whole{1,1};
    for i = 2:dimension_number
        value = [value; objective_function(sample_point_whole{1,i})];
        sample = [sample; sample_point_whole{1,i}];
    end
    for i = 1:second_order_count
        value = [value; objective_function(point1DtoND(sample_point{2,i}, x_0, K(1:i+1)))];
        sample = [sample; point1DtoND(sample_point{2,i}, x_0, K(1:i+1))];
    end
  %% find f_0_new  
    [f_0_new, ia] = min([value; f_0_new]);
    if(ia <= length(value))
        x_0_new = sample(ia,:);
    end
    new_value = f_0_new - HDMR(x_0_dimension) + f{2,second_order_count}.func(x_0_dimension);
    if (f_0_new > f_0_old)
        f_0_new = f_0_old;
        x_0_new = x_0_old;
        new_value = 0;
    end
    TEST_ACCURACY
    newbound = shrinkspace(f_0_old,f_0_new,f_0_predict,x_0_new,x_0_old,...
    lower_bound,upper_bound,lower_bound_current,upper_bound_current);
    lower_bound_current(K(1:dimension_number)) = newbound(1,K(1:dimension_number));
    upper_bound_current(K(1:dimension_number)) = newbound(2,K(1:dimension_number));
%     TEST_ACCURACY
%     [NOE, f{2,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}, HDMR] = Addpoint(NOE, x_0_new, HDMR, objective_function, dimension_number, K, lower_bound_current, upper_bound_current, f{2,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}, x_0_old);% add points along the new cut line
%     sample_point_whole{2, second_order_count} = point1DtoND(sample_point{2,second_order_count}, x_0_old, K(1:dimension_number));
%     HDMR = @(x)HDMR_function_2(f,f_0,x,dimension_number,second_order_count,K);
%     count = 0;
%     err_count = 0;
%     while count == err_count
%         try
%              [NOE, f{2,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}, HDMR] = Addpoint(NOE, x_0_new, HDMR, objective_function, dimension_number, K, lower_bound_current, upper_bound_current, f{2,second_order_count}, sample_point{2, second_order_count}, sample_value{2, second_order_count}, x_0_old);% add points along the new cut line
%              sample_point_whole{2, second_order_count} = point1DtoND(sample_point{2,second_order_count}, x_0_old, K(1:dimension_number));
%         catch MyErr
%             err_count = err_count + 1;
%         end
%         count = count + 1;
%     end
    [ f{2,second_order_count}, sample_point, sample_value] = ReconstructComponentFunciton3( sample_point, sample_value, x_0_new(K(1:second_order_count + 1)), new_value, f_0_old, f_0_new, second_order_count);
    NOE = NOE + 1;
%     if detect(f_0_old,f_0_new) == 1
%         pause; 
%     end
    f_0_new
    NOE
end
f_0_new
NOE


%% evaluate HDMR functions' performance
% which is not useful now
% S = lhsdesign(100,dimension);
% S = S.*repmat(upper_bound - lower_bound,100,1) + repmat(lower_bound,100,1);
% for i = 1:100
%     point_temp = S(i,:);
%     point_HDMR = point_temp(K(1:min(5:end)));
%     obj_value(i) =  objective_function(point_temp);
%     HDMR_value(i) = HDMR(point_HDMR); 
% end
% mean = sum(obj_value)/100;
% temp1 = 0;temp2 = 0;
% for i = 1:100
%     temp1 = temp1 + (obj_value(i) - HDMR_value(i)).^2;
% end
% for i = 1:100
%     temp2 = temp2 + (obj_value(i) - mean).^2;
% end
% R_squre = 1 - temp1/temp2;
% 
% 
% STD = std(obj_value);
% temp1 = 0;
% for i = 1:100
%     temp1 = temp1 + norm(obj_value(i) - HDMR_value(i));
% end
% RAAE = temp1/100/STD;
% 
% RMAE = max(abs(obj_value - HDMR_value))/STD;
% NOE;
% R_squre
% RAAE
% RMAE
% 
% dimension = 2;
% S = lhsdesign(100,dimension);
% S = S.*repmat(upper_bound_current(1:dimension) - lower_bound_current(1:dimension),100,1) + repmat(lower_bound_current(1:dimension),100,1);
% for i = 1:100
%     point_temp = S(i,:);
%     point_HDMR = point_temp;
%     obj_value(i) =  objective_function(point1DtoND(point_temp, x_0_new, K(1:dimension)));
%     HDMR_value(i) = HDMR(point_HDMR); 
% end
% mean = sum(obj_value)/100;
% temp1 = 0;temp2 = 0;
% for i = 1:100
%     temp1 = temp1 + (obj_value(i) - HDMR_value(i)).^2;
% end
% for i = 1:100
%     temp2 = temp2 + (obj_value(i) - mean).^2;
% end
% R_squre = 1 - temp1/temp2;
% 
% 
% STD = std(obj_value);
% temp1 = 0;
% for i = 1:100
%     temp1 = temp1 + norm(obj_value(i) - HDMR_value(i));
% end
% RAAE = temp1/100/STD;
% 
% RMAE = max(abs(obj_value - HDMR_value))/STD;
% NOE;
% R_squre
% RAAE
% RMAE

