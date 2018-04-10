function [NOE, pre, sample_point_temp, sample_value_temp, HDMR] = Addpoint(NOE, x_0_new, HDMR, objective_function, dimension_number, K, lower_bound, upper_bound, pre0, sample_point_temp, sample_value_temp, x_0_old)% add points along the new cut line
  
%Addpoint function adds some points along the new cut line to make the
%prediction of the value near the new cut center more accurate.
%   此处显示详细说明
for i = 1:dimension_number
    t = rand();
    test_point = t.*(upper_bound(K(i)) - lower_bound(K(i))) + lower_bound(K(i));
    test_point = point1DtoND(test_point, x_0_new, K(i));
    test_value = objective_function(test_point) - HDMR(test_point(K(1:dimension_number)));
    NOE = NOE + 1;
    while (abs(test_value/(objective_function(test_point) + eps)) > 0.01)
        sample_point_temp = [sample_point_temp; test_point(K(1:dimension_number))];
        sample_value_temp = [sample_value_temp; test_value + pre0.func(test_point(K(1:dimension_number)))];
        [sample_point_temp,ia,~] = unique(sample_point_temp,'rows');
        sample_value_temp = sample_value_temp(ia);
        [dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, ones(1, dimension_number),repmat(1e-3, 1, dimension_number), repmat(100, 1, dimension_number));
        pre.dmodel = dmodel;
        pre.perf = perf;
        pre.func = @(x)predictor(x, pre.dmodel);
        HDMR = @(x)HDMR(x) - pre0.func(x) + pre.func(x);
        pre0 = pre;
        t = rand();
        test_point = t.*(upper_bound(K(i)) - lower_bound(K(i))) + lower_bound(K(i));
        test_point = point1DtoND(test_point, x_0_new, K(i));
        test_value = objective_function(test_point) - HDMR(test_point(K(1:dimension_number)));  
        NOE = NOE + 1;
    end
    sample_point_temp = [sample_point_temp; test_point(K(1:dimension_number))];
    sample_value_temp = [sample_value_temp; test_value + pre0.func(test_point(K(1:dimension_number)))];
    [sample_point_temp,ia,~] = unique(sample_point_temp,'rows');
    sample_value_temp = sample_value_temp(ia);
    [dmodel, perf] = dacefit(sample_point_temp, sample_value_temp, @regpoly0, @corrgauss, ones(1, dimension_number),repmat(1e-3, 1, dimension_number), repmat(100, 1, dimension_number));
    pre.dmodel = dmodel;
    pre.perf = perf;
    pre.func = @(x)predictor(x, pre.dmodel);
    HDMR = @(x)(HDMR(x) - pre0.func(x) + pre.func(x));
        
end

end

