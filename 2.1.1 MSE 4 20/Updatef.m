function [f_new] = Updatef(sample_point_whole, sample_value, dimension_number, x_0, K, f)
%% load sample point
point = sample_point_whole{1,1};
valu = sample_value{1,1};
count = 1;
for i = 2:dimension_number
    if size(sample_value{1,i},1) > 0
        point = [point;sample_point_whole{1,i}];
        valu = [valu;sample_value{1,i}];
    end
    for j = 1:i - 2
        if size(sample_value{2,count},1) > 0
            point = [point;sample_point_whole{2,count}];
            valu = [valu;sample_value{2,count}];
        end
        count = count + 1;
    end
end
[m, n] = size(point);
theta = repmat(10, n, 1); lob = repmat(1e-1, n, 1); upb = repmat(20, n, 1);
[dmodel, perf] = dacefit(point, valu, @regpoly0, @corrgauss, theta, lob, upb);
f{1,1}.dmodel = dmodel;
f{1,1}.perf = perf;
f{1,1}.mse = @(x)-Get_MSE(@predictor, pre.dmodel, x);
f{1,1}.func = @(x)predictor(x,pre.dmodel);
for i = 1:dimension_number
    f{1,i}.func = @(x)predictor(point1DtoND(x,x_0,K(i)),f{1,1}.dmodel);
end
f_new = f;
end





