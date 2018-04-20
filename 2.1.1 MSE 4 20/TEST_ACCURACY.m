dimension0 = dimension_number;
totalnumberof_test_points = 5000;
S = lhsdesign(totalnumberof_test_points,dimension0);
S = S.*repmat(upper_bound_current(K(1:dimension0)) - lower_bound_current(K(1:dimension0)),totalnumberof_test_points,1) + repmat(lower_bound_current(K(1:dimension0)),totalnumberof_test_points,1);
for i = 1:totalnumberof_test_points
    point_temp = S(i,:);
    point_HDMR = point_temp;
    obj_value(i) =  objective_function(point1DtoND(point_temp, x_0_new, K(1:dimension0)));
    HDMR_value(i) = HDMR(point_HDMR); 
end
mean = sum(obj_value)/totalnumberof_test_points;
temp1 = 0;temp2 = 0;
for i = 1:totalnumberof_test_points
    temp1 = temp1 + (obj_value(i) - HDMR_value(i)).^2;
end
for i = 1:totalnumberof_test_points
    temp2 = temp2 + (obj_value(i) - mean).^2;
end
R_squre = 1 - temp1/temp2;


STD = std(obj_value);
temp1 = 0;
for i = 1:totalnumberof_test_points
    temp1 = temp1 + norm(obj_value(i) - HDMR_value(i));
end
RAAE = temp1/totalnumberof_test_points/STD;

RMAE = max(abs(obj_value - HDMR_value))/STD;
NOE;
R_squre
RAAE
RMAE

