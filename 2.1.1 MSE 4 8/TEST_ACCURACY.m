dimension0 = dimension_number;
S = lhsdesign(500,dimension0);
S = S.*repmat(upper_bound_current(K(1:dimension0)) - lower_bound_current(K(1:dimension0)),500,1) + repmat(lower_bound_current(K(1:dimension0)),500,1);
for i = 1:500
    point_temp = S(i,:);
    point_HDMR = point_temp;
    obj_value(i) =  objective_function(point1DtoND(point_temp, x_0_new, K(1:dimension0)));
    HDMR_value(i) = HDMR(point_HDMR); 
end
mean = sum(obj_value)/500;
temp1 = 0;temp2 = 0;
for i = 1:500
    temp1 = temp1 + (obj_value(i) - HDMR_value(i)).^2;
end
for i = 1:500
    temp2 = temp2 + (obj_value(i) - mean).^2;
end
R_squre = 1 - temp1/temp2;


STD = std(obj_value);
temp1 = 0;
for i = 1:500
    temp1 = temp1 + norm(obj_value(i) - HDMR_value(i));
end
RAAE = temp1/100/STD;

RMAE = max(abs(obj_value - HDMR_value))/STD;
NOE;
R_squre
RAAE
RMAE

