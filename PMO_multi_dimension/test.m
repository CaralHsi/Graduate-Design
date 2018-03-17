for i = -5:0.5:5
    for j = -5:0.5:5
        Z(:,)
 
hold on;
[X,Y] = meshgrid(-5:.1:5);
[m, n] = size([X]);
Z = (X.^2 + Y - repmat(11,m,n)).^2 + (X + Y.^2 - repmat(7,m,n)).^2;
mesh(X, Y, Z);

hold on;
plot3(x_0_new(1),x_0_new(2),f_0_new,'*');

hold on;
c = 1;
for i = -5:.1:5
    y(c) = f{1,2}.func(i) + f_0_new;
    c = c + 1;
end
plot3(-5:.1:5, repmat(x_0_new(2),1,101), y);

figure;

[X,Y] = meshgrid(-5:.1:5);
[m, n] = size([X]);
for i = 1:101
    for j = 1:101
        Z(i, j) = pre.func([X(i ,j) Y(i, j)]);
    end
end
mesh(X, Y, Z);
hold on;
plot3(test_point(1),test_point(2),pre.mse(test_point),'*');
