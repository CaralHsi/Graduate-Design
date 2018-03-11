function [ output_args ] = loadfunc( input_args,x)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
[m n] = size(x);
c = [-6.089,-17.164,-34.054,-5.914,-24.721,-14.986,-24.100,-10.708,-26.662,-22.179];
switch input_args
    case 1
        output_args = x(:,2).^2 + x(:,1).*x(:,3) + x(:,1) - repmat(4,m,1);
    case 2
        % dimension 10 lb 2.1 ub 9.9
        SUM = 0;
        for i = 1:n
            SUM = SUM + (log(x(:,i) - repmat(2,m,1))).^2 + (log(repmat(10,m,1) - x(:,i))).^2;
        end
        SUM = SUM - prod(x')'.^(0.2);
        output_args = SUM;
    case 3
        SUM = 0;
        for i = 1:n
            SUM = SUM + x(:,i).*(repmat(c(i),m,1) + log((x(:,i))./sum(x,2)));
        end
        output_args = SUM;
    case 4
        SUM = 0;
        for i = 1:n
            temp = 0;
            for j = 1:n
                temp = temp + exp(x(:,j));
            end                
            SUM = SUM + exp(x(:,i)).*(repmat(c(i),m,1) + x(:,i) - log(temp));
        end
        output_args = SUM;
    case 5
        output_args = x(:,1).^2 + x(:,2).^2 + x(:,1).*x(:,2) - 14*x(:,1)-...
            16*x(:,2) + (x(:,3) - repmat(10,m,1)).^2+ 4*(x(:,5) - repmat(5,m,1)).^2 +... 
            (x(:,5) - repmat(3,m,1)).^2 + 2*(x(:,6) - repmat(1,m,1)).^2 + ...
            5*(x(:,7)).^2 + 7*(x(:,8) - repmat(11,m,1)).^2 +...
            2*(x(:,9) - repmat(10,m,1)).^2 + (x(:,10) - repmat(7,m,1)).^2 + ...
            +45;
    case 6
        % F16 16 lb ub [-1 1]
        SUM = 0;
        a  = ...
        [1 0 0 1 0 0 1 1 0 0 0 0 0 0 0 1;
        0 1 1 0 0 0 1 0 0 1 0 0 0 0 0 0;
        0 0 1 0 0 0 1 0 1 1 0 0 0 1 0 0;
        0 0 0 1 0 0 1 0 0 0 1 0 0 0 1 0;
        0 0 0 0 1 1 0 0 0 1 0 1 0 0 0 1;
        0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0;
        0 0 0 0 0 0 0 1 0 1 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1;
        0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
        for i = 1:n
            for j = 1:n
                SUM = SUM + a(i,j) * (x(:,i).^2 + x(:,i) + 1).*(x(:,j).^2 + ...
                x(:,j) + 1);
            end
        end
        output_args = SUM;
    case 7
        % 3 dimensions lb 0 ub 1
        alpha = [1 1.2 3 3.2]';
        A = [3 10 30; 0.1 10 35; 3 10 30; 0.1 10 35];
        P = 10^(-4)*[3689 1170 2673; 4699 4387 7470; 1091 8732 5547; 381 ...
            5743 8808];
        SUM = zeros(m,1);
        K = [3 1 2];
        x = x(:,K);
        for k = 1:m
            for i = 1:4
                SUM(k) = SUM(k) - alpha(i) * exp(-(A(i,:)*(x(k,:) - P(i,:))'.^2)');
            end
        end
        output_args = SUM;
    case 8
        % SUR-T1-14 function, n = 10 20 30
        n = 10;
        sum_temp = 0;
        for i = 1:n-1
            sum_temp = sum_temp + (n - i) * (x(:,i).^2 - x(:,i + 1)).^2;
        end
        SUM = (x(:,1) - ones(m,1)).^2 + (x(:,n) - ones(m,1)).^2 + ...
            n * sum_temp;
        output_args = SUM;
    case 9
        % Rosenbrock dimension 10 lb -5 ub 5
        SUM = zeros(m,1);
        for i = 1:9
            SUM(:,i) = 100 * (x(:,i+1) - x(:,i)).^2 + (x(:,i)-1).^2;
        end
        output_args = sum(SUM,2);
        
    case 10
        % Trid dimension 10 lb ub [-100,100]
        n = 10;
        SUM1 = zeros(m,1);
        SUM2 = zeros(m,1);
        for j = 1:n;
            SUM1 = SUM1+(x(:,j)-ones(m,1)).^2;    
        end
        for j = 2:n;
            SUM2 = SUM2+x(:,j).*x(:,j-1);    
        end
        output_args = SUM1-SUM2;
    case 11
        % Girewank [-600 600]
        n = 10;
        fr = 4000;
        s = zeros(m,1);
        p = ones(m,1);
        for j = 1:n
            s = s + x(:,j).^2;
        end
        for j = 1:n
            p = p .* cos(x(:,j)./sqrt(j));
        end
        output_args = s./fr - p + ones(m,1);
    case 12
        %Ackley  10 (20 30)[-30,30]
        n = 20;
        a = repmat(20,m,1);
        b = 0.2; c = 2*pi;
        s1 = zeros(m,1); s2 = zeros(m,1);
        for i = 1:n;
            s1 = s1 + x(:,i).^2;
            s2 = s2 + cos(c.*x(:,i));
        end
        output_args = -a .* exp(-b .* sqrt(1 / n .* s1)) - exp(1 / n .* s2) + a + repmat(exp(1),m,1);
    case 13
        % Rastrigin dimension  [-5.12,5.12]
        n = repmat(10*20, m, 1); 
        s = zeros(m, 1);
        for j = 1:20
            s = s + (x(:,j).^2 - 10 .* cos(2 .* pi .* x(:,j))); 
        end
        output_args  = n + s;
    case 14
         % SUR-T1-16 dimension 20 lb ub [-2 5]
        f = zeros(m, 1);
        for i = 1:5
            f = f + (x(:,i) + 10 .* x(:,i+5)).^2 + 5 .* (x(:,i+10) - x(:,i+15)).^2 + ...
                (x(:,i+5) - 2 .* x(:,i+10)).^4 + 10 .* (x(:,i) - x(:,i+15)).^4;
        end
        output_args = f;
    case 15
        % Powell dimension 20 lb ub[-4,5]
        n = 20;
        SUM = zeros(m, 1);
        for i = 1:n/4
            t1 = x(:,4*i-3) + 10 .* (x(:,4*i-2));
            t2 = sqrt(5) .* (x(:,4*i-1) - x(:,4*i));
            t3 = (x(:,4*i-2) - 2 .* (x(:,4*i-1))).^2;
            t4  = sqrt(10) .* (x(:,4*i-3) - x(:,4*i)).^2;
            SUM = SUM + t1.^2 + t2.^2 + t3.^2 + t4.^2;
        end;
        output_args = SUM;
    case 16
        % perm dimension 20 [-20 20]
        n = 20;
        b = 0.5;
        SUM = zeros(m, 1);
        for k = 1:n;
            s_in = zeros(m, 1);
            for j = 1:n
                s_in = s_in + (j^k + b) .* ((x(:,j)./j).^k - ones(m, 1));
            end
            SUM = SUM + s_in.^2;
        end
        output_args = SUM;

%%
% f=x(1).*x(2)*100+x(3).*x(4)*100+x(5).*x(6)*100+x(7).*x(8)*100+x(9).*x(10)*100;

% 
% 
% end

% for i = 1:10
%     temp(i) = (log(x(i)-2))^2 + (log(10-x(i)))^2;
% end
% f = sum(temp) + (x(1) * x(2)* x(3)* x(4)* x(5)* x(6)* x(7)* x(8)* x(9)* x(10))^0.2;

%% PUR-T1-13
% for i = 1:10
%     temp(i) = i^3 * (x(i)-1)^2;
% end
% f = sum(temp)^3;

%%
% f = x(2)^2 + 2*x(1) * x(3) + x(1) + x(1) -x(3)- 4;
% f = x(1)^2 + x(1) * x(2) + 15 * x(2)^2+ x(1) +  x(2) * x(3) * x(1) + 5 * x(3)^2  - x(3) * x(1);
% f = x(1)^2 + x(1) * x(2) + x(2)^2+ x(1);
% f = 4 * x(1)^2 - 2.1 * x(1)^2 + 11/3 *x(1)^6 + x(1) * x(2) - 4 * x(2)^2 + 4 * x(2)^4;
% f = sin(x(1)+x(2)) + (x(1)-x(2))^2-1.5 * x(1)+2.5*x(2)+1;
% f = (x(1)-1)^2+(x(1) - x(2))^2 +(x(2)-x(3))^4;
% cd('ac range')
% Z = x(5:10);
% X1 = x(1:2);
% X2 = x(3);
% X3 = x(4);
% i0=[.25 1 1 .5 .05 45000 1.6 5.5 55 1000];
% [Y1,Y2,Y3,Y4,Y12,Y14,Y21,Y23,Y24,Y31,Y32,Y34,G1,G2,G3,C,Twist_initial,...
%     x_initial,L_initial,R_initial,ESF_initial,Cf_initial,Lift_initial,tc_initial,M_initial,h_initial,...
%     T_initial]=system_analyis(Z,X1,X2,X3,i0);
% f = -Y4;
% cd ..

% f=0;
% for i=1:9
%     f=f+(100*(x(i+1)-x(i).^2).^2+(x(i)-1).^2);
% end

%% 
% a(:,2)=10.0*ones(4,1);
% for j=1:2;
%    a(2*j-1,1)=3.0; a(2*j,1)=0.1; 
%    a(2*j-1,3)=30.0; a(2*j,3)=35.0; 
% end
% c(1)=1.0;c(2)=1.2;c(3)=3.0;c(4)=3.2;
% p(1,1)=0.36890;p(1,2)=0.11700;p(1,3)=0.26730;
% p(2,1)=0.46990;p(2,2)=0.43870;p(2,3)=0.74700;
% p(3,1)=0.10910;p(3,2)=0.87320;p(3,3)=0.55470;
% p(4,1)=0.03815;p(4,2)=0.57430;p(4,3)=0.88280;
% s = 0;
% for i=1:4;
%    sm=0;
%    for j=1:3;
%       sm=sm+a(i,j)*(x(j,:)-p(i,j)).^2;
%    end
%    s=s+c(i).*exp(-sm);
% end
% f = -s;

%% levy
% n = 4;
% for i = 1:n; z(i) = 1+(x(i,:)-1)/4; end
% s = sin(pi*z(1)).^2;
% for i = 1:n-1
%     s = s+(z(i)-1).^2*(1+10*(sin(pi*z(i)+1)).^2);
% end 
% f = s+(z(n)-1).^2*(1+(sin(2*pi*z(n))).^2);

%% Dixon
% n = 25;
% s1 = 0;
% for j = 2:n;
%     s1 = s1+j*(2*x(j)^2-x(j-1))^2;    
% end
% f = s1+(x(1,:)-1).^2;

%% airfoil
% cd airfoil
% N = 90;
% alfa = 2;
% dz = 0;
% f = Xfoilfun2(x',N,alfa,dz);
% cd ..

%% Zakharov [-5,10]
% n = 20;
% for i = 1:n
%   temp1(i) = x(i)^2;
%   temp2(i) = 0.5 * i * x(i);
%   temp3(i) = 0.5 * i * x(i);
% end
% f = sum(temp1) + sum(temp2)^2 + sum(temp3)^4;
            
            
end


end

