function [S,V] = sensitivityAnalyse(dim, HDMR, lower_bound, upper_bound,K)
%     function id = sensitivityAnalyse(structMatrix, rbfHdmrModel, dim)
% func is the function handle
% dim is the dimension of the problem
% s is the number of function phi
% lb, ub are the lower and upper boundary of the space respectively
% id is the index of the top 2 sensitive variables 

nSample = 3000;
x = lhsdesign(dim,nSample);
lb = -1 * ones(dim,1);
ub = ones(dim,1);
for i = 1:dim
    x(i,:) = lb(i) + (ub(i) - lb(i)) .* x(i,:);
end
tmp1 = zeros(dim,nSample);
tmp2 = zeros(dim,nSample);
y = zeros(1,nSample);
for j = 1:nSample
    for i = 1:dim
        if ~any(K == i)
            y(j) = HDMR(lower_bound + (upper_bound - lower_bound).*(x(:,j)' + ones(1,dim))/2);
            % y(j) = blackboxFuc(x(:,j));
            tmp1(i,j) = y(j) * Phi1(x(i,j));
            tmp2(i,j) = y(j) * Phi2(x(i,j));
        else
            tmp1(i,j) = 0;
            tmp2(i,j) = 0;      
        end
    end
end
c1 = zeros(1,dim);
c2 = zeros(1,dim);
V = zeros(1,dim);
for i = 1:dim
    c1(i) = sum(tmp1(i,:)) / nSample;
    c2(i) = sum(tmp2(i,:)) / nSample;
    V(i) = sqrt(c1(i)^2 + c2(i)^2);
end
S = V/sum(V);
% cp(1) = S(1);
% for i = 2:dim
%     cp(i) = cp(i-1) + S(i);
% end
% tmp = rand;
% for i = 1:dim
%     if tmp < cp(i)
%         id = i;
%         break;
%     end
% end
% 
% [v, I] = sort(S,'descend');
% id = I(1:2);
end



