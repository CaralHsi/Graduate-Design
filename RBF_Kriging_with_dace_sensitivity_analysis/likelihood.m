function [NegLnLike, R, U] = likelihood( theta )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
global ModelInfo
S = ModelInfo.S;
Y = ModelInfo.Y;
[m,n] = size(S);
theta = theta; % theta = 10.^theta;don't know why!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Calculate distances D between points
mzmax = m*(m-1) / 2;        % number of non-zero distances
ij = zeros(mzmax, 2);       % initialize matrix with indices
D = zeros(mzmax, n);        % initialize m[dmodel, perf] = ... dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb)atrix with distances
ll = 0;
for k = 1 : m-1
  ll = ll(end) + (1 : m-k);
  ij(ll,:) = [repmat(k, m-k, 1) (k+1 : m)']; % indices for sparse matrix
  D(ll,:) = repmat(S(k,:), m-k, 1) - S(k+1:m,:); % differences between points
end
if  min(sum(abs(D),2) ) == 0
    error('Multiple design sites are not allowed')
end
D_square = D.^2;
% calculate correlation matrix R
R = constructR(theta,D_square,m);
one = ones(m,1);
%Choleskey factorization
[U, p] = chol(R);
if p > 0
    NegLnLike = 1e4;
else
    %Sum Lns of diagonal to find ln(det(R))
    LnDetR = 2 * sum(log(abs(diag(U))));
    
    %Use back-substitution of Cholesky instead of inverse
    mu = (one' * (U\(U'\Y)))/(one'*(U\(U'\one)));
    SigmaSqr = ((Y - one * mu)'*(U\(U'\(Y - one * mu))))/n;
    NegLnLike = -1 * (-(n/2) * log(SigmaSqr) - 0.5 * LnDetR); 

end

