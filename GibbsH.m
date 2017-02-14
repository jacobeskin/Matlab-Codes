function [Bmean, Bstd, R, SIG] = GibbsH(X, simulaatioita, burnin)
% Homoskedastic Gibbs sampling for b-terms in y = xb+e. Input is 
% datamatrix X, number of simulations, burn-in number for MCMC. 
% Gives out mean and standard deviation of b, mean of the residual 
% and mean of the variance.
tic
[n,m] = size(X);
if mod(n,2)~=0 % make X have even number of rows
    X = X(2:end,:);
    [n,m] = size(X);
end

y = X(:,1); % separate vector y from datamatrix X
x = X(:,2:m); % separate matrix x from datamatrix X
b0 = (x'*x)\(x'*y); % b from OLS
yols = x*b0; % vector y calculated with b0

% Lasketaan jokin jakauma beetoille ja siita niiden odotusarvovektori seka
% kovarianssimatriisi

% Calculate a distribution for b:s and their expectation valuevector and 
% covariance matrix

beetat = zeros(n/2,m-1); % vector for b:s
for i=1:(n/2) 
    yi = y(i:(i+n/2),:);
    xi = x(i:(i+n/2),:);
    bi = (xi'*xi)\(xi'*yi);
    beetat(i,:) = bi';
end
r = mean(beetat); %  distribution of b:s should be ~N(r,T) 
T = cov(beetat);
sigma = ((y-x*b0)'*(y-x*b0))/(n-m-1); % sigma-term
Q = chol(inv(T)); % inv(T)=Q'*Q
q = Q*r';
beettamatriisi = zeros(simulaatioita, m-1); % this will have the mean of b 
residuaalimatriisi = zeros(simulaatioita, n); % this will have the residuals
Varianssit = zeros(simulaatioita, 1);  % this will have the variance terms
% Homoscedastic sampling
for j=1:simulaatioita
    b = inv(x'*x+sigma*(Q'*Q))*(x'*y+sigma*Q'*q); % exp.val. for b w/ condition sigma
    bj = mvnrnd(b,sigma*inv(x'*x+sigma*(Q'*Q))); % generate new b from distribution
    beettamatriisi(j,:) = bj; % save b
    e = y-x*bj'; % calcualte error
    residuaalimatriisi(j,:) = e'; % save residuals
    sigma = (e'*e)/chi2rnd(n);% update sigma
    Varianssit(j,1) = sigma; % save sigma
end

Beetta = beettamatriisi((burnin+1):simulaatioita,:);
Bmean = mean(Beetta);
Bstd = std(Beetta);
ERR = residuaalimatriisi((burnin+1):simulaatioita,:);
R = mean(ERR);
Var = Varianssit((burnin+1):end,:);
SIG = mean(Var);

toc
end
