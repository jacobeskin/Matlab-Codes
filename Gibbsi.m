function [Bmean, Bstd, Smean, Sstd, Rmean, Varianssi] = Gibbsi(X, simulaatioita, burnin)
% Heteroscedastic Gibbs sampling for b-terms in y=xb+e. Input is datamatrix X, where 1st 
% column is vector y, rest are the matrix x, number of Monte Carlo runs and burn-in for 
% MCMC. Output is mean and standard deviation of b:s, mean and standard deviation of 
% the variance, mean and variance of the residual. 

tic
% Break up matrix X and calculate the starting values from regular OLS
[n,m] = size(X);
if mod(n,2)~=0 % make X have even number of rows
    X = X(2:end,:);
    [n,m] = size(X);
end
y = X(:,1); 
x = X(:,2:m); 
b0 = (x'*x)\(x'*y);

% Calculate a distribution for b:s and their expectation value vector and 
% covariance matrix
beetat = zeros(n/2,m-1);
for i=1:(n/2) 
    yi = y(i:(i+n/2),:);
    xi = x(i:(i+n/2),:);
    bi = (xi'*xi)\(xi'*yi);
    beetat(i,:) = bi';
end
r = mean(beetat); % distribution of b:s should be something like ~N(r,T) 
T = cov(beetat);
sigma = ((y-x*b0)'*(y-x*b0))/(n-m-1); % calculate variance
Q = chol(inv(T)); % inv(T)=Q'*Q
q = Q*r';
beettamatriisi = zeros(simulaatioita, m-1); % for calculation of b:s mean (MCMC)
sigmavektori = zeros(simulaatioita,1); % save variances here
residuaali = y-x*b0;

    %Heteroscedastic sampling
    rval = 4; % starting value for rval which determines new elements of V
    E = (residuaali./sigma).*residuaali+rval;
    vi = zeros(n,1);
    for k=1:n 
        vi(k,1) = (1/chi2rnd(rval+1,1,1)).*E(k,1); % starting values vi
    end
    V = diag(vi);
    Var = zeros(simulaatioita,n);
    error = zeros(simulaatioita,n);
    for j=1:simulaatioita
        xs = x'*(V\x);
        xys = x'*(V\y);
        SIG = sigma*(Q'*Q);
        b = ((xs+SIG)\(xys+sigma*(Q'*q))); %exp. val. for b w/ condition sigma and V
        bj = mvnrnd(b, sigma./(xs+SIG)); % generate new b
        beettamatriisi(j,:) = bj; % save it
        e = y-x*bj'; % calculate error
        error(j,:) = e';
        sigma = ((e.*e)'*(1./vi))/chi2rnd(n); % update variance
        sigmavektori(j,1) = sigma;
        E = (e./sigma).*e+rval;
        for k=1:n 
            vi(k,1) = (1/chi2rnd(rval+1,1,1)).*E(k,1); % update vi
        end  
        V = diag(vi); % update V
        Sigma = ones(n,1)*sigma;
        Var(j,:) = (V*Sigma)';
        rval = gamrnd(8,2,1,1); % update rval 
    end

% Calculate mean and standard deviation of b and variance and other output
Beetta = beettamatriisi((burnin+1):end,:);
Bmean = mean(Beetta);
Bstd = std(Beetta);
S = sigmavektori((burnin+1):end,1);
Smean = mean(S);
Sstd = std(S);
Var = Var((burnin+1):end,:);
Varianssi = mean(Var);
ERR = error((burnin+1):end,:);
Rmean = mean(ERR);
toc
end
        
        
        
    
    


    
    
    
    
    





