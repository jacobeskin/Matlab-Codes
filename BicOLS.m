function [P, p] = BicOLS(V, l)
% Performs BIC (Bayesian Information Criterion) test on datamatrix V
% of size nx1 in order to determine the best number of lags for VAR 
% model. l is the number of maximum lags the test done for. Gives out
% value p for optimal amount of lags.

[n,m] = size(V);
if m>1
    error('Too many columns in V');
end
    
P = zeros(l,1); % vector for different amount of lags
X = zeros(n-l,l+1); % Matrix for holding the original and lagged data

% The formula is BIC(p)=log(e'*e/n)+(p+1)*log(n)/n.

% Making the correct length first column
W = V(l+1:end,1);
X(:,1) = W; 

for i=1:l
    Vi = V(1:end-i,1); % lagging V
    Vi = Vi(l-i+1:end,1); % cutting it into the right length
    X(:,i+1) = Vi; % save to X
    Y = X(:,1:i+1); 
    y = Y(:,1); 
    x = Y(:,2:end); 
    bi = (x'*x)\(x'*y);
    r = y-x*bi;
    P(i,1) = log((r'*r)/i) + (i+1)*log(i)/i; % save to BIC vector 
end

minP = min(P);

for j = 1:l
    if P(j:1)==minP 
        p = j;
    end
    break
end




    
    
    
    
