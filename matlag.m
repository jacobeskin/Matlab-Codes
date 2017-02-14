function Y = matlag(X, k)
% Lags matrix X k times

[n,m]=size(X);
Y = zeros(n-k,m);
for i=1:m
    Y(:,i) = X(1:end-k,i);
end
end

