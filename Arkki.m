function [K, b, V, E, R] = Arkki(X)
% Performs ARCH for datavector X, gives out number of lags b for errors b
% and coefficients K

tic
[n,~] = size(X);

% Perform BIC for vector X for optimal lags
[~, p] = BicOLS(X, 5);

% Create matrix Y from data X and its lags
Y = [X(p+1:end,1) zeros(n-p,p)];
for i=1:p
    Xi = matlag(X,i); % lag X
    Xi = Xi(p-i+1:end,1); % cut to right length
    Y(:,i+1) = Xi; % save
end

% Calculate error vector and its variance vector 
[~, ~, ~, ~, Rmean, Varianssi] = Gibbsi(Y, 5500, 500);

% ARCH is sigma(t) = B1*error(t-1)^2+..+Bb*error(t-b)^2
% Use BIC to determine number of lags for errors  
Var = Varianssi';
R = Rmean'.*Rmean';
[~, b] = BicOLS([Var R], 1);

% Create datamatrix D from sigma and lagged values 
[n,~] = size(Var);
D = zeros(n-b, b+1);
D(:,1) = Var(b+1:end, :);

for j=1:b % lags of errors to matrix D 
    Ri = matlag(R,j);
    Ri = Ri(b-j+1:end,:);
    D(:,j+1) = Ri;
end

V = D(:,1);
R = R(b+1:end,1);

% Calculate coefficients for lagged values 
[Bmean, ~, ~, ~, ~] = Gibbsi(D, 5500, 500);

% Formula is now D(:,1)=D(:,2:end)*Bmean'. Denote this by K
K = Bmean';
E = D(:,2);
toc
end
