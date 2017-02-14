function [D, B] = GarkkiV(X, a, b)
% Performs GARCH for times series vector X with a number of lags 
% of the error term and b number of lags of the variance term.
% Gives out B for the equation sigma=[errors,variances]*B and 
% datamatrix D.
tic
[n,~] = size(X);

% Form matrix Y from vector X ja its lags a
Y = [X(2:end,1) zeros(n-1,1)];
for i=1:1
    Xi = matlag(X,i); % laging X with function matlag()
    Xi = Xi(1-i+1:end,1); % cut it to right length
    Y(:,i+1) = Xi; % save
end

% Calcualte error vector and its variance vector using GibbsiH() 
[~, ~, ~, ~, Rmean, Varianssi] = Gibbsi(Y, 5500, 500);

% Garch is sigma(t) = B1*error(t-1)+..+Bb*error(t-b)+Q1*sigma(t-1)+..+Qq*sigma(t-q)
Var = Varianssi';
R = Rmean'.*Rmean';

% Form datamatrix D from sigma and the lagged values
[n,~] = size(Var);
D = zeros(n-max(a,b), a+b+1);
D(:,1) = Var(max(a,b)+1:end, :);

for j=1:a % substitute first lags of errors to D 
    Ri = matlag(R,j);
    Ri = Ri(max(a,b)-j+1:end,:);
    D(:,j+1) = Ri;
end

for k=1:b % substitute then lags of variances to D
    Vi = matlag(Var,k);
    Vi = Vi(max(a,b)-k+1:end,:);
    D(:,1+a+k) = Vi;
end

% Calculate coefficients for lagged values
[Bmean, ~, ~, ~] = GibbsH(D, 10000, 1000);

% Equation is now D(:,1)=D(:,2:end)*Bmean'. Denote it with B
B = Bmean';
toc
end









    
    





