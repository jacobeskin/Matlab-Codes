function [archB, garchB] = GArkkiOLS(X,a,b,c)
% Takes in datavector X, does ARCH and GARCH. a is X's own lag,
% b is lag of the errors, c is lag of the variances. Everything is
% done with OLS. Outputs are the coefficients for ARCH (udiff = K1*u(t-1).^2+...Kb*u(t-b).^2)
% in archB-vector, coefficients of GARCH (udiff =
% K1*u(t-1).^2+...Kb*u(t-b).^2+C1*udiff(t-1)+...+Cc*udiff(t-c)) 
% in garchB-vector. Plots results. 
% WARNING for now b>=c!! NEEDS FURTHER DEVELOPMENT 
tic

% Muodostetaan matriisi M vektorista X ja sen lageista a
[mat,~] = size(X);
M = [X(a+1:end,1) zeros(mat-a,1)];
for i=1:a
    Xi = matlag(X,i); % lagataan X:aa
    Xi = Xi(a-i+1:end,1); % tehdaan oikean mittainen
    M(:,i+1) = Xi; % tallennetaan 
end

% Lasketaan OLS y = x*B+u ja otetaan talteen errorit ja lasketaan errorin
% muutoksen suuruus udiff =(u(t)-u(t-1))^2
y = M(:,1);
x = M(:,2:end);
B = (x'*x)\(x'*y);
u = y-x*B;
udif = u(2:end)-u(1:end-1);
udiff = udif.^2; % errorin muutoksen suuruus^2 laskettu
[n,~] = size(udiff); % udiffin pituus
u = u(2:end,:); % saman mittaiseksi kuin udiff
u = u.^2; % error potenssiin 2
AR = x*B;

% ARCH, udiff = K1*u(t-1).^2+...Kb*u(t-b).^2
D = zeros(n-b, b); % luodaan matriisi errorien lageille
for j=1:b % sijoitetaan errorin lagit matriisiin D 
    ui = matlag(u,j);
    ui = ui(b-j+1:end,:);
    D(:,j) = ui;
end
udiffARCH = udiff(b+1:end,1); % lyhennet‰‰n udiffi‰
K = (D'*D)\(D'*udiffARCH); % lasketaan beettakertoimet udiff = D(lagit)*K
archB = K;

ARCH = D*K; % Lasketaan mallin m‰‰r‰‰m‰ arvo udifille plottaamista varten

% GARCH (udiff =
% K1*u(t-1).^2+...Kb*u(t-b).^2+C1*udiff(t-1)+...+Cc*udiff(t-c))    ¥
A = zeros(n-max(b,c), b+c); % luodaan matriisi errorin ja udiffin lageille riippuen siit‰ kumman lageja on enemm‰n
if b<c % varianssin l‰gej‰ enemm‰n
    for j=1:b % A matriisiin errororien l‰geist‰
        ui = matlag(u,j);
        ui = ui(b+c-j:end,:);
        A(:,j) = ui;
    end
    for j=1:c % A matriisiin udiffin l‰git
        udi = matlag(udiff,j);
        udi = udi(c-j+1:end,:);
        A(:,b+j) = udi;
    end
end

if c<b % errorin l‰gej‰ enemm‰n
    for j=1:b % A matriisiin errororien l‰geist‰
        ui = matlag(u,j);
        ui = ui(b-j+1:end,:);
        A(:,j) = ui;
    end
    for j=1:c % A matriisiin udiffin l‰git
        udi = matlag(udiff,j);
        udi = udi(c+b-j:end,:);
        A(:,b+j) = udi;
    end
end

if b==c % l‰gej‰ saman verran
   for j=1:b % A matriisiin errororien l‰geist‰
       ui = matlag(u,j);
       ui = ui(b-j+1:end,:);
       A(:,j) = ui;
   end
   for j=1:c % A matriisiin udiffin l‰git
       udi = matlag(udiff,j);
       udi = udi(c-j+1:end,:);
       A(:,b+j) = udi;
   end
end
udiffGARCH = udiff(max(b,c)+1:end,1);
C = (A'*A)\(A'*udiffGARCH); % lasketaan beettakertoimet udiffGARCH = A*C
garchB = C;

GARCH = A*C; % Lasketaan mallin m‰‰r‰‰m‰ arvo udifille plottaamista varten 

% Tehd‰‰n kaikesta saman mittaista plottausta varten
[AA,~] = size(AR);
[ar,~] = size(ARCH);
[gar,~] = size(GARCH);
[uu,~] = size(udiff);

m = 10;

X = X(m:end,:);
AR = AR(m-mat+AA:end,1);
ARCH = ARCH(m-mat+ar:end,1);
GARCH = GARCH(m-mat+gar:end,1);
udiff = udiff(m-mat+uu:end,1);
z = 1:length(udiff);

% plotataan tulokset
subplot(3,1,1)
plot(z,X,'b',z,AR,'r')
grid MINOR
title('Hinnan kehitys')
subplot(3,1,2)
plot(z,udiff,'b',z,ARCH,'r')
grid MINOR
title('volatiliteetti ARCHilla')
subplot(3,1,3)
plot(z,udiff,'b',z,GARCH,'r')
grid MINOR
title('volatiliteetti GARCHilla')
toc
end
