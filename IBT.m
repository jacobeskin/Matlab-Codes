function [ProfitLoss, CumSumV, Hmatrix, Pmatrix, NumberOfTrades, NN] = IBT(tiedosto, startDate, endDate, beetta, TN, kustannus, Mittaus)
% Backtest for simple strategy built on the idea of the Ising-model.
% Downloads price data, measures required parameters, calculates and plot P/L.

% "Tiedosto" is the excel file that has the ticker names for the instruments, 
% SPY being the first, "startDate" and "endDate" defien the time period, format
% is dd/mm/yyyy. "beetta" is a parameter of the Ising model, can be measured 
% or determined in some other way. "TN" is propability limit which has to be 
% crossed in order to open a position. "Kustannus" is transaction cost. 
% "Mittaus" is the number of datapoints for measuring correlation coefficients.

% "ProfitLoss", "CumSum", "Hmatrix", "Pmatrix" are matricies and vecors for
% investigating the results and error handling. "NumberOfTrades" is well, the
% number of trades. "NN" is number of trades made by the control strategy.

% Correlations are measured dynamically. "Mittaus" also effects the number of
% trading days in the data. Plots 2 lines, one is P/L of Ising model the other 
% P/L of the control model where positions are opened purely based on previous 
% days close.

% Functions and scripts used, in order, are:
%
% tic
% xlsread() 
% numel()
% getYahooDailyData() 3RD PARTY SCRIPT
% struct2cell()
% table2array()
% zeros()
% isequal
% error()
% intersect()
% sort()
% length()
% sign()
% corrcoef()
% sum()
% cumcum()
% figure
% plot()
% legend()
% toc


tic 


%---------- Fetch data and make it easier to handle ----------


[X, x] = xlsread(tiedosto);                                 
dim = numel(x);                                             % # of instruments
D = getYahooDailyData(x, startDate, endDate, 'dd/mm/yyyy'); % Get price data
D = struct2cell(D);                                        


%---------- For speed, create some matrices ready for use ----------


SDATE = table2array(D{1}(:,1));    % Spider DATE info
OPEN = zeros(length(SDATE), dim);  % OPEN price matrix
CLOSE = zeros(length(SDATE), dim); % CLOSE price matrix
Y = zeros(size(OPEN));             % This will be change in price matrix
YY = zeros(size(Y));               % This will be relative change in price matrix
Z = length(SDATE);                 % Number of datapoints


% Prices of instruments into matricies
for i=1:dim
    
    % Make usre dates are correct
    
    DATE = table2array(D{i}(:,1)); % Instrument i DATE info
    equal = isequal(SDATE,DATE);   
    % If there is a problem
    if equal==0 
        error('IBT:BadDate', ['Days do not match for instrument ''' x{i} '''. Fix data!']); % error message
    end
    
    OPEN(:,i) = table2array(D{i}(:,2));    % Instrument i unadjusted open
    CLOSE(:,i) = table2array(D{i}(:,5));    % Instrument i unadjusted close
    
end


%---------- Build datasets for parameter measuring and backtesting ----------


% Check if enough datapoints 
if Z<=Mittaus
    error('IBT:MoreDates', ['Not enough days for trading and backtesting, smaller "Mittaus" or more days!']); 
end

% Changes in prices and "direction" matrix
for i=1:dim
    Y(:,i) = CLOSE(:,i)-OPEN(:,i);   % Change in prices
    YY(:,i) = CLOSE(:,i)./OPEN(:,i); % Relative change in prices 
end
K = sign(Y);        % Matrix whose elements are + or -
    
 
% Correlation coefficeint matrix for testing
CORR = corrcoef(YY(1:Mittaus,:)); % Correlation coefficients from beginning of data

ZZ = Z-Mittaus+1;         % Length of datavector for testing
K1 = K(Mittaus+1:end,:);  % plus-minus data for backtesting
OPEN1 = OPEN(Mittaus+1:end,:);   % New Open- price matrix
CLOSE1 = CLOSE(Mittaus+1:end,:); % New Close- price matrix


%---------- Ising model strategy ----------


% For speed, define matricies that are to be used
PL = zeros(ZZ,dim-1);     % P/L
PLK = zeros(ZZ, dim-1);   % Control P/L
[s1, s2] = size(K1);      % New dimensions of pricedata
PP = zeros(size(PL));     % P's of traded positions for later review
HH = zeros(size(PL));     % H's of traded positions for later review
N = 0;                    % # of trades 
KM = 0;
KO = 0;

% ising model, j goes throug days, i instruments
for j=1:s1-1
    
    for i=2:s2
        
        H = CORR(i,2:end)*K1(j,2:end)'*K1(j,i)-1+K1(j,1)*CORR(1,i)*K1(j,i); % The Hamiltonian
        P = 1/(1+exp((-2)*beetta*H));                                       % "probability" for price to be in this direction
        
        % P/L for control strategy
        if K1(j,i)<0
            PLK(j,i-1) = OPEN1(j+1,i)-CLOSE1(j+1,i)-2*kustannus;
            KM = KM+1;
        else if K1(j,i)>0
            PLK(j,i-1) = CLOSE1(j+1,i)-OPEN1(j+1,i)-2*kustannus;
            KO = KO+1;
            end
        end
        
        % P/L of opened position
        if (P>=TN)
            PP(j,i-1) = P;
            HH(j,i-1) = H;
            N = N+1;
            if K1(j,i)<0       % If this day was -1
            PL(j,i-1) = OPEN1(j+1,i)-CLOSE1(j+1,i)-2*kustannus; % P/L for short
            else if K1(j,i)>0  % Jos tämä päivä ollut +1
            PL(j,i-1) = CLOSE1(j+1,i)-OPEN1(j+1,i)-2*kustannus; % P/L for long
                end
            end
        end
    end
    
    % Update correlation coefficinents
    CORR = corrcoef(YY(j+1:Mittaus+j,:));
end


% ---------- Results and plots ----------


% Vectors for total P/L for each day and cumulative sum of P/L for Ising and control strategies
PLTOT = sum(PL,2);     
PLKTOT = sum(PLK,2);
PLKUM = cumsum(PLTOT); 
PLKKUM = cumsum(PLKTOT);

% Plots

x = length(PLKUM); % X-axis length
z = 1:x;           

figure                                                                         
plot(z, PLKUM, '-ro', z, PLKKUM, '-bx')                                        
legend('Ising','Control','Location','northoutside','Orientation','horizontal') 


% Outputs
CumSumV = PLKUM;
ProfitLoss = PL;
Hmatrix = HH;
Pmatrix = PP;
NumberOfTrades = N;
NN = KO+KM;

toc 


end

