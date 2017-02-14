function [varargout] = ETFS(startDate, endDate, HalfLife, STD)
% Retrieves price data for ETF's in ETFS.xlxs and performs Johansen 
% cointegration test for all combinations of 2, 3 and 4 instruments,
% in an attempt to make stationary portfolios. Gets price data from Yahoo 
% finance. Outputs are excel files of the instruments and their weights, and 
% and plot of the time series. Rejects null hypotesis at 95%. Inputs are 
% start date and end date for the price data, "half life" of the decay to 
% average price, standard deviation from the standard price.
% USES 3rd PARTY SCRIPTS FOR RETRIEVING DATA AND PERFORMING COINTEGRAION TEST

tic
[A, B] = xlsread('ETFS.xlsx'); % Instruments listed in column B 
dim = numel(B);           % dim is number of instruments

% Get the price data into D
D = getYahooDailyData(B, startDate, endDate, 'dd/mm/yyyy');

% D from struct to cell
D = struct2cell(D); 

% Longest time series we have, need for comparing data to make sure vectors are same length
[X, ~] = size(table2array(D{1})); 
dimmax = X;

ll = 0;
z = 1:dimmax;      % need this later for plotting
montakoparia = 0;  % keeps count of stationary pairs
montakokolme = 0;  % keeps count of stationary triplets 
montakonelja = 0;  % keeps count of stationary quadruplets

%---------- First pairs, then triplets, then quadruplets ----------%


for i=1:dim-1 
    
    % Easier to handel arrays then tabels
    di = table2array(D{i}); % price data of i-th instrument
    Ci = di(:,7);           % Ci is adjusted closing price- vector
    
    % If time series is too short, skip it
    if (length(Ci)~=dimmax)
        ll = ll+1;
        continue
    end
    
    for j=(i+1):dim
        
        dj = table2array(D{j});
        Cj = dj(:,7);
        
        if (length(Cj)~=dimmax)
            continue
        end
        
        % Johansen test
        R = johansen([Ci Cj], 0, 1);
        r = struct2cell(R); 
        
	% We are interested in r(2) that has the weigts, r(3) and r(4) that 
	% have the trace statistic and eigen satistic, and r(4) and r(5)
	% that have the critical values for trace and eigen statistics
        
        
        TS = r{3}(2);    % Trace statisitc
        ES = r{4}(2);    % Eigen statistic
        TCV = r{5}(2,3); % Trace statistic critical values 99%
        ECV = r{6}(2,3); % Eigen statistic critical values 99%
        ki = r{2}(1,1);  % Instrument i weight
        kj = r{2}(2,1);  % Instrument j weight

        
        % Condition for cointegration
        if (TS>TCV) && (ES>ECV) % Trace statistic and Eigen statistic > 99% critical values
            
	    % "Half life" of the time series from LAMBDA, defined as 
            % deltaY(t) = LAMBDA*Y(t-1)
            
            I = ki*Ci; % Time series i with weights
            J = kj*Cj; % Time series j with weights
            IJ = I+J;  % Time series made from the above two 
            LIJ = IJ(1:end-1);         % Lagged IJ
            DIJ = IJ(2:end)-LIJ;       % IJ minus its lag
            LAMBDA = regress(DIJ,LIJ); % Perform OLS
            HL = -log(2)/LAMBDA;       % "Half life"
            
            % Condition from desired half life and standard deviation 
            if (LAMBDA<0) && (HL<HalfLife) && (STD<=std(IJ))
            names = strcat(char(B(i)), 'vs', char(B(j))); % Tickers of instruments of the pair
            O = cell(2);   
            O(1,1) = B(i); 
            O(2,1) = B(j); 
            O{1,2} = ki;
            O{2,2} = kj;
            xlswrite(names, O); % Output excel file
            montakoparia = montakoparia+1;
            
            % Plot and save time series
            f = figure('visible', 'off');
            plot(z,IJ, z, mean(IJ), 'r', z, mean(IJ)+std(IJ), 'b', z, mean(IJ)-std(IJ), 'b');
            title(names);
            saveas(f, names, 'jpg');
            
            end
        end
        % Same for triplets
        for k=j+1:dim
            
            dk = table2array(D{k});
            Ck = dk(:,7);
            
            if (length(Ck)~=dimmax) 
                continue
            end
            
            R = johansen([Ci Cj Ck], 0, 1);
            r = struct2cell(R);
            
            TS = r{3}(3);    
            ES = r{4}(3);    
            TCV = r{5}(3,3); 
            ECV = r{6}(3,3); 
            ki = r{2}(1,1);  
            kj = r{2}(2,1);  
            kk = r{2}(3,1);
            
            if (TS>TCV) && (ES>ECV)
                
                I = ki*Ci;
                J = kj*Cj;
                K = kk*Ck;
                IJK = I+J+K;
                LIJK = IJK(1:end-1);
                DIJK = IJK(2:end)-LIJK;
                LAMBDA = regress(DIJK, LIJK);
                HL = -log(2)/LAMBDA;
                
                if (LAMBDA<0) && (HL<HalfLife) && (STD<=std(IJK))
                names = strcat(char(B(i)), 'vs', char(B(j)), 'vs', char(B(k)));
                O = cell(3,2);
                O(1,1) = B(i);
                O(2,1) = B(j);
                O(3,1) = B(k);
                O{1,2} = ki;
                O{2,2} = kj;
                O{3,2} = kk;
                xlswrite(names, O);
                montakokolme = montakokolme+1;
                
                f = figure('visible', 'off');
                plot(z, IJK, z, mean(IJK), 'r', z, mean(IJK)+std(IJK), 'b', z, mean(IJK)-std(IJK), 'b');
                title(names);
                saveas(f, names, 'jpg');
                end
                
            end
            % Same for quadruplets
            for o = k+1:dim
            
                do = table2array(D{o});
                Co = do(:,7);
            
                if (length(Co)~=dimmax) 
                    continue
                end
            
                R = johansen([Ci Cj Ck Co], 0, 1);
                r = struct2cell(R);
            
                TS = r{3}(4);    
                ES = r{4}(4);    
                TCV = r{5}(4,3); 
                ECV = r{6}(4,3); 
                ki = r{2}(1,1);  
                kj = r{2}(2,1);  
                kk = r{2}(3,1);
                ko = r{2}(4,1);
            
                if (TS>TCV) && (ES>ECV)
                
                    I = ki*Ci;
                    J = kj*Cj;
                    K = kk*Ck;
                    OO = ko*Co;
                    IJKO = I+J+K+OO;
                    LIJKO = IJKO(1:end-1);
                    DIJKO = IJKO(2:end)-LIJKO;
                    LAMBDA = regress(DIJKO, LIJKO);
                    HL = -log(2)/LAMBDA;
                
                    if (LAMBDA<0) && (HL<HalfLife) && (STD<=std(IJKO))
                    names = strcat(char(B(i)), 'vs', char(B(j)), 'vs', char(B(k)), 'vs', char(B(o)));
                    O = cell(4,2);
                    O(1,1) = B(i);
                    O(2,1) = B(j);
                    O(3,1) = B(k);
                    O(4,1) = B(o);
                    O{1,2} = ki;
                    O{2,2} = kj;
                    O{3,2} = kk;
                    O{4,2} = ko;
                    xlswrite(names, O);
                    montakonelja = montakonelja+1;
                
                    f = figure('visible', 'off');
                    plot(z, IJKO, z, mean(IJKO), 'r', z, mean(IJO)+std(IJO), 'b', z, mean(IJO)-std(IJO), 'b');
                    title(names);
                    saveas(f, names, 'jpg');
                    end
                
                end
            
            end
            
        end
        
        
        
    end
    
end

% Print out number of each portfolios found
disp(strcat('Pareja löytyi',{' '}, num2str(montakoparia),'!'));
disp(strcat('Triplettejä löytyi',{' '}, num2str(montakokolme), '!'));
disp(strcat('Kvadrupletteja löytyi',{' '}, num2str(montakonelja),'!'));


figure('visible', 'on');
ll
toc
end

        
    
    


