%% Parameters
rho = 0.02; b = 0.5; mu =0.05; sigma = 0.010;

%% Other parameters (FRED: real annual time series of GDP 1996 - 2024)
% rho = 0.02; b = 0.5; mu =0.008; sigma = 0.036;%RUSSIA
% rho = 0.02; b = 0.5; mu =0.02958; sigma = 0.039272;%PERU 
% rho = 0.02; b = 0.5; mu =0.049793; sigma = 0.032046; %CHINA
%  rho = 0.02; b = 0.5; mu =0.00879; sigma = 0.011444; %BRAZIL

%% Calculating Bond Price and Bond yield
tauVector = 1:50;
YtVector  = 1:1:100;
B      = zeros(length(YtVector),length(tauVector));
yt_tau = zeros(length(YtVector),length(tauVector));

for i = 1:length(YtVector)
    Yt = YtVector(i);
    for j = 1:length(tauVector)
        tau = tauVector(j);
        % Bond price per maturity
        B(i,j) = B_fun(tau,rho,b,mu,sigma,Yt);
        % Bond yields
        yt_tau(i,j) = -(1/tau)*log( B(i,j) );
    end
end

%% Graph: Bond yield
BondYield = yt_tau;
% (1)Slope
    %recessions (low "w"): lowest decile: 101/10 = 10.1 (first 50 rows)
    EY_rec = mean(BondYield(1:10,:));
        %Q1
        EY_Q1 = mean(BondYield(1:25,:));    
    %normal times (high Y): highest decile: 101/10 = 10.1 (last 50 rows)
    EY_nt = mean(BondYield(end-10+1:end,:));
        %Q4
        EY_Q4 = mean(BondYield(75:end,:));    
    %Full sample
    EY_average = mean(BondYield);

        figure('Name','Bond Yields')
    subplot(2,2,1)
        plot(tauVector(1:50),100*EY_rec(1:50),'r',...
            tauVector(1:50),100*EY_average(1:50),'k--',...
            tauVector(1:50),100*EY_nt(1:50),'b','LineWidth',2.5)
        leg1 =legend('Low $Y_t$','Full sample','High $Y_t$','Location','southeast')
            xlabel('Maturity ($\tau$)','Interpreter','latex')
            ylabel('Bond Yield (%)')
            title('Bond Yield')
            grid; 
            set(leg1,'Box','off')
            ax = gca; 
            ax.FontSize = 14; 


% Save the figure
hfig=gcf;
set(hfig,'PaperOrientation','landscape');
set(hfig,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(hfig, '-dpdf', strcat('P','_Fig1.pdf'));    
