%% Open folder (Mac)
clear
pad = '/Users/Rutger-Jan/Desktop/KF_demo'; % change this to where you've saved this folder
cd(pad);
format short

%% Open folder (Windows)
%clear 
%pad = 'C:/Users/61849rla/Desktop/KF_demo';
%cd(pad);
%format short

%% Load data
load('CPI_AUC_SA_NSA.mat')
logreturn  = 100 * 2 * log( CPI(2:end) ./ CPI(1:end-1) )';
mu         = mean(logreturn);
y          = logreturn -mu; % demeaned log-returns
T          = length(y);
dates2     = dates(2:end)'+datenum('31/12/1899','dd/mm/yyyy'); % Correction because Excel dates start counting at 1 Jan 1900

%% Plot data 

figure
plot(dates2,mu+y,'o-','MarkerSize',20,'MarkerEdgeColor','k')
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/8/2016','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
%title('Annualised inflation (semi-annual data), $y_t=2\times \log(CPI_{t}/CPI_{t-1})$','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%% Run the filter at some arbitary values and plot the result

[xi,P,predictedxi,predictedP] = KalmanFilter([1;1;30],y);
% [1;0.5;1/100] are phi, sigma_eta^2 and sigma_epsilon^2
% Play around with the values here to see the effect on the output!

figure
plot(dates2,mu+y,'o','MarkerSize',20,'MarkerEdgeColor','k')
hold on
plot(dates2,mu+xi,'r','LineWidth',2)
hold on
plot(dates2,mu+xi+2*sqrt(P),'b-','LineWidth',1)
hold on
plot(dates2,mu+xi-2*sqrt(P),'b-','LineWidth',1)
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/1/2017','dd/mm/yyyy')
axis([startrow,endrow,-3,7])
%title('Inflation (filtered)','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')


%% Get the correct parameters by Maximimum likelihood estimation (fminunc)

% Pick starting values roughly consistent with the data
dy = y(2:end)-y(1:end-1);
startingvalues=[1;1/3*var(y);2/3*var(y)];

% Do the optimisation
clearvars options
options  =  optimset('fminunc');
[ML_parameters,ML_LogL]=fminunc('NegativeLogLikelihood',  startingvalues,options,y )

%% Get the correct parameters by Maximimum likelihood estimation (fmincon)

clearvars options
options  =  optimset('fmincon')
lb = [0;0;0];
ub = [1;2;2];
[ML_parameters,ML_LogL]=fmincon('NegativeLogLikelihood', startingvalues,[],[],[],[],lb,ub,[],options,y )


%% Calculate standard errors and display results

format short

ML_std = sqrt( diag ( inv(  fdhess6('NegativeLogLikelihood',ML_parameters,y) )));

[ML_parameters,ML_std]

%% Run the filter and smoothed at the estimated parameters

[xi,P,predictedxi,predictedP] = KalmanFilter(ML_parameters,y);

[smoothedxi,smoothedP] = KalmanSmoother(ML_parameters,y);

%% The Kalman filter inproves the RMSE vs a random walk

sqrt( var( y(2:end) - predictedxi(2:end) ) )
sqrt( var( y(2:end) - y(1:end-1)         ) )

%% This matches this expression
sqrt( mean( predictedP ) + 0.9184^2 )

%% The relative improvement is around 15%
sqrt( var( y(2:end) - predictedxi(2:end) ) ) / sqrt( var( y(2:end) - y(1:end-1) ) )

%% Plot the filtered result
figure
plot(dates2,mu+y,'o','MarkerSize',20,'MarkerEdgeColor','k')
hold on
plot(dates2,mu+xi,'r','LineWidth',2)
hold on
plot(dates2,mu+xi+2*sqrt(P),'b-','LineWidth',1)
hold on
plot(dates2,mu+xi-2*sqrt(P),'b-','LineWidth',1)
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/1/2017','dd/mm/yyyy')
axis([startrow,endrow,-3,7])
%title('Inflation (filtered)','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')


%% Plot the smoothed result
figure
plot(dates2,mu+y,'o','MarkerSize',20,'MarkerEdgeColor','k')
hold on
plot(dates2,mu+smoothedxi,'r','LineWidth',2)
hold on
plot(dates2,mu+smoothedxi+2*sqrt(smoothedP),'b-','LineWidth',1)
hold on
plot(dates2,mu+smoothedxi-2*sqrt(smoothedP),'b-','LineWidth',1)
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/1/2017','dd/mm/yyyy')
axis([startrow,endrow,-3,7])
%title('Inflation (smoothed)','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%%
smoothedP(end)

%% Plot the uncertainty
figure
plot(dates2,sqrt(P),'b','LineWidth',1)
hold on
plot(dates2,sqrt(smoothedP),'r','LineWidth',1)
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/1/2017','dd/mm/yyyy')
axis([startrow,endrow,0,1])
%title('Inflation (smoothed)','Interpreter','latex')
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
