%% Clear
clear
pad = '/Users/Rutger-Jan/Desktop/Markov_Switching' 
% Change the above line to wherever you have saved this folder
cd(pad)

%% 'True' constants to generate the data (see slide 13, Lecture 4)

T      = 10^3; % length of the dataset
p11    = 0.99; % persistence of state 1
p22    = 0.95; % persistence of state 2
P      = [ p11 , 1-p22 ; 1-p11 , p22]; % transition matrix
mu     = [0;0]; % the observations have mean zero in either state
sigma  = [1/3;2/3]; % the standard dev in state 2 is twice as big
xi(:,1)= [1;0]; % we start off in state 1

%% Generate state vector xi(:,t) for each t
for i=2:T
    stay_in_1           = heaviside( p11 - rand ); % this returns a '1' if we remain in state 1
    stay_in_2           = heaviside( p22 - rand ); % this returns a '1' if we remain in state 2
    randomtransitionmatrix = [ stay_in_1 , 1-stay_in_2 ; 1-stay_in_1 , stay_in_2 ]
    xi(:,i)             = randomtransitionmatrix * xi(:,i-1); % converts previous state into next state
end

% The data that you generate will be different from that in Lecture 4!

figure
plot(xi(2,:),'r','Linewidth',2)
axis([-inf inf -1/10 11/10])
hold off
%ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
title('Simulated high volatility state')
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%% Generate noisy signals of the state, i.e. realisations of y_t
for i=1:T
    if xi(1,i)>0 
        state = 1;
    else
        state = 2;
    end
    y(i) = normrnd( mu(state), sigma(state) ) ;
    % if the state is i, we draw an observation from normrnd(mu(i),sigma(i))
end

%% Display the data together with the state

figure
plot(y,'k','Linewidth',0.3)
hold on
plot(xi(2,:),'r','Linewidth',3)
axis([0 500 -inf inf])
hold off
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
legend('Generated data','High volatility state')

%% Run the Hamilton filter

[ filteredxi , predictedxi ] = Hamilton_filter(p11,p22,mu,sigma,y);

figure
plot(y,'k','Linewidth',0.1)
hold on
plot(xi(2,:),'r','Linewidth',3)
hold on
plot(filteredxi(2,:),'g','Linewidth',2)
axis([-inf inf -11/10 11/10])
hold off
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
legend('Data','High volatility state','Filtered probability')

%% In general, we do not know the parameters. So we maximise the log likelihood

% Starting values based roughly on the data
% starting values = [p_{11},p_{22},[mu_1;mu_2],[sigma_1;sigma_2] ]
startingvalues = [0.8;0.8;[mean(y);mean(y)];[1/2*std(y);2*std(y)]]

% It is useful to take the starting points *not* symmetric in sigma_1 and
% sigma_2. Also, initialise the p_{11} and p_{22} as relatively persistent

% Get a reference point for the 
LogLikelihood(startingvalues,y)

clearvars options
options  =  optimset('fminunc'); % This sets the options at their standard values
options  =  optimset(options , 'MaxFunEvals' , 10^10) ; % extra iterations
options  =  optimset(options , 'TolFun'      , 1e-10); % extra precision
options  =  optimset(options , 'TolX'        , 1e-10); % extra precision

% Run the ML optimisation and store the result:
[ML_parameters,ML_LogL]=fminunc('LogLikelihood',startingvalues,options,y)

%% Optimise using fmincon

clearvars options
options  =  optimset('fmincon') % This sets the options at their standard values
options  =  optimset(options , 'MaxFunEvals' , 10^10) ; % extra iterations
options  =  optimset(options , 'TolFun'      , 1e-10); % extra precision
options  =  optimset(options , 'TolX'        , 1e-10); % extra precision

lb = [0;0;-100;-100;0;0]; % Lower bound for the parameter vector [p_{11},p_{22},mu_1,mu_2,sigma_1,sigma_2]
ub = [1;1;100;100;100;100]; % Upper bound 

[ML_parameters,ML_LogL]=fmincon('LogLikelihood', startingvalues,[],[],[],[],lb,ub,[],options,y )

%% Calculate standard errors and display results

% Standard errors are calculated by calculating the Hessian of the
% LogLikelihood at the ML parameters. Then invert, take diagional, square
% root.

format short
ML_std = sqrt( diag ( inv(  fdhess6('LogLikelihood',ML_parameters,y) )))

% display true parameters, estimated, and standard error
[[p11;p22;mu;sigma],ML_parameters,ML_std]

%% Now run the filter at the estimated parameters

p11   = ML_parameters(1,1);
p22   = ML_parameters(2,1);
mu    = ML_parameters(3:4,1);
sigma = ML_parameters(5:6,1);

[ filteredxi_ML , predictedxi_ML ] = Hamilton_filter(p11,p22,mu,sigma,y);

figure
plot(y,'k','Linewidth',0.1)
hold on
plot(xi(2,:),'r','Linewidth',3)
hold on
plot(filteredxi(2,:),'b','Linewidth',2)
hold on
plot(filteredxi_ML(2,:),'g-','Linewidth',2)
axis([-inf inf -11/10 11/10])
hold off
%ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
legend('Data','High volatility state','Filtered probability (true parameters)','Filtered probability (ML parameters)')

%% Run smoother at ML parameters

p11   = ML_parameters(1,1);
p22   = ML_parameters(2,1);
mu    = ML_parameters(3:4,1);
sigma = ML_parameters(5:6,1);

[ smoothedxi_ML ] = Hamilton_smoother(p11,p22,mu,sigma,y);

figure
plot(y,'k','Linewidth',0.1)
hold on
plot(xi(2,:),'r','Linewidth',3)
hold on
plot(filteredxi_ML(2,:),'g-','Linewidth',2)
hold on
plot(smoothedxi_ML(2,:),'b','Linewidth',2)
axis([-inf inf -11/10 11/10])
hold off
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
legend('Data','High volatility state','Filtered probability (ML parameters)','Smoothed probability (ML parameters)')

