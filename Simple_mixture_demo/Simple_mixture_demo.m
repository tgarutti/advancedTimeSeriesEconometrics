%% Open folder (Mac)
clear
pad = '/Users/Rutger-Jan/Desktop/Simple_mixture_demo'; % change this to where you've saved this folder
cd(pad);
format short

%% Open folder (Windows)
%clear 
%pad = 'C:/Users/61849rla/Desktop/Simple_mixture_demo';
%cd(pad);
%format short

%% Load data
load('CPI_AUC_SA_NSA.mat')
logreturn  = 100 * 2 * log( CPI(2:end) ./ CPI(1:end-1) )';
Dy         = logreturn(2:end) - logreturn(1:end-1);
T          = length(Dy);
dates2     = datenum('31/12/1899','dd/mm/yyyy')+dates(2:end)';

figure
plot(dates2(1:end),logreturn,'-o','MarkerSize',20,'MarkerEdgeColor','k')
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/8/2016','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
ylabel('$y_t$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
title('Inflation')
set(gca,'FontSize',30)

figure
plot(dates2(2:end),Dy,'-o','MarkerSize',20,'MarkerEdgeColor','k')
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/8/2016','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
ylabel('$y_t - y_{t-1}$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
title('Is inflation a random walk?')
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%% Pick starting values based (somewhat) on the data
p         = [0.9;0.1];
mu        = [mean(Dy,2);mean(Dy,2)];
sigma     = [0.9*std(Dy);1.1*std(Dy)];

%% True (numerical) likelihood: sum (over the data) of a log of a sum (over states)
L = @(v) (-1) * sum( log( abs(v(1)) * normpdf( Dy(1,:) , v(2), abs(v(4)) ) + abs(( 1 - v(1) )) * normpdf( Dy(1,:) , v(3), abs(v(5)) ) ) );

%% This maximisation goes wrong (fminunc)

startingvals=[p(1),mu(1),mu(2),sigma(1),sigma(2)]';
L(startingvals)

%%
[ML_parameters,ML_LogL]=fminunc(L, startingvals )

%% Here we force the estimates to stay within bounds (fmincon)

startingvals=[p(1),mu(1),mu(2),sigma(1),sigma(2)]';
clearvars options
options  =  optimset('fmincon');
lb = [0;-100;-100;0;0];
ub = [1;100;100;100;100];
[ML_parameters,ML_LogL]=fmincon(L,startingvals,[],[],[],[],lb,ub,[],options)

% We have forced the parameter estimates to behave, but the likelihood is worse!
% Why is this?

%% Set the iterations at zero
iteration = 0 ;
T = length(Dy) ;

%% Set the initial group probabilities randomly
%for i=1:T
%    pstar(1,i)  = rand();
%    pstar(2,i)  = 1-pstar(1,i);
%end

%% Starting values
p         = [0.8;0.2];
mu        = [mean(Dy,2);mean(Dy,2)];
sigma     = [0.95*std(Dy);1.05*std(Dy)];
iteration = 0 ;
T = length(Dy) ;
clear EM_estimates
        
%% Now let's try EM

for k=1:100 % there is no k in this code, this just runs this section 100 times
    
    % Increase the iteration counter
    iteration = iteration + 1;
    
    % Expectation step
    for i=1:T
        Pjoint1(i)  = p(1) * normpdf( Dy(i), mu(1), sigma(1) ) ;
        Pjoint2(i)  = p(2) * normpdf( Dy(i), mu(2), sigma(2) ) ;
        pstar(1,i)  = Pjoint1(i) / (Pjoint1(i)+Pjoint2(i)) ; 
        pstar(2,i)  = Pjoint2(i) / (Pjoint1(i)+Pjoint2(i)) ;
    end
    
    % Maximisation step
    p     = [ mean(pstar(1,:)) ; mean(pstar(2,:)) ] ;
    mu    = [ sum( pstar(1,:) .* Dy ) / sum( pstar(1,:) ) ; sum( pstar(2,:) .* Dy ) / sum( pstar(2,:) ) ];
    sigma = sqrt([ sum( pstar(1,:) .* ( Dy - mu(1) ).^2 ) / sum( pstar(1,:) ) ; sum( pstar(2,:) .* ( Dy - mu(2) ).^2 ) / sum( pstar(2,:) ) ]);

    % Save the parameter vector
    EM_estimates(:,iteration) = [p(1); mu(1,1);mu(2,1);sigma(1,1);sigma(2,1)];

    % Close loop over time
end

%% Convergence of EM estimates

plot(EM_estimates','LineWidth',2)
legend('p','mu1','mu2','sigma1','sigma2')
title('Convergence of EM estimates')
xlabel('EM iterations')
ylabel('Parameter estimates')
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')

%% Comparsion of EM and ML ----> very similar!

[ML_parameters,[p(1);mu(1);mu(2);sigma(1);sigma(2)]]

%% Plot the probability of high volatility state
figure
plot(dates2(2:end),Dy,'o','MarkerSize',20,'MarkerEdgeColor','k')
hold on
plot(dates2(2:end),pstar(2,:),'r','LineWidth',2)
%hold on
%plot(dates2(1:end),pstar_employment(2,:),'g','LineWidth',2)
dateFormat = 'mmmyy';
datetick('x',dateFormat)
startrow=datenum('1/1/1985','dd/mm/yyyy')
endrow=datenum('1/8/2016','dd/mm/yyyy')
axis([startrow,endrow,-inf,inf])
ylabel('$\%$','FontSize',24,'FontWeight','bold','Color','k','Interpreter','latex') % y-axis label
set(gca,'FontSize',30)
set(gca,'FontName','Times New Roman')
