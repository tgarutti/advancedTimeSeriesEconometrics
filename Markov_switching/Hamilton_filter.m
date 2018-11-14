function [ filteredxi , predictedxi ] = Hamilton_filter(p11,p22,mu,sigma,y)

% Extract length of data
T = length(y);

% We can "hack" into it by forcing p11 + p22 = 1, which makes this a finite
% mixture model (of course this implies the parameter p22 is no longer
% identified)
% p22 = 1 - p11;

% Build transition matrix from p11 and p22
P   = [ p11 , 1-p22 ; 1-p11 , p22];

% Initialise at long-run equilibrium probabilities of being in state 1 or 2
predictedxi(:,1) = [ (1-p22) / (2-p11-p22) ; (1-p11)/(2-p11-p22) ];

% Run the Hamilton filter
for i=1:T
   likelihood(:,i)   = [ normpdf(y(1,i),mu(1),sigma(1)) ; normpdf(y(1,i),mu(2),sigma(2)) ];
   filteredxi(:,i)   = predictedxi(:,i) .* likelihood(:,i) / ([1,1]*(predictedxi(:,i).*likelihood(:,i)) );
   predictedxi(:,i+1)= P * filteredxi(:,i) ;
end

% Delete the last prediction, because we want filteredxi and
% predictedxi to have the same length
predictedxi = predictedxi(:,1:T);

% Close the function
end

