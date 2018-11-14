function [ xi,P,predictedxi,predictedP ] = KalmanFilter( parameter_vector,y )
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior (roughly)

% Extract lenght of the data
T = size(y,2);

% Extract the stuff we need from the input arguments
F = parameter_vector(1,1);
Q = parameter_vector(2,1);
R = parameter_vector(3,1);

% The Kalman filter for AR(1)+noise
for t=1:T
    % Diffuse initialisation
    if t==1
    % Update
    xi(t)           = y(t);
    P(t)            = R;
    predictedxi(t) = y(t);
    predictedP(t)   = R;
    %predictedxi(t)  = 0;
    %predictedP(t)   = 10^8;
    % Predict
    else
    predictedxi(t)  = F * xi(t-1);
    predictedP(t)   = F * P(t-1) * F' + Q;
    % Update
    xi(t)  = predictedxi(t)  + predictedP(t) * 1/( predictedP(t) + R ) * ( y(t) - predictedxi(t) );
    P(t)   = predictedP(t)   - predictedP(t) * 1/( predictedP(t) + R ) * predictedP(t);
    end
    % Close the loop over time
end
% Close the function
end

