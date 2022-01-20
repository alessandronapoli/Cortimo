clearvars; close all; clc;

Fs = 1e3;
t = 0:1/Fs:6;
f0 = 30;
f1 = 5;

idealSig = 17*sin(2*pi*f0*t) + 5.*sin(2*pi*f1*t) ;
signal1 = 17.*randn(1,length(t)) + 17*sin(2*pi*f0*t) + 15*randn(1,length(t)) .* 3.*sin(2*pi*f1*t) ;
 
plot(t,signal1)
hold on
plot(t,idealSig)
legend('Noisy','Ideal')


%% Create some correlated features

training_length = 200;
numFeatures = 20;
num_observations = 1;

Z = signal1 + 30*randn(numFeatures,length(t));

% Using the Kalman Filter Model
X_t = signal1(:,1:training_length);
X_1 = X_t(:,1:end-1);
X_2 = X_t(:,2:end);
Z_t = Z(1:numFeatures,1:training_length);

% Calculate the 4 necessary parameters for the Kalman
A = (X_2*X_1')*inv(X_1*X_1');
H = (Z_t*X_t')*inv(X_t*X_t');
W = ((X_2-A*X_1)*(X_2-A*X_1)')./(training_length-1);
Q = (Z_t-H*X_t)*(Z_t-H*X_t)' ./ (training_length);

%% Let's run the prediction now
x_apr = zeros(size(signal1));
x_apo = zeros(size(signal1));
%Pk_apr = zeros(3,3);
Pk_apr = zeros(num_observations,num_observations);
%Pk_apo = zeros(3,3);
Pk_apo = zeros(num_observations,num_observations);
KK = zeros(size(signal1));

z = Z;

for step=training_length:length(t)-1
    % TIME UPDATE EQUATIONS
    % Step 1
    x_apr(:,step) = A*x_apo(:,step-1);  % Estimate from previous time point
    
    % Step 2
    % Pk_apr(step) = A*Pk_apo(step-1)*A' + W; % Covariance matrix of the estimated error
    Pk_apr = A*Pk_apo*A' + W; % Covariance matrix of the estimated error
    
    % MEASUREMENT UPDATE EQUATIONS
    % Step 3
    % Filter gain
    %KK(step) = Pk_apr(step)*H'*inv(H*Pk_apr(step)*H' + Q); 
    KK = Pk_apr*H'*inv(H*Pk_apr*H' + Q);
    %%%%% KK = Pk_apr*H'/(H*Pk_apr*H' + Q);
    
    
    % A posteriore estimate
    % x_apo(step) = x_apr(step) + KK(step)*(z(step-H*x_apr(step)));
    x_apo(:,step) = x_apr(:,step) + KK*(z(:,step)-H*x_apr(:,step));
    
    % Step 4
    % Calculate the covariance matrix of the error estimate
   % Pk_apo = (ones()-KK(step)*H)*Pk_apr(step);
   % Pk_apo = (ones(3,3)-KK*H)*Pk_apr;
    Pk_apo = (ones(num_observations,num_observations)-KK*H)*Pk_apr;

end

%% Plot Kalman results
figure
plot(t,x_apo)
hold on
plot(t,idealSig)
plot(t,signal1)
legend('Rebuilt','Ideal','Observed')


%% 




