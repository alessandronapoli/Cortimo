clearvars; close all; clc;

Fs = 1e3;
t = 0:1/Fs:600;
f0 = 30;
f1 = 5;
nChs = 10;

% idealSig = 17*sin(2*pi*f0*t) + 5.*sin(2*pi*f1*t) ;
signal1 = 17.*randn(1,length(t)) + 17*sin(2*pi*f0*t) + 15*randn(1,length(t)) .* 3.*sin(2*pi*f1*t) ;

rawData = 5.*randn(nChs,length(signal1)) .* ones(nChs,length(signal1)) .* signal1;
rawData = rawData'; % Column matrix. rows are datapoints, columns are channels

% plot(t,signal1)
% hold on
% plot(t,idealSig)
% legend('Noisy','Ideal')


%% Generate some spike count data
spks = cell(1,size(rawData,2));
spkT = cell(1,size(rawData,2));

for i=1:size(rawData,2)
[spks{i},spkT{i}] = findpeaks(rawData(:,i),Fs,'MinPeakProminence',50);
end

% Count spikes in time bins
deltaTime = 0.05;

nBins = round((size(rawData,1)/Fs) / deltaTime);
feats = zeros(nBins,nChs);
Bin=1;
for t=0:deltaTime:(size(rawData,1)-2)/Fs

    for ch=1:nChs
            tp = length(find(spkT{ch}>t & spkT{ch}<=t+deltaTime));
            feats(Bin,ch) = tp;      
    end
    Bin = Bin+1;
end

%% Identify training data and cursor trajectory data


% Let's say my x position correlates well with signal1
% Let's say my y position correlates well with signal1 in phase quadrature
% What is to be expected is a back and forth motion along a line

% Let's start by using trainingTime seconds worth of data
trainingTime = 300;
nTsamps = trainingTime * Fs;
% X = randn(length(signal1),1) + signal1';
% Y = randn(length(signal1),1) + cos(signal1)';

% Build a slower sham trajectory
% X = lowpass(randn(length(signal1),1) + signal1',10,Fs);
% Y = lowpass(randn(length(signal1),1) + cos(signal1)',10,Fs);

% X = lowpass(randn(length(signal1),1) + (cos(2*pi*t*5) .* signal1'),3,Fs);
% Y = lowpass(randn(length(signal1),1) + (cos(2*pi*t*3) .* cos(signal1)'),3,Fs);
X = lowpass(randn(length(signal1),1) + 3.*(cos(2*pi*t*2) .* signal1'),3,Fs);
Y = lowpass(randn(length(signal1),1) + 3.*(cos(2*pi*t*1) .* cos(signal1)'),3,Fs);


Xtr = X(1:nTsamps);
% Ytr = randn(nTsamps,1) - signal1(1:nTsamps)' + cos(signal1(1:nTsamps))';
Ytr = Y(1:nTsamps);

training_length = nTsamps;

%%
%%
%%
%% Traditional Kalman Filter approach
%%

%%%%***     Build the Kalman Filter Model 
%%%%***     X_t = A*X_(t-1) + W_t;
%%%%***     Y_t = C*X_t + Q_t;


% Kinematic state of the Kalman model
% x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs],ones(size(Xtr,1),1)];
x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs]];

% The measurements in this specific case only occur every time bin deltaTime
% Before using the system we should bring all the variables to the same
% sampling frequency, that is the trajectory information can be
% downsampled.
downsampFactor = deltaTime / (1/Fs);
% Downsampled now
X_tPOS = downsample(x_t(:,1:2),downsampFactor)';
X_t = [X_tPOS;[0,diff(X_tPOS(1,:)).*deltaTime];[0,diff(X_tPOS(2,:)).*deltaTime]];

% Measurmement Matrix
nS = (trainingTime * (Fs/downsampFactor)); % nS time bins correspond to training data
Y_tr = feats(1:nS,:)';

% These are one-step shifted state matrices
X_1 = X_t(:,1:end-1);
X_2 = X_t(:,2:end);

% Calculate the coefficient matrices
% H = (Z_t*X_t')*inv(X_t*X_t'); % Slightly different nomenclature used elsewhere
% C = (Y_t * X_t') * (inv(X_t * X_t'));
C = (Y_tr * X_t') * (pinv(X_t * X_t'));

Q = (1/nS) .* ((Y_tr - (C * X_t)) * (Y_tr - (C * X_t))'); % Q = (Z_t-H*X_t)*(Z_t-H*X_t)' ./ (training_length);

% A = (X_2 * X_1') * inv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');
A = (X_2 * X_1') * pinv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');

W = (1/(nS-1)) .* (X_2 - A * X_1) * ((X_2 - A * X_1)');
% W = ((X_2-A*X_1)*(X_2-A*X_1)')./(training_length-1);

%% Let's run the prediction now
% Pk_apr = zeros(size(X_t,1),size(X_t,1));
% Pk_apo = ones(size(X_t,1),size(X_t,1));
Pk_apo = diag(10.*ones(1,size(X_t,1)));
x_apr = zeros(size(X_t,1),size(feats,1));
x_apo = zeros(size(X_t,1),size(feats,1));

% Generate the complete vector observations
Y_t = feats';

% x_apo(:,nS) = X_t(:,nS);

for step = nS+1:size(feats,1)
    
    % TIME UPDATE EQUATIONS
    % Step 1
    x_apr(:,step) = A*x_apo(:,step-1);  % Estimate from previous time point
    
    % Step 2
    Pk_apr = A*Pk_apo*A' + W;  % Covariance matrix of estimated error
    
    % MEASUREMENT UPDATE EQUATIONS
    try
    % Step 3
    % KK = Pk_apr*H'*inv(H*Pk_apr*H' + Q); % If using different notation
    % KK = Pk_apr*C'*inv(C*Pk_apr*C' + Q); % Filter Gain
    KK = Pk_apr*C'*pinv(C*Pk_apr*C' + Q); % Filter Gain
    
    
    % Step 4
    % x_apo(step) = x_apr(step) + KK(step)*(z(step-H*x_apr(step)));
    x_apo(:,step) = x_apr(:,step) + KK*(Y_t(:,step)-C*x_apr(:,step));  % A posteriori estimate

    % Step 5
    Pk_apo = (ones(size(X_t,1),size(X_t,1))-KK*C)*Pk_apr;
    
   catch
        Pk_apo = diag(1.*ones(1,size(X_t,1)));
        x_apo(:,step) = x_apr(:,step);

    end
    
end


% x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs],ones(size(Xtr,1),1)];
% x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs]];

% The measurements in this specific case only occur every time bin deltaTime
% Before using the system we should bring all the variables to the same
% sampling frequency, that is the trajectory information can be
% downsampled.
downsampFactor = deltaTime / (1/Fs);

% Downsampled now
XPOS = downsample(X,downsampFactor)';
YPOS = downsample(Y,downsampFactor)';

% ideal_X_t = [X,Y,[0;diff(X).*Fs],[0;diff(Y).*Fs],ones(size(X,1),1)];
% ideal_X_t = [X,Y,[0;diff(X).*Fs],[0;diff(Y).*Fs]];
iX_t = [XPOS;YPOS;[0,diff(XPOS).*deltaTime];[0,diff(YPOS).*deltaTime]];

% iX_t = downsample(ideal_X_t,downsampFactor)';
iX_t = iX_t(:,1:end-1);

%% Plot Kalman results
figure
tt = (0:(size(iX_t,2)-1))*deltaTime;
subplot(2,1,1)
plot(tt,iX_t(1,:));
hold on
plot(tt,x_apo(1,:));
legend('Ideal','Rebuilt');

subplot(2,1,2)
plot(tt,iX_t(2,:));
hold on
plot(tt,x_apo(2,:));
legend('Ideal','Rebuilt');


%%
%%
%%
%% Linear Model
%%
% linearmodX = fitlm(X_t(1,:),Y_tr);
% linearmodY = fitlm(X_t(2,:),Y_tr);
linearmodX = fitlm(Y_tr',X_t(1,:)','linear');
linearmodY = fitlm(Y_tr',X_t(2,:)','linear');

% Predict
XXX = predict(linearmodX,Y_t(:,nS+1:size(Y_t,2))');
YYY = predict(linearmodY,Y_t(:,nS+1:size(Y_t,2))');

% Visualize linear regression results
figure
subplot(2,1,1)
plot(tt,iX_t(1,:));
hold on
plot(tt,[X_t(1,:),XXX']);
legend('Ideal','Rebuilt');

subplot(2,1,2)
plot(tt,iX_t(2,:));
hold on
plot(tt,[X_t(2,:),YYY']);
legend('Ideal','Rebuilt');


%%
%%
%%
%% Velocity only Kalman filter
%%

%%%%***     Build the Kalman Filter Model 
%%%%***     X_t = A*X_(t-1) + W_t;
%%%%***     Y_t = C*X_t + Q_t;


% Kinematic state of the Kalman model
x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs],ones(size(Xtr,1),1)];
% The constant 1s are added to the matrix to have a fixed offset, which
% corresponds to the baseline firing rate.

downsampFactor = deltaTime / (1/Fs);
% Downsampled now
X_tPOS = downsample(x_t(:,1:2),downsampFactor)';
X_t = [X_tPOS;[0,diff(X_tPOS(1,:)).*deltaTime];[0,diff(X_tPOS(2,:)).*deltaTime];ones(1,size(X_tPOS,2))];

% Measurmement Matrix
nS = (trainingTime * (Fs/downsampFactor)); % nS time bins correspond to training data
Y_tr = feats(1:nS,:)';

% These are one-step shifted state matrices
X_1 = X_t(:,1:end-1);
X_2 = X_t(:,2:end);


% Calculate the coefficient matrices
% H = (Z_t*X_t')*inv(X_t*X_t'); % Slightly different nomenclature used elsewhere
% C = (Y_t * X_t') * (inv(X_t * X_t'));
C = (Y_tr * X_t') * (pinv(X_t * X_t'));

Q = (1/nS) .* ((Y_tr - (C * X_t)) * (Y_tr - (C * X_t))'); % Q = (Z_t-H*X_t)*(Z_t-H*X_t)' ./ (training_length);

% A = (X_2 * X_1') * inv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');
AComplete = (X_2 * X_1') * pinv(X_1 * X_1'); %  A = (X_2*X_1')*inv(X_1*X_1');

% W = ((X_2-A*X_1)*(X_2-A*X_1)')./(training_length-1);

% So far there is no difference from the traditional Kalman approach. Now
% let's introduce some constraints on the form of A and W matrices using
% some simple kinematic considerations, that is integrated velocity
% precisely gives position

A = diag(ones(1,size(x_t,2)));

A(1,3) = deltaTime;
A(2,4) = deltaTime;
A(3,3) = AComplete(3,3);
A(3,4) = AComplete(3,4);
A(4,3) = AComplete(4,3);
A(4,4) = AComplete(4,4);

WComplete = (1/(nS-1)) .* (X_2 - A * X_1) * ((X_2 - A * X_1)');
W = zeros(size(WComplete));
W(3,3) = WComplete(3,3);
W(3,4) = WComplete(3,4);
W(4,3) = WComplete(4,3);
W(4,4) = WComplete(4,4);


%% Let's run the prediction now
% Pk_apr = zeros(size(X_t,1),size(X_t,1));
% Pk_apo = ones(size(X_t,1),size(X_t,1));
Pk_apo = diag(1.*ones(1,size(X_t,1)));
x_apr = zeros(size(X_t,1),size(feats,1));
x_apo = zeros(size(X_t,1),size(feats,1));

% Generate the complete vector observations
Y_t = feats';

% x_apo(:,nS) = X_t(:,nS);

for step = nS+1:size(feats,1)
    % TIME UPDATE EQUATIONS
    % Step 1
    x_apr(:,step) = A*x_apo(:,step-1);  % Estimate from previous time point
    
    % Step 2
    Pk_apr = A*Pk_apo*A' + W;  % Covariance matrix of estimated error
    
    % MEASUREMENT UPDATE EQUATIONS
    try
    % Step 3
    % KK = Pk_apr*H'*inv(H*Pk_apr*H' + Q); % If using different notation
    % KK = Pk_apr*C'*inv(C*Pk_apr*C' + Q); % Filter Gain
    KK = Pk_apr*C'*pinv(C*Pk_apr*C' + Q); % Filter Gain
    
    % Step 4
    % x_apo(step) = x_apr(step) + KK(step)*(z(step-H*x_apr(step)));
    x_apo(:,step) = x_apr(:,step) + KK*(Y_t(:,step)-C*x_apr(:,step));  % A posteriori estimate

    % Step 5
    Pk_apo = (ones(size(X_t,1),size(X_t,1))-KK*C)*Pk_apr;

    catch
        Pk_apo = diag(1.*ones(1,size(X_t,1)));
        x_apo(:,step) = x_apr(:,step);
    end
end


% x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs],ones(size(Xtr,1),1)];
% x_t = [Xtr,Ytr,[0;diff(Xtr).*Fs],[0;diff(Ytr).*Fs]];

% The measurements in this specific case only occur every time bin deltaTime
% Before using the system we should bring all the variables to the same
% sampling frequency, that is the trajectory information can be
% downsampled.
downsampFactor = deltaTime / (1/Fs);

% Downsampled now
XPOS = downsample(X,downsampFactor)';
YPOS = downsample(Y,downsampFactor)';

% ideal_X_t = [X,Y,[0;diff(X).*Fs],[0;diff(Y).*Fs],ones(size(X,1),1)];
% ideal_X_t = [X,Y,[0;diff(X).*Fs],[0;diff(Y).*Fs]];
iX_t = [XPOS;YPOS;[0,diff(XPOS).*deltaTime];[0,diff(YPOS).*deltaTime]];

% iX_t = downsample(ideal_X_t,downsampFactor)';
iX_t = iX_t(:,1:end-1);

%% Plot Kalman results
figure
tt = (0:(size(iX_t,2)-1))*deltaTime;
subplot(2,1,1)
plot(tt,iX_t(3,:));
hold on
plot(tt,x_apo(3,:));
legend('Ideal','Rebuilt');

subplot(2,1,2)
plot(tt,iX_t(4,:));
hold on
plot(tt,x_apo(4,:));
legend('Ideal','Rebuilt');
title('VELOCITY')

