clc; clear vars; 

Fs = 1000;
 
t = 0:1/Fs:600;
f0 = 8;
f1 = 23;
f2 = 2;
f3 = 4;
chs = 20;



x = 5.* sin(2*pi*f0*t) + sin(2*pi*f1*t);
y = cos(2*pi*f2*t) + 5 .* cos(2*pi*f3*t);

% Trajectories
xy = [x',y'];

% Create some observed signals from which we will derive some features
Obs = zeros(size(x,2),chs);
cc = chs/2;
Obs(:,1:cc) = (x .* cos(2*pi*330*t))' + 1.* randn(length(t),cc);
Obs(:,cc+1:size(Obs,2)) = (y .* cos(2*pi*430*t))' + 1.* randn(length(t),size(Obs,2)-cc);


%% Generate some spike count data: FEATURES
spks = cell(1,size(Obs,2));
spkT = cell(1,size(Obs,2));

for i=1:size(Obs,2)
% [spks{i},spkT{i}] = findpeaks(Obs(:,i),Fs,'MinPeakProminence',5);
[spks{i},spkT{i}] = findpeaks(Obs(:,i),Fs);
end

% Count spikes in time bins
deltaTime = 0.1;

nBins = round((size(Obs,1)/Fs) / deltaTime);

% These are the spike counts in deltaTime bins
feats = zeros(nBins,chs);
Bin=1;
for t=0:deltaTime:(size(Obs,1)-2)/Fs

    for ch=1:chs
            tp = length(find(spkT{ch}>t & spkT{ch}<=t+deltaTime));
            feats(Bin,ch) = tp;      
    end
    Bin = Bin+1;
end



%% Identify training data and cursor trajectory data

% Let's start by using trainingTime seconds worth of data
trainingTime = 400;
% Downsample ideal signals to match the time bin width selected for the
% features
ds = Fs*deltaTime;
X_t = downsample(xy,ds);
nTsamps = round(trainingTime * (1/deltaTime));

Xtraining = X_t(1:nTsamps,:);

%nTsamps = trainingTime * Fs;

%Xtr = xy(1:nTsamps,:);

%tdown = 0:deltaTime:(length(X_t)-1)*deltaTime;
%plot(tdown,X_t(:,1))
%figure
%plot(0:1/Fs:(length(x)-1)/Fs,x)

%% Traditional Kalman Filter approach
%%
%%%%***     Build the Kalman Filter Model 
%%%%***     X_t = A*X_(t-1) + W_t;
%%%%***     Y_t = C*X_t + Q_t;

% Full state matrix
Vxy = [[0,0];diff(X_t)./deltaTime];
Axy = [[0,0];diff(Vxy)./deltaTime];
x_t = [X_t,Vxy,Axy];

%x_t = [X_t,Vxy,ones(length(X_t),1)]; % This is suggested by Chestek for
% reducing baseline variations

% Let's try to detrend the features
% It turns out this is probably fundamental in catching the variability and
% thus the correlation between the features and the trajectory signal
feats = detrend(feats);

% Train the filter
fakeKalman = fitKalmanContinuous(x_t(1:nTsamps,:)',feats(1:nTsamps,:)',nTsamps);

% Run the filter simulating RT processing (one time sample at a time)
len = size(x_t,1)-nTsamps;
predictedXY = zeros(size(x_t,2),len);

for i=1:len-1
    fakeKalman = rtmyKalman(fakeKalman,feats(nTsamps+i,:));
    %predictedXY(1,i) = fakeKalman.x_apo(1);
    %predictedXY(2,i) = fakeKalman.x_apo(2);
    predictedXY(:,i) = fakeKalman.x_apo;
end


%% Plot the results
%%% FIGURE 1 shows you the difference between predicted and ideal
%%% trajectories
predictionTime = nTsamps*deltaTime:deltaTime:(length(x_t)-1)*deltaTime;
figure(1)
subplot(2,1,1)
plot(predictionTime,predictedXY(1,:));
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),1))
legend('Predicted X','Ideal X')
title('X position')
subplot(2,1,2)
plot(predictionTime,predictedXY(2,:));
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),2))
legend('Predicted Y','Ideal Y')
title('Y position')
xlabel('TIME (s)')


%% Plot the raw data used for feature extraction
%%% FIGURE 2 shows you how the features (only 2 channels used) compare to
%%% the trajectory signals
figure(2)
subplot(2,1,1)
plot(predictionTime,[0;feats(nTsamps+1:length(feats),1)])
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),1))
legend('Feature 1','Ideal X')
title('X position')
subplot(2,1,2)
plot(predictionTime,[0;feats(nTsamps+1:length(feats),size(feats,2))])
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),2))
legend('Feature Last','Ideal Y')
title('Y position')

%% Plot Velocities and their predictions
%% Plot the results
%%% FIGURE 3 shows you the difference between predicted velocities and ideal
%%% velocities
predictionTime = nTsamps*deltaTime:deltaTime:(length(x_t)-1)*deltaTime;
figure(3)
subplot(2,1,1)
plot(predictionTime,predictedXY(3,:));
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),3))
legend('Predicted X','Ideal X')
title('X Velocity')
subplot(2,1,2)
plot(predictionTime,predictedXY(4,:));
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),4))
legend('Predicted Y','Ideal Y')
title('Y Velocity')
xlabel('TIME (s)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LINEAR MODEL %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Train the linear model
myLinearMod = fitLinearContinuous(X_t(1:nTsamps,:)',feats(1:nTsamps,:)');

% Run the predictions on the testing data
XXX = predict(myLinearMod.X,feats(nTsamps+1:size(feats,1),:));
YYY = predict(myLinearMod.Y,feats(nTsamps+1:size(feats,1),:));

%%% Plot
%%% Figure 4
predictionTime = nTsamps*deltaTime:deltaTime:(length(x_t)-1)*deltaTime;
figure(4)
subplot(2,1,1)
plot(predictionTime,[0;XXX]);
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),1))
legend('Predicted X','Ideal X')
title('X Position')
subplot(2,1,2)
plot(predictionTime,[0;YYY]);
hold on
plot(predictionTime,x_t(nTsamps+1:length(x_t),2))
legend('Predicted Y','Ideal Y')
title('Y Position')
xlabel('TIME (s)')