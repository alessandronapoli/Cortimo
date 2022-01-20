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

predictions = zeros(size(X_t));

newTrainedKalman = fitKalmanContinuous(X_t(:,1:nS), Y_tr, nS);

updatedKalman = rtmyKalman(newTrainedKalman, feats(nS,:));
predictions(:,1) = updatedKalman.x_apo;
j = 2;
for i = nS + 1:length(feats(:,1))
    updatedKalman = rtmyKalman(updatedKalman, feats(i,:));
    predictions(:,j) = updatedKalman.x_apo;
    j = j+1;
end

