clearvars; close all; clc;

% Let's implement a Kalman Filter approach like the one described by the
% Ohio group paper. Details are also summarized in the Signal Prcessing
% Guidelines.docx file.

% Electrode number
en = 192;
% Acuisition system frequency
Fs = 10e3;
% BCI analysis time bin
BCI_win = 0.02; 
% Design the passband
Hd1 = mypassband(Fs,250,3000);

duration = 60;
% 3 joints moving linearly back and forth
% Position is normalized for each joint ROM
motiontime= 0:BCI_win:duration;
j1 = (sin(2*pi*0.25*motiontime)./2+0.5)*100;
j2 = (sin(2*pi*0.125*motiontime)./2+0.5)*100;
j3 = (sin(2*pi*0.5*motiontime)./2+0.5)*100;

motion = [j1;j2;j3];

actime = 0:1/Fs:duration;
sig = 80 .* randn(192,length(actime));
sig(1:5:end,:) = sig(1:5:end,:) + 50.*sin(2*pi*actime*3000);
sig(1:30:end,:) = sig(1:30:end,:) + 50.*sin(2*pi*actime*0.25);
sig(1:22:end,:) = sig(1:22:end,:) + 50.*sin(2*pi*actime*0.5);
sig(1:46:end,:) = sig(1:46:end,:) + 50.*sin(2*pi*actime*0.125);

filtsig = zeros(size(sig));
% Filter
for i=1:size(sig,1)
    filtsig(i,:) = filtfilt(Hd1.sosMatrix,Hd1.scaleValues,sig(i,:));
end

%% Extract firing Rate
% For each channel find the threshold and count spikes
th = rms(sig,2) * 1.5;

%sp = zeros(192,1);
% Get the spikes in terms of samples
sp = cell(size(sig,1),1);
for i=1:size(sig,1)
    
    % sp(i,:) = filtsig(i,filtsig(i,:)>th(i));
    sp{i} = find(filtsig(i,:) > th(i));

end

%% Run computation using BCI_win time windows
% Count total number of spikes per channel and average spectral power in
% each time bin

spcount = zeros(size(sig,1),length(motiontime));
av_power = zeros(size(sig,1),length(motiontime));
% Calculate the number of samples in a time bin at the original Fs
winlength = BCI_win * Fs;

for i=2:length(motiontime)-1
    
    
    %tp = zeros(size(sig,1),1);
    % Get the time window in samples
    id = [motiontime(i-1)*Fs+1, motiontime(i)*Fs];
    
    for ch=1:size(sp,1)
        % Count spikes in time bins
        spcount(ch,i) = length(find((sp{ch,:}*(1/Fs)<motiontime(i)) & ...
            (sp{ch,:}*(1/Fs)>motiontime(i-1))));
        
        av_power(ch,i) = rms(filtsig(ch,id(1):id(2))).^2; 
        
    end
    
    % spcount(i) = sum(tp);
    
   % av_power = 
end

% Try to express av_power in dB
av_power = 10.*log10(av_power);


% % Calculate the spectra
% for ch=1:size(sp,1)
%     
%         % Calculate spectral power in time bins
%         % spectra(:,ch) = pwelch(filtsig(ch,:)',winlength,0,[],Fs,'power',...
%         %     'onesided');
%         
%         av_spectra(ch,1) = mean(pwelch(filtsig(ch,:)',winlength,0,[],Fs,'power',...
%             'onesided');
%         
% end


%% Concatenate the feature vector
features = [spcount;av_power];
features = features(:,2:end);

% At each time stamp the feature vector is ampped onto a joint angle vector
% and errors between current position and target position are calculated

% Create a noisy version of the motion vector
j1_n = j1 + 5.*randn(size(j1));
j2_n = j2 + 5.*randn(size(j2));
j3_n = j3 + 5.*randn(size(j3));
% j1_n = (j1_n ./ norm(j1_n));

% Multi-dimensional vectors
j_t = [j1; j2; j3];

j_n = [j1_n; j2_n; j3_n];

%  For now include all the features
% % Calculate the error vector
% c_t = (j_t-j_n) ./ norm(j_t-j_n);
% 
% 
% %R = ;
% 
% % Find the features that best correlate with the error
% for tst=1:length(motiontime)
%     
%    R{tst} = corrcoef(features,c_t);
%     
% end


%% Select training data. Use first 200 data points in the signals

training_length = 200;
numFeatures = 192;

num_observations = 3;

% Using the Kalman Filter Model
X_t = j_n(:,1:training_length);
X_1 = X_t(:,1:end-1);
X_2 = X_t(:,2:end);
Z_t = features(1:numFeatures,1:training_length);

% Calculate the 4 necessary parameters for the Kalman
A = (X_2*X_1')*inv(X_1*X_1');
H = (Z_t*X_t')*inv(X_t*X_t');
W = ((X_2-A*X_1)*(X_2-A*X_1)')./(training_length-1);
Q = (Z_t-H*X_t)*(Z_t-H*X_t)' ./ (training_length);

%% Let's run the prediction now
x_apr = zeros(size(j_n));
x_apo = zeros(size(j_n));
Pk_apr = zeros(num_observations,num_observations);
Pk_apo = zeros(num_observations,num_observations);
KK = zeros(size(j_n));

z = features(1:numFeatures,:);

for step=training_length:length(motiontime)-1
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
    Pk_apo = (ones(num_observations,num_observations)-KK*H)*Pk_apr;


end

