%%%%%% TEST DECODING TECHNIQUES WITH MOCK-UP TRAINING DATA

clearvars; clc;

addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('Aux Funcs'));
% rmpath(['eeglab',filesep,'plugins']);
%% Open NS file containing Blackrock data files
% Stores continuous data
[file,path] = uigetfile('*.*','Select raw data file');
filename = fullfile(path,file);
%openNSx('report','read',filename,'sample', 'p:short', 's:75');
% openNSx('report','read',filename,'sample', 'p:short');

% Open corresponding NEV file to load the events (extracted spikes)
NEV = openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');
comments = NEV.Data.Comments;
events = NEV.Data.Spikes.TimeStamp;
electrode = NEV.Data.Spikes.Electrode;
Unit = NEV.Data.Spikes.Unit;

% Open NES file
openNSx('report','read','uV',filename,'t:0:40','min','p:short');
% openNSx('report','read','uV',filename,'p:short','s:2');

% openNSx('report','read','uV',filename,'p:short');

% Depending on the data format the structure names are different
forma = filename(end);
%rawData = NS5.Data;
tp = eval(['NS',forma]);
rawData = tp.Data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%     IMPORTANT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The .NSx raw time is computed assuming that the NSP had been reset at the
% beginning of the file! This is always the case when using CENTRAL to run
% the recordings. In fact, every time that a new recording is saved, the
% NSP gets reset. With this assumption, we can match the NSP time stamps
% with data points stored in the raw files.
% Although this might not always be the case, it should work as long as the
% training data files are matched with the correct .NSx raw files. In other
% words the NSP time values might not be "unique" in the ttraining files,
% but they are if matched with the correct raw files.
rawtime = (0:length(rawData)-1)*(1/tp.MetaTags.SamplingFreq); % in seconds


Fs = double(tp.MetaTags.SamplingFreq);


% Drop sampling rate
rawtime = downsample(rawtime,5);
rawData = downsample(rawData',5);
Fsdown = Fs/5;

clear tp;

%% Sham training motion
movement = cos(2*pi*rawtime(1:size(rawData,1))*0.2)';
i1 = find(movement>0.5);
i2 = find(movement<-0.5);
i0 = find(movement>-.5 & movement<0.5);
movement(i0) = 0;
movement(i1) = 1;
movement(i2) = 2;

% Build a sham-data vector to see if the classifiers can be trained
chs2 = [1:3,10:16];
shamData = ones(size(rawData)) .* cos(2*pi*rawtime(1:size(rawData,1))*0.2)'+...
    randn(size(rawData));
shamData(i0) = shamData(i0) + randn(size(shamData(i0)));
shamData(i1,chs2) = shamData(i1,chs2) + 3.*sin(2*pi*rawtime(i1)*22)' + ...
    cos(2*pi*rawtime(i1)*80)'+ randn(size(shamData(i1,chs2)));
shamData(i2,chs2) = shamData(i2,chs2) + 4.*sin(2*pi*rawtime(i2)*93)' + randn(size(shamData(i2,chs2)));

%% Break down the data into 1second time windows with 500ms overlap.
win = 1 * Fsdown;
thisstep = 0.5 * Fsdown;
co = 1;
epochsraw = zeros(win,size(rawData,2),100);
epochssham = zeros(win,size(shamData,2),100);

epoMov = zeros(100,1);

for step=1:thisstep:size(rawData,1)/4
    
epochssham(:,:,co) = shamData((step):(step+win-1),:);    
epochsraw(:,:,co) = rawData((step):(step+win-1),:);
epoMov(co) = movement(step+win-1);
co = co + 1;

end



fr = 1:10:100;
chs = [1:16,18:20,23:30];

%% Use movement as motion labels!!!
% Extract corresponding features
trainingrawData = zeros(size(epochsraw,3),length(fr)*length(chs));
trainingshamData = zeros(size(epochssham,3),length(fr)*length(chs));
% trainingshamData = zeros(size(epochssham,3),length(fr)*1);
trainingLabels = zeros(size(epochsraw,3),1);

for trial=1:size(epochsraw,3)
    tp = epochsraw(:,chs,trial);
    feat = periodogram(tp,[],fr,Fsdown);
    
    trainingrawData(trial,:) = feat(:);
    trainingLabels(trial) = epoMov(trial);
    
    tp = epochssham(:,chs,trial);
    feat = periodogram(tp,[],fr,Fsdown);
    trainingshamData(trial,:) = feat(:);

end

%%%%% Train models %%%%
[Mod1,Acc1] = Classifiers([trainingrawData,trainingLabels],5);

[Mod2,Acc2] = Classifiers([trainingshamData,trainingLabels],5);

%% Predict
decoderoutput2 = ones(size(epoMov)) .* 120;
decoderoutput1 = ones(size(epoMov)) .* 120;

for tt=1:length(epoMov)
        
    decoderoutput1(tt) = predict(Mod1{7},trainingrawData(tt,:));

    decoderoutput2(tt) = predict(Mod2{7},trainingshamData(tt,:));

end

plot(epoMov,'-*')
hold on
plot(decoderoutput1,'--o')
plot(decoderoutput2,'--.')
errId1 = find(epoMov-decoderoutput1 ~= 0);
errId2 = find(epoMov-decoderoutput2 ~= 0);
legend('ideal labels','Raw Data','Sham Data');


%% Try with a Kalman filter

