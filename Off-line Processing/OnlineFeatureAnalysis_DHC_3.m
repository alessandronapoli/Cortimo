function [ff,ElbData,ElbZdata,HanData,HanZdata] = OnlineFeatureAnalysis_DHC_3(ff)

% Offline analysis for Gamma breach project
% This analysis focuses on analyzing EEG features extracted from NSP data.
% The features are synchronized with training data derived from the Cortimo
% suite.

% Features used in this analysis are power spectra calculated in 
% 3Hz-frequency bands from 0Hz to 1000Hz and variable number of channels based 
% on the user selection.
% Sampling frequency is usually 2000 sps and analysis time intervals for the
% power calculation is 500ms. These parameters can vary and are all stored
% in the feature files.
% clearvars; close all; clc;
dd = pwd;
for i=1:2
idcs   = strfind(dd,filesep);
dd = dd(1:idcs(end)-1);
end

addpath([dd,filesep,'Off-line Analysis',filesep,'Aux Funcs']);

[fileList,pathT] = uigetfile('*.bin','Select file with EEG data','MultiSelect',...
    'on');

% dataElbow = [];
% dataHand = [];
% dataElbow2 = [];
% dataHand2 = [];
% dataElbowZ = [];
% dataHandZ = [];
% dataElbowZ2 = [];
% dataHandZ2 = [];

ElbowSpec = cell(1,3);
ElbowZscores = cell(size(ElbowSpec));

HandSpec = cell(1,3);
HandZscores = cell(size(HandSpec));

% ElbowS = cell(1,3);
% HandS = cell(1,3);
% ElbowZ = cell(1,3);
% HandZ = cell(1,3);


ElbowavSpec = [];
HandavSpec = [];
ElbowavZ = [];
HandavZ = [];


if(iscell(fileList))
    LEN = length(fileList);
    
else
    LEN = 1;
    tp{1} = fileList;
    fileList = tp;
end

    VRTrainingLabels = [];
    EEGfeats = [];
    
for file=1:LEN
    
    fileT = fileList{file};
    
filenameFeat = fullfile(pathT,fileT);

% Open both files, the training data file and the EEG feature file
filenameTrain = strrep(filenameFeat,"FreqFeatures","ArmVR");
try
    [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile(filenameTrain)];
catch
    [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile_noDecoderOutput(filenameTrain)];

end

try
    [EEGfeats] = [EEGfeats; loadEEGFeatureFile(filenameFeat)];
catch
    [EEGfeats] = [EEGfeats; loadEEGFeatureFile_old(filenameFeat)];

end

%%
tp = EEGfeats{1,4};

if(size(EEGfeats{1},2)>4)
        % NEW file extraction

    fr = str2num(tp(5:end));
Fs = EEGfeats{1,2}(5);
Feature_TS = zeros(size(EEGfeats,1),1);
Training_TS = VRTrainingLabels.NSPTime;
% FeatMatCycle = zeros(size(EEGfeats,1),1);
Ch = str2num(EEGfeats{1,3}(13:end));
% ChLabels = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
%     'F7','F8','T3','T4','T5','T6','Z-Gr','Fz','Cz','Pz','A1','A2','C2',...
%     'C6','FC4','CP4','CP2','CP6','FC2','FC6'};

ChLabels = num2cell(Ch);

    TimeUpdate = EEGfeats{1,2}(4);
    FFTtimeduration = EEGfeats{1,2}(3);
    
else
    
    % OLD extraction methods
    fr = str2num(tp(5:end));
Fs = EEGfeats{1,2}(4);
Feature_TS = zeros(size(EEGfeats,1),1);
Training_TS = VRTrainingLabels.NSPTime;
% FeatMatCycle = zeros(size(EEGfeats,1),1);
Ch = str2num(EEGfeats{1,3}(13:end));
ChLabels = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T3','T4','T5','T6','Z-Gr','Fz','Cz','Pz','A1','A2','C2',...
    'C6','FC4','CP4','CP2','CP6','FC2','FC6'};
     TimeUpdate = EEGfeats{1,2}(3);
    FFTtimeduration = 1;
    
    
end

Frange = [fr(1),fr(length(fr))+1];
Fstep = diff(fr(1:2));


ff.usersettings.Frange = Frange;
ff.usersettings.Fstep = Fstep;
ff.usersettings.EEGFFTWin = FFTtimeduration;
ff.usersettings.EEGtimewin = TimeUpdate;
ff.usersettings.ch_to_freq = num2str(Ch);



for i=1:size(EEGfeats,1)
Feature_TS(i) = EEGfeats{i,2}(2);       % in seconds
% FeatMatCycle(i) = EEGfeats{i,2}(1);
end

TrainMatCycle = VRTrainingLabels.MatlabCycle;
%plot(fr,10*log10(EEGfeats{100,5}(:,1)))

% Turn the data table to a matrix


% Grab and average results
% Rest
% idsR = find(VRTrainingLabels{:,4} <0);


%% Synchronize the data streams
% tt = intersect(VRTrainingLabels{:,3},TS);
% firstelement = tt(1);
% lastelement = TS(end);
% 
% idfirst = find(VRTrainingLabels{:,3} == firstelement,1);
% idlast = find(VRTrainingLabels{:,3} == lastelement,1);

% len = (idlast-idfirst);

% timeS = unique(VRTrainingLabels{:,3});

% New time vector
% VRtime = firstelement:1/50.025:lastelement; % This needs to be verified

% ----------------------------------------------------------------------- % 
% This is not a very precise way to sync the data streams and pull data
% epochs out of the continuous data. Given the slow (1 second of data 
% refreshed every 500 ms in this case) data rate of the EEG features, being
% not too precise in frame identification shouldn't affect performance too
% much. Other algorithms to have a better resolution in timing data has also 
% been developed for the Cortimo App. See SynchTrainingLabels.m for that purpose.

b = find(VRTrainingLabels.Trial == 0);

% [el,Fids,Tids] = intersect(FeatMatCycle,TrainMatCycle);
[el,Fids,Tids] = intersect(round(Feature_TS,6),round(Training_TS,6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* In case of multiple recordings/training files,the baselines should
%* correspond to each data session that we would like to process. In other words,
%* the Z scores for each data session should be calculated based on
% baseline values specific for that session



% Select 5 seconds worth of data right at the start of the training
% session
stID = Fids(1);
numofwins = round(5 / (TimeUpdate/1000)); % Take 5 seconds worth of baseline data

BLSpectrum = EEGfeats(stID:stID+numofwins-1,5);
% muBLSp = cellfun(@mean,BLSpectrum);
% stBLSp = cellfun(@std,BLSpectrum);
BLSpec = zeros(size(BLSpectrum{1,1},1),size(BLSpectrum{1,1},2),length(BLSpectrum));
for i=1:length(BLSpectrum)

    BLSpec(:,:,i) = BLSpectrum{i}; 
    
end

% These are the parameters that can be used to create the Z scores
muBLSp = mean(BLSpec,3);
stBLSp = std(BLSpec,[],3);


% Let's use the labels to group the trials
% Elbow first
idsElbow = cell(1,3);
idsElbow{1} = find(VRTrainingLabels{:,8} == 0);
idsElbow{2} = find(VRTrainingLabels{:,8} == 1);
idsElbow{3} = find(VRTrainingLabels{:,8} == 2);

% Hand
idsHand = cell(1,3);
idsHand{1} = find(VRTrainingLabels{:,7} == 0);
idsHand{2} = find(VRTrainingLabels{:,7} == 1);
idsHand{3} = find(VRTrainingLabels{:,7} == 2);

%% ELBOW LABELS
% Find the corresponding feature vectors
% cycles = unique(VRTrainingLabels{:,1}(idsElbow0,1));

% ElbowSpec = cell(1,length(idsElbow));
% ElbowZscores = cell(size(ElbowSpec));
newTraining_TS = round(Training_TS,6);
newFeature_TS = round(Feature_TS,6);

for cond=1:length(idsElbow)
   
    idsE = idsElbow{cond};
    
    if(~isempty(idsE))
    [el,TrainId,FeatId] = intersect(newTraining_TS(idsE),newFeature_TS);
    ElbowSpec{cond} = zeros([size(EEGfeats{FeatId(1),5}),length(FeatId)]);
    ElbowZscores{cond} = zeros([size(EEGfeats{FeatId(1),5}),length(FeatId)]);
       % Spec = cell2mat(EEGfeats(FeatId,5));
       % ElbowSpec{cond} = [ElbowSpec{cond}; Spec];
        
        for trial=1:length(FeatId)
        ElbowSpec{cond}(:,:,trial) = EEGfeats{FeatId(trial),5};    
        ElbowZscores{cond}(:,:,trial) = (EEGfeats{FeatId(trial),5} - muBLSp) ./ stBLSp;
        % ElbowZscores{cond} = [ElbowZscores{cond}; Spec2];
        end
        
    end
    
end


%% HAND LABELS
% Find the corresponding feature vectors
% cycles = unique(VRTrainingLabels{:,1}(idsElbow0,1));

% HandSpec = cell(1,length(idsHand));
% HandZscores = cell(size(HandSpec));

for cond=1:length(idsHand)
   
   idsH = idsHand{cond};
    
   if(~isempty(idsH))
    [el,TrainId,FeatId] = intersect(newTraining_TS(idsH),newFeature_TS);
    HandSpec{cond} = zeros([size(EEGfeats{FeatId(1),5}),length(FeatId)]);
    HandZscores{cond} = zeros([size(EEGfeats{FeatId(1),5}),length(FeatId)]);
       % Spec = cell2mat(EEGfeats(FeatId,5));
       % ElbowSpec{cond} = [ElbowSpec{cond}; Spec];
        
        for trial=1:length(FeatId)
        HandSpec{cond}(:,:,trial) = EEGfeats{FeatId(trial),5};    
        HandZscores{cond}(:,:,trial) = (EEGfeats{FeatId(trial),5} - muBLSp) ./ stBLSp;
        % ElbowZscores{cond} = [ElbowZscores{cond}; Spec2];
        end
        
   end
    
end


%% Average the matrices now

for cond=1:length(ElbowSpec)
    
    ElbowavSpec{cond} = mean(ElbowSpec{cond},3);
    ElbowavZ{cond} = mean(ElbowZscores{cond},3);
    HandavSpec{cond} = mean(HandSpec{cond},3);
    HandavZ{cond} = mean(HandZscores{cond},3);
    
end




end

%% Get the R-squared values out here before averaging the matrices
RsqMatElbS = cell(1,3);
RsqMatHanS = cell(1,3);

RsqMatElbZ = cell(1,3);
RsqMatHanZ = cell(1,3);

% Elbow
% Rest VS. Flex
if(~isempty(ElbowSpec{2}))
RsqMatElbS{1} = calc_rsqu(ElbowSpec{1},ElbowSpec{2},1);
RsqMatElbZ{1} = calc_rsqu(ElbowZscores{1},ElbowZscores{2},1);
end

% Rest VS. Extend
if(~isempty(ElbowSpec{3}))
RsqMatElbS{2} = calc_rsqu(ElbowSpec{1},ElbowSpec{3},1);
RsqMatElbZ{2} = calc_rsqu(ElbowZscores{1},ElbowZscores{3},1);
end

if(~isempty(ElbowSpec{2}) && ~isempty(ElbowSpec{3}))
% Flex VS. Extend
RsqMatElbS{3} = calc_rsqu(ElbowSpec{2},ElbowSpec{3},1);
RsqMatElbZ{3} = calc_rsqu(ElbowZscores{2},ElbowZscores{3},1);
end

% Hand
% Rest VS. Flex
if(~isempty(HandSpec{2}))
RsqMatHanS{1} = calc_rsqu(HandSpec{1},HandSpec{2},1);
RsqMatHanZ{1} = calc_rsqu(HandZscores{1},HandZscores{2},1);
end

% Rest VS. Extend
if(~isempty(HandSpec{3}))
RsqMatHanS{2} = calc_rsqu(HandSpec{1},HandSpec{3},1);
RsqMatHanZ{2} = calc_rsqu(HandZscores{1},HandZscores{3},1);
end

% Flex VS. Extend
if(~isempty(HandSpec{2}) && ~isempty(HandSpec{3}))
RsqMatHanS{3} = calc_rsqu(HandSpec{2},HandSpec{3},1);
RsqMatHanZ{3} = calc_rsqu(HandZscores{2},HandZscores{3},1);
end

labelsR = {'Rest VS. Flex','Rest VS. Extend','Flex VS. Extend'};
figure(1)
for cond=1:length(RsqMatElbS)
    if(~isempty(RsqMatElbS{cond}))
subplot(3,1,cond)
% surf(1:size(RsqMatElbS{cond},2),fr,10*log10(RsqMatElbS{cond}))
surf(1:size(RsqMatElbS{cond},2),fr,RsqMatElbS{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(['ELBOW  ',labelsR{cond}]);
colormap(jet);
% caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Spectrum R2';
    end
end

figure(2)
for cond=1:length(RsqMatHanS)
        if(~isempty(RsqMatHanS{cond}))
subplot(3,1,cond)
% surf(1:size(RsqMatHanS{cond},2),fr,10*log10(RsqMatHanS{cond}))
surf(1:size(RsqMatHanS{cond},2),fr,RsqMatHanS{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(['HAND  ',labelsR{cond}]);
colormap(jet);
% caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Spectrum R2';
        end
end


figure(3)
for cond=1:length(RsqMatElbZ)
    if(~isempty(RsqMatElbZ{cond}))
subplot(3,1,cond)
% surf(1:size(RsqMatElbS{cond},2),fr,10*log10(RsqMatElbS{cond}))
surf(1:size(RsqMatElbZ{cond},2),fr,RsqMatElbZ{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(['ELBOW  ',labelsR{cond}]);
colormap(jet);
% caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Z-Score R2';
    end
end

figure(4)
for cond=1:length(RsqMatHanZ)
    if(~isempty(RsqMatHanZ{cond}))
subplot(3,1,cond)
% surf(1:size(RsqMatHanS{cond},2),fr,10*log10(RsqMatHanS{cond}))
surf(1:size(RsqMatHanZ{cond},2),fr,RsqMatHanZ{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(['HAND  ',labelsR{cond}]);
colormap(jet);
% caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Z-Score R2';
    end
end



%% Plot the average spectra

% figure
% subplot(2,1,1)
% plot(fr,10*log10(ElbowavSpectrum0))
% hold on
% plot(fr,10*log10(ElbowavSpectrum1))
% plot(fr,10*log10(ElbowavSpectrum2))
% legend('Rest', 'Flex', 'Extend')
% title('Elbow Motion')
% subplot(2,1,2)
% plot(fr,10*log10(HandavSpectrum0))
% hold on
% plot(fr,10*log10(HandavSpectrum1))
% plot(fr,10*log10(HandavSpectrum2))
% legend('Rest', 'Flex', 'Extend')
% title('Hand Motion')
% 
labels1 = {'Elbow Rest','Elbow Flex','Elbow Extend'};
labels2 = {'Hand Rest','Hand Flex','Hand Extend'};

figure(5)
for cond=1:length(ElbowavSpec)
    if(~isempty(ElbowavSpec{cond}))
subplot(3,1,cond)
surf(1:size(ElbowavSpec{cond},2),fr,10*log10(ElbowavSpec{cond}))
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(labels1{cond});
colormap(jet);
caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Mean Spectrum [dB]';
    end
end


figure(6)
for cond=1:length(HandavSpec)
    if(~isempty(HandavSpec{cond}))
subplot(3,1,cond)
surf(1:size(HandavSpec{cond},2),fr,10*log10(HandavSpec{cond}))
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(labels2{cond});
colormap(jet);
caxis([-20 30]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Mean Spectrum [dB]';
    end
end


figure(7)
for cond=1:length(ElbowavZ)
    if(~isempty(ElbowavZ{cond}))
subplot(3,1,cond)
surf(1:size(ElbowavZ{cond},2),fr,ElbowavZ{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(labels1{cond});
colormap(jet);
caxis([-3 3]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Mean Z Score Spectrum';
    end
end


figure(8)
for cond=1:length(HandavZ)
    if(~isempty(HandavZ{cond}))
subplot(3,1,cond)
surf(1:size(HandavZ{cond},2),fr,HandavZ{cond})
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
title(labels2{cond});
colormap(jet);
caxis([-3 3]); % Change limits
axis('tight')
c = colorbar;
c.Visible = 'on';
c.Label.String = 'Mean Z Score Spectrum';
    end
end


%% Feed the data to the classification learner

% clasLab = [0,1,2];
ElSp = cell(size(ElbowSpec));
HaSp = cell(size(HandSpec));
ElZ = cell(size(ElbowSpec));
HaZ = cell(size(HandSpec));


% dataElbow = [];
% dataHand = [];
% dataElbow2 = [];
% dataHand2 = [];
% dataElbowZ = [];
% dataHandZ = [];
% dataElbowZ2 = [];
% dataHandZ2 = [];


% Reshape the matrices for machine learning training
% Concatenate the frequency bins and the channels

% C = permute(ElbowS{1},[1 2 3]);
% ChSelect = [5,6,23,26];
% FreqSelect = [1:200];
% ChSelect = [1:16,18:20,23:30]; % Exclude 17 (ground), 21(A1), 22(A2)
% FreqSelect = [1:100];
ChSelect = Ch;
FreqSelect = find(fr == fr);
clasLab = [0,1,2];

for cond=1:length(ElbowSpec)
            if(~isempty(ElbowSpec{cond}))

    E1 = zeros(size(ElbowSpec{cond},3),length(ChSelect)*length(FreqSelect));
    E2 = zeros(size(ElbowZscores{cond},3),length(ChSelect)*length(FreqSelect));

    for trial=1:size(ElbowSpec{cond},3)
        tp = ElbowSpec{cond}(FreqSelect,ChSelect,trial);
        E1(trial,:) = tp(:);
        
        tp = ElbowZscores{cond}(FreqSelect,ChSelect,trial);
        E2(trial,:) = tp(:);
  
    end
            else
                E1 = [];
                E2 = [];
            end
    
    if(~isempty(HandSpec{cond}))

    H1 = zeros(size(HandSpec{cond},3),length(ChSelect)*length(FreqSelect));
    H2 = zeros(size(HandZscores{cond},3),length(ChSelect)*length(FreqSelect));

    for trial=1:size(HandSpec{cond},3)
        tp = HandSpec{cond}(FreqSelect,ChSelect,trial);
        H1(trial,:) = tp(:);
        tp = HandZscores{cond}(FreqSelect,ChSelect,trial);
        H1(trial,:) = tp(:);
        
    end
    
    else
        H1 = [];
        H2 = [];
    end
    
     ElSp{cond} = [E1, clasLab(cond).*ones(size(E1,1),1)];
     ElZ{cond} = [E2, clasLab(cond).*ones(size(E2,1),1)];
     
     HaSp{cond} = [H1, clasLab(cond).*ones(size(H1,1),1)];
     HaZ{cond} = [H2, clasLab(cond).*ones(size(H2,1),1) ];
     
     
     
%     % lll = size(ElbowS{cond},3);
%     E1 = reshape(ElbowSpec{cond}(FreqSelect,ChSelect,:),[],size(ElbowSpec{cond}(FreqSelect,ChSelect,:),3),1)';
%     E2 = reshape(ElbowZ{cond}(FreqSelect,ChSelect,:),[],size(ElbowZ{cond}(FreqSelect,ChSelect,:),3),1)';
%     H1 = reshape(HandS{cond}(FreqSelect,ChSelect,:),[],size(HandS{cond}(FreqSelect,ChSelect,:),3),1)';
%     H2 = reshape(HandZ{cond}(FreqSelect,ChSelect,:),[],size(HandZ{cond}(FreqSelect,ChSelect,:),3),1)';
%     
%     ElSp{cond} = [E1, clasLab(cond).*ones(size(E1,1),1)];
%     ElZ{cond} = [E2, clasLab(cond).*ones(size(E2,1),1)];
%     
%     HaSp{cond} = [H1, clasLab(cond).*ones(size(H1,1),1)];
%     HaZ{cond} = [H2, clasLab(cond).*ones(size(H2,1),1) ];
%     
% %     ElSp{cond} = [ElbowSp{cond}, clasLab(cond).*ones(size(ElbowSpec{cond},1),1)];
% %     ElZ{cond} = [ElbowZscores{cond}, clasLab(cond).*ones(size(ElbowSpec{cond},1),1)];
% %     
% %     HaSp{cond} = [HandSpec{cond}, clasLab(cond).*ones(size(HandSpec{cond},1),1)];
% %     HaZ{cond} = [HandZscores{cond}, clasLab(cond).*ones(size(HandSpec{cond},1),1) ];
end

% Concatenate matrices from different training files
dataElbow = [ElSp{1}; ElSp{2}; ElSp{3}];
dataElbowZ = [ElZ{1}; ElZ{2}; ElZ{3}];

dataHand = [HaSp{1}; HaSp{2}; HaSp{3}];
dataHandZ = [HaZ{1}; HaZ{2}; HaZ{3}];


%% Special case for MyoPro single mode use
% When in single mode, the MyoPro requires the user to focus on a single
% muslce group activation for activating the motors in one direction, the
% opposite motion is evoked by rest, that is signals below threshold. In
% this context, the decoder is a two-class classification problem.

%%% Specific for Subject 0001:
% Elbow was in extension mode
% Hand was in close mode
modElSp{1} = ElSp{1};
modElSp{1}(:,size(modElSp{1},2)) = 1;
modElZ{1} = ElZ{1};
modElZ{1}(:,size(modElZ{1},2)) = 1;

modHaSp{1} = HaSp{1};
modHaSp{1}(:,size(modHaSp{1},2)) = 2;
modHaZ{1} = HaZ{1};
modHaZ{1}(:,size(modHaZ{1},2)) = 2;


dataElbow2 = [modElSp{1}; ElSp{2}; ElSp{3}];
dataElbowZ2 = [modElZ{1}; ElZ{2}; ElZ{3}];

dataHand2 = [modHaSp{1}; HaSp{2}; HaSp{3}];
dataHandZ2 = [modHaZ{1}; HaZ{2}; HaZ{3}];



ElbData = dataElbow;
ElbZdata = dataElbowZ;
HanData = dataHand;
HanZdata = dataHandZ;



%% Try extracting the features that have the highest discriminant power




