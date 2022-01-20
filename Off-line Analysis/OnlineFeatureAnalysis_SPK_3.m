function [tVec,ff,ElbData,ElbZdata,HanData,HanZdata] = OnlineFeatureAnalysis_SPK_3(ff)

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

[fileList,pathT] = uigetfile('*.bin','Select training file with Spike data','MultiSelect',...
    'on');


ElbData = [];
ElbZdata = [];
HanData = [];
HanZdata = [];
tVec = [];
tVecElb = [];
tVecHan = [];


%%%%%%%%%%%%%%%%%%%%%%%
if fileList==0
  % user pressed cancel
  return;
end
%%%%%%%%%%%%%%%%%%%%%%


% dataElbow = [];
% dataHand = [];
% dataElbow2 = [];
% dataHand2 = [];
% dataElbowZ = [];
% dataHandZ = [];
% dataElbowZ2 = [];
% dataHandZ2 = [];

ElbowS = cell(1,3);
ElbowZscores = cell(size(ElbowS));

HandS = cell(1,3);
HandZscores = cell(size(HandS));

% ElbowS = cell(1,3);
% HandS = cell(1,3);
% ElbowZ = cell(1,3);
% HandZ = cell(1,3);


ElbowavS = [];
HandavS = [];
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
    SPKfeats = [];
    
for file=1:LEN
    
    fileT = fileList{file};
    
filenameFeat = fullfile(pathT,fileT);

% Open both files, the training data file and the EEG feature file
filenameTrain = strrep(filenameFeat,"spikeFeatures","ArmVR");
try
    [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile(filenameTrain)];
catch
    [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile_noDecoderOutput(filenameTrain)];

end

try
    [SPKfeats] = [SPKfeats; loadSPKFeatureFile(filenameFeat)];
catch
    [SPKfeats] = [SPKfeats; loadSPKFeatureFile(filenameFeat)];

end

end


if(isempty(SPKfeats) || isempty(VRTrainingLabels))
    return;
end
%%

ff.usersettings.spktimewin = SPKfeats.AnalysisTimewinDur(1);
cco = SPKfeats(1,5:size(SPKfeats,2)).Properties.VariableNames;
cco = extractBetween(cco,'_','_');
ff.usersettings.ch_to_spi = strjoin(cco);

Training_TS = VRTrainingLabels.NSPTime;
Feature_TS = SPKfeats.NSPts;

TrainMatCycle = VRTrainingLabels.MatlabCycle;

spks = table2array(SPKfeats(:,5:size(SPKfeats,2)));

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



% % Select 5 seconds worth of data right at the start of the training
% % session
stID = Fids(1);
numofwins = 50; % Take 5 seconds worth of baseline data
% 
muSPK = mean(spks(stID:numofwins-1,:));
stSPK = std(spks(stID:numofwins-1,:));

% BLSpectrum = SPKfeats(stID:stID+numofwins-1,5);
% % muBLSp = cellfun(@mean,BLSpectrum);
% % stBLSp = cellfun(@std,BLSpectrum);
% BLSpec = zeros(size(BLSpectrum{1,1},1),size(BLSpectrum{1,1},2),length(BLSpectrum));
% for i=1:length(BLSpectrum)
% 
%     BLSpec(:,:,i) = BLSpectrum{i}; 
%     
% end
% 
% % These are the parameters that can be used to create the Z scores
% muBLSp = mean(BLSpec,3);
% stBLSp = std(BLSpec,[],3);


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
    ElbowS{cond} = spks(FeatId,:);
    ElbowZscores{cond} = (spks(FeatId) - muSPK) ./ stSPK;
        
    tVecElb = [tVecElb;el];
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
    HandS{cond} = spks(FeatId,:);
    HandZscores{cond} = (spks(FeatId) - muSPK) ./ stSPK;
        
       % Spec = cell2mat(EEGfeats(FeatId,5));
       % ElbowSpec{cond} = [ElbowSpec{cond}; Spec];
       
    tVecHan = [tVecHan;el];

   end
    
end


%% Average the matrices now

for cond=1:length(ElbowS)
    
    ElbowavS{cond} = mean(ElbowS{cond});
    ElbowavZ{cond} = mean(ElbowZscores{cond});
    HandavS{cond} = mean(HandS{cond});
    HandavZ{cond} = mean(HandZscores{cond});
    
end


%% Feed the data to the classification learner

% clasLab = [0,1,2];
ElSp = cell(size(ElbowS));
HaSp = cell(size(HandS));
ElZ = cell(size(ElbowS));
HaZ = cell(size(HandS));

%%%
clasLab = [0,1,2];

for cond=1:length(ElbowS)
            if(~isempty(ElbowS{cond}))
    
     ElSp{cond} = [ElbowS{cond}, clasLab(cond).*ones(size(ElbowS{cond},1),1)];
     ElZ{cond} = [ElbowZscores{cond}, clasLab(cond).*ones(size(ElbowZscores{cond},1),1)];
     
     HaSp{cond} = [HandS{cond}, clasLab(cond).*ones(size(HandS{cond},1),1)];
     HaZ{cond} = [HandZscores{cond}, clasLab(cond).*ones(size(HandZscores{cond},1),1) ];
     
            end
end

% Concatenate matrices
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

tVec{1} = tVecElb;
tVec{2} = tVecHan;


end
%% Try extracting the features that have the highest discriminant power




