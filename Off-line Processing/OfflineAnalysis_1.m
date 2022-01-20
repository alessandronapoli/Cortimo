% OFFLINE ANALYSIS


% This script runs all the necessary functions to analyze the raw
% Blackrock files, synchronize them with the corresponding Cortimo training
% data files, extract a set of features use them to train a set of
% classification models and find the optimal feature set / classifier set
% combination set. Each of the intermediate steps will generate dedicated
% files for further analysis. Exisiting files will be overwritten.

% 1) Select the Raw Blackrock data files for a specific session. Select the
% corresponding Cortimo training binary files for the specific session.

% 2) Grab all the epochs and generate to a set of features files using the 
% different paremeters
% specified in ParameterSweep_1.m This will place all the output files in
% the Features folder

% 3) Grab all the features .mat files in the Features folder and identify
% the optimal features with the optimal classification algorithms. The
% final results will be stored in the Results folder.

% The classifiers trained and validated are those in the file Classifiers.m
% Different ones can be introduced here

clearvars; close all; clc;

% STEP 1
% This will grab the raw Blackrock files and Cortimo training files,
% synchronize them and generate new files that are the epoched version of
% the continuous raw data ready to be used for feature extraction
addpath(genpath('C:\Users\aless\OneDrive - Thomas Jefferson University and its Affiliates\Jefferson\Cortimo Project\MatlabCode\Off-line Analysis'));

% Set the analysis parameters


DiscCont = "Discrete"; % Set to Discrete for discrete classifiers

% Parameters for creating new training labels from position data 
% "Target Position", "VR Arm Position", "MyoPro Position"
% ContinuousLabelVR =  "VR Arm Position";
ContinuousLabelVR =  "Target Position";

% Set either to "Ball Location", "Target Location", "Ball Direction"
ContinuousLabel2D = "Ball Location";

% Possible parameters are arrays with the following structure
% [-1,0,1,2] where:
% -1 is Rest
%  0 is Hold
%  1 is Flex
%  2 is Extend

% If any of the above parameters are not used, replace the value in the array
% with []

%DiscLabelsHand =  [-1,0,1,2];
%DiscLabelsElbow = [-1,0,1,2];
DiscLabelsHand = [];
DiscLabelsElbow = [];

[progress,isCont] = RawDataToEpochs_4(DiscCont,ContinuousLabelVR,ContinuousLabel2D,DiscLabelsHand,...
    DiscLabelsElbow);

 delete(progress.msgBox);
 
%%
% STEP 2
% This might take some time, depending on the number of parameters that are
% included in the parameter sweep phase. This step includes the feature
% extraction step.
% ParameterSweep_1;

%%%%%% LFP  %%%%%%
chsLFP = [1:16,18:20,23:30];
% Parameter 1
TimeWins = [0.3,0.5,0.8,1,2];
% Parameter 2
overlaps = [0,0.25,0.5,0.75];
% Parameter 3
freqBins = [3,5,10,20,50];
% Parameter 4
% freqRanges = {[0,800],[0,500],[0,250],[0,100],...
%      [2,800],[2,500],[2,250],[2,100]};
freqRanges = {[0,250],[0,100],[2,500],[2,250],[2,100]};


%%%%%% SPIKES   %%%%%%
% Spike parameters
spkChs = 1:128;
spkTWs = [0.05,0.1];

% This performs paramter sweep for all features
%ParameterSweepAllFeats_1(chsLFP,TimeWins,overlaps,freqBins,freqRanges,...
%    spkChs,spkTWs,isCont);

OnlyUseTheseLabels = [0,1,2];

% This performs parameter sweep for LFPs and Spikes separately
ParameterSweep_2(chsLFP,TimeWins,overlaps,freqBins,freqRanges,spkChs,...
    spkTWs,isCont,OnlyUseTheseLabels);


%%
% STEP 3
% This might take some time depending on the time necessary to train the
% classification models. Reducing the number of features will speed up this
% step as well.
thisClassifier = 1;

Find_Optimal_Classifiers_1(isCont,thisClassifier);


%% Look for the best features
% FindBestFeats([]);


%% Look for the optimal classifiers using the best features
% Find_Optimal_Classifiers_reducedFeats();

