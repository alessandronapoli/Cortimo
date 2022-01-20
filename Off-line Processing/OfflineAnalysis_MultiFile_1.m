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

dirToProcess = uigetdir(pwd,'Select the Folder containing the folders to be processed');
thesefolders = dir(dirToProcess);

thesefolders = thesefolders(3:end);

run([dirToProcess,filesep,'ParameterList_1']);
% Parameters are:
% 1) DiscCont
% 2) ContinuousLabelVR
% 3) ContinuousLabel2D
% 4) DiscLabelsHand
% 5) DiscLabelsElbow
% 6) OnlyTheseLabels

spkChs = 33:96;
spkTWs = 0.05;

ThisClassifier = 1;
isCont = 0;



for folder=1:length(thesefolders)

    tp = thesefolders(folder).name;
    
       dirs = [thesefolders(folder).folder,filesep,tp];
       FL = dir(dirs);
       FL = FL(3:end);
       

       
        for ii=1:length(FL)
       
            if(contains(FL(ii).name,'.nev'))
                
                % Identify corresponding training files. Let's start by
                % loading the VR arm training data sets
                TrainIDs = find(contains({FL.name},'_ArmVR_TRAINING_.bin') == 1);
                
                TrainingFiles = {[FL(TrainIDs).folder,filesep,FL(TrainIDs).name]};
                
                % Set to Discrete for discrete classifiers
                DiscCont = ParameterList{ii,1};
                % Parameters for creating new training labels from position data 
                % "Target Position", "VR Arm Position", "MyoPro Position"
                ContinuousLabelVR = ParameterList{ii,2};
                % Set either to "Ball Location", "Target Location", "Ball Direction"
                ContinuousLabel2D = ParameterList{ii,3};
                % Possible parameters are arrays with the following structure
                % [-1,0,1,2] where:
                % -1 is Rest
                %  0 is Hold
                %  1 is Flex
                %  2 is Extend
                DiscLabelsHand = ParameterList{ii,4};
                DiscLabelsElbow = ParameterList{ii,5};
                % Only use specific labels. These are training-dependent 
                OnlyUseTheseHandLabels = ParameterList{ii,6};
                OnlyUseTheseElbowLabels = ParameterList{ii,7};
                
                % Extract Epochs from data
        [progress,isCont] = RawDataToEpochs_Multiple_4(DiscCont,...
            ContinuousLabelVR,ContinuousLabel2D,DiscLabelsHand,...
            DiscLabelsElbow,[FL(ii).folder,filesep,FL(ii).name],TrainingFiles);

 delete(progress.msgBox);


            end
        end
        
        
end


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
%%% Spike parameters
spkChs = 30:100;
spkTWs = 0.05;
OnlyUseTheseHandLabels = [-1,1,2];
OnlyUseTheseElbowLabels = [-1,1,2];

% This performs parameter sweep for LFPs and Spikes separately
ParameterSweep_2(chsLFP,TimeWins,overlaps,freqBins,freqRanges,spkChs,...
    spkTWs,isCont,OnlyUseTheseHandLabels,OnlyUseTheseElbowLabels);


%% STEP 3
% Raster plots of repeated trials for each movement
% Hand first

filelist = dir('Features');
filelist = filelist(3:end);

MeanRestHan = zeros(length(filelist),length(spkChs));
MeanFlexHan = zeros(length(filelist),length(spkChs));
MeanExtHan = zeros(length(filelist),length(spkChs));

MeanRestElb = zeros(length(filelist),length(spkChs));
MeanFlexElb = zeros(length(filelist),length(spkChs));
MeanExtElb = zeros(length(filelist),length(spkChs));


for fi=1:length(filelist)

    tp = filelist(fi).name;
    
    if(contains(tp,'Han'))
       dirs = [filelist(fi).folder,filesep,filelist(fi).name];
       FL = dir(dirs);
       FL = FL(3:end);
       
        for ii=1:length(FL)
            
        load([FL(ii).folder,filesep,FL(ii).name])
        
        
% % Rest / Hold
% figure(1)
% id = find(feats(:,size(feats,2)) == 0);
% plot(feats(id,1:128))
% 
% % Flex
% figure(2)
% id = find(feats(:,size(feats,2)) == 1);
% plot(feats(id,1:128))
% 
% % Extend
% figure(3)
% id = find(feats(:,size(feats,2)) == 2);
% plot(feats(id,1:128))

% New channel IDs
idC = spkChs - spkChs(1) + 1;

% Average trials
id = find(feats(:,size(feats,2)) == OnlyUseTheseHandLabels(1));
AvRest = mean(feats(id,idC));

id = find(feats(:,size(feats,2)) == OnlyUseTheseHandLabels(2));
AvFlex = mean(feats(id,idC));

id = find(feats(:,size(feats,2)) == OnlyUseTheseHandLabels(3));
AvExtend = mean(feats(id,idC));


MeanRestHan(ii,:) = AvRest;
MeanFlexHan(ii,:) = AvFlex;
MeanExtHan(ii,:) = AvExtend;


% figure
% plot(spkChs,AvRest)
% hold on
% plot(spkChs,AvFlex)
% plot(spkChs,AvExtend)
% title(['Average Spike counts in ',FL(ii).name(5:8), 's time bins']);
% xlabel('Channel Number');
% ylabel('Spike Count in time bins')
% legend('Hand Hold','Hand Flexion','Hand Extension');
% hold off
% 
% figure
% hold on
% grouped = [AvRest;AvFlex;AvExtend]';
% bar([spkChs;spkChs;spkChs]',grouped)
% legend('Hand Hold','Hand Flexion','Hand Extension');
% xlabel('Channel Number');
% ylabel('Spike Count in time bins')
% title(['Average Spike counts in ',FL(ii).name(5:8), 's time bins']);
% hold off

        end
    end
    
    if(contains(tp,'Elb'))
    dirs = [filelist(fi).folder,filesep,filelist(fi).name];
       FL = dir(dirs);
       FL = FL(3:end);
       
        for ii=1:length(FL)
            
        load([FL(ii).folder,filesep,FL(ii).name])
% % Rest / Hold
% figure(1)
% id = find(feats(:,size(feats,2)) == 0);
% plot(feats(id,1:128))
% 
% % Flex
% figure(2)
% id = find(feats(:,size(feats,2)) == 1);
% plot(feats(id,1:128))
% 
% % Extend
% figure(3)
% id = find(feats(:,size(feats,2)) == 2);
% plot(feats(id,1:128))

idC = spkChs - spkChs(1) + 1;

% Average trials
id = find(feats(:,size(feats,2)) == OnlyUseTheseElbowLabels(1));
AvRest = mean(feats(id,idC));

id = find(feats(:,size(feats,2)) == OnlyUseTheseElbowLabels(2));
AvFlex = mean(feats(id,idC));

id = find(feats(:,size(feats,2)) == OnlyUseTheseElbowLabels(3));
AvExtend = mean(feats(id,idC));


MeanRestElb(ii,:) = AvRest;
MeanFlexElb(ii,:) = AvFlex;
MeanExtElb(ii,:) = AvExtend;

% figure
% %plot(AvRest)
% stem(spkChs,AvRest,'LineStyle','--','Marker','*')
% hold on
% %plot(AvFlex)
% stem(spkChs,AvFlex,'LineStyle','--','Marker','square')
% %plot(AvExtend)
% stem(spkChs,AvExtend,'LineStyle','--','Marker','diamond')
% title(['Average Spike counts in ',FL(ii).name(5:8), 's time bins']);
% xlabel('Channel Number');
% ylabel('Spike Count in time bins')
% legend('Elbow Hold','Elbow Flexion','Elbow Extension');
% 
% hold off
% 
% 
% figure
% hold on
% grouped = [AvRest;AvFlex;AvExtend]';
% bar([spkChs;spkChs;spkChs]',grouped)
% legend('Elbow Hold','Elbow Flexion','Elbow Extension');
% xlabel('Channel Number');
% ylabel('Spike Count in time bins')
% title(['Average Spike counts in ',FL(ii).name(5:8), 's time bins']);
% hold off

        end

    end
    
    
end

MeanFeatsHan = [mean(MeanRestHan);mean(MeanFlexHan);mean(MeanExtHan)];
MeanFeatsElb = [mean(MeanRestElb);mean(MeanFlexElb);mean(MeanExtElb)];


dest = [pwd,filesep,'AverageFeatures'];
 if(exist(dest,'dir') ~= 7)
       mkdir(dest);   %create the directory
 end
         
        save([dest,filesep,'MeanFeatHand','.mat'],'MeanFeatsHan','-v7.3');
        save([dest,filesep,'MeanFeatElbow','.mat'],'MeanFeatsElb','-v7.3');
        
        
        
%% Summary Plots


figure
%plot(AvRest)
stem(spkChs,MeanFeatsHan(1,:),'LineStyle','--','Marker','*')
hold on
%plot(AvFlex)
stem(spkChs,MeanFeatsHan(2,:),'LineStyle','--','Marker','square')
%plot(AvExtend)
stem(spkChs,MeanFeatsHan(3,:),'LineStyle','--','Marker','diamond')
title(['Average Spike counts in ',num2str(spkTWs), 's time bins']);
xlabel('Channels');
ylabel('Spike Count in time bins')
legend('Hand Hold','Hand Flexion','Hand Extension');
hold off


figure
stem(spkChs,MeanFeatsElb(1,:),'LineStyle','--','Marker','*')
hold on
%plot(AvFlex)
stem(spkChs,MeanFeatsElb(2,:),'LineStyle','--','Marker','square')
%plot(AvExtend)
stem(spkChs,MeanFeatsElb(3,:),'LineStyle','--','Marker','diamond')
title(['Average Spike counts in ',num2str(spkTWs), 's time bins']);
xlabel('Channels');
ylabel('Spike Count in time bins')
legend('Elbow Hold','Elbow Flexion','Elbow Extension');
hold off



%%
% STEP 3
% This might take some time depending on the time necessary to train the
% classification models. Reducing the number of features will speed up this
% step as well.
%Find_Optimal_Classifiers_1(isCont,ThisClassifier,1);


%% Look for the best features
%FindBestFeats(spkChs,OnlyUseTheseLabels);


%% Look for the optimal classifiers using the best features
%Find_Optimal_Classifiers_reducedFeats();

