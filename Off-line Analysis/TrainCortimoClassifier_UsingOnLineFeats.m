
% OFFLINE ANALYSIS BAED ON FEARTURES EXTRACTED AT CORTIMO RUN TIME


function [hh] = TrainCortimoClassifier_UsingOnLineFeats(ff) 

% This script runs all the necessary functions to analyze the raw
% Blackrock files, synchronize them with the corresponding Cortimo training
% data files, extract a set of features use them to train a set of
% classification models. Exisiting files will be overwritten.

% 1) Select the Raw Blackrock data files for a specific session. Select the
% corresponding Cortimo training binary files for the specific session.

% 2) Grab all the epochs and generate a set of features files using the 
% Use the user selected parameters to extract specific features and run the
% classification using these parameters.

% 3) Grab all the features .mat files in the Features folder and identify
% the optimal features with the optimal classification algorithms. The
% final results will be stored in the Results folder.

% The classifiers trained and validated are those in the file Classifiers.m
% Different ones can be introduced here

% clearvars; close all; clc;

try
% STEP 1
% This will grab the raw Blackrock files and Cortimo training files,
% synchronize them and generate new files that are the epoched version of
% the continuous raw data ready to be used for feature extraction
progress.msgBox = waitbar(0.05,'PLEASE WAIT','Name','DATA PROCESSING',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

[tLFP,ff,ElbDataL,ElbZdataL,HanDataL,HanZdataL] = OnlineFeatureAnalysis_LFP_3(ff);

[tSPK,fff,ElbDataSS,ElbZdataSS,HanDataSS,HanZdataSS] = OnlineFeatureAnalysis_SPK_3(ff);


% STEP 2
% Get apramteters from usersettings, extract features and get them ready
% for classification training

waitbar(0.6,progress.msgBox);

% 0 is best classifier
% 1 is LDA
% if(ff.usersettings.save_freq == 1)
if(~isempty(HanDataL) || ~isempty(ElbDataL))
    
    %if(ff.usersettings.save_neural_spikes == 1)
     if(~isempty(HanDataS) || ~isempty(ElbDataS))
         
        ElbDataM = [];
        HanDataM = [];
        
         % Elbow times are in column one, while hand is column 2
         %for i=1:size(tSPK,2)
         
         % ELBOW
        [val,SPKi,LFPi] = intersect(tSPK{1,1},...
            tLFP{1,1});
        % BUILD a new feature matrix with the common features between
        % spikes and LFPs that get updated at the speed of the slowest time
        % update, given by the LFPs and introduce a new summary metric
        % which accounts for the spike variability within that longer time
        % windows containing a few original spike windows
        step = round(ff.usersettings.EEGtimewin / ff.usersettings.spktimewin); 
        sta = find(SPKi-step>0);
        nIDxs = SPKi(sta:length(SPKi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % This implementation takes into account all the spike windows that occur between
% % consecutive slower LFP updates. This might wind up generating too many
% features, so for now we will switch to the average of the spike windows
% %         newF = zeros(length(nIDxs),(size(ElbDataS,2)-1)*step);
% %         for sttt=1:length(nIDxs)
% %         pp = ElbDataS([nIDxs(sttt)-step+1:nIDxs(sttt)],1:size(ElbDataS,2)-1)';
% %         newF(sttt,:) = pp(:);
% %         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newF = zeros(length(nIDxs),size(ElbDataS,2)-1);
    for sttt=1:length(nIDxs)
         pp = ElbDataS([nIDxs(sttt)-step+1:nIDxs(sttt)],1:size(ElbDataS,2)-1)';
          newF(sttt,:) = mean(pp);
    end
        ElbDataM = [ElbDataL(LFPi(sta:length(LFPi)),(1:size(ElbDataL,2)-1)),...
            ElbDataS(nIDxs,1:size(ElbDataS,2)-1),newF,ElbDataL(LFPi(sta:length(LFPi)),size(ElbDataL,2))];
        
        % HAND
        [val,SPKi,LFPi] = intersect(tSPK{1,2},...
            tLFP{1,2});
        % BUILD a new feature matrix with the common features between
        % spikes and LFPs that get updated at the speed of the slowest time
        % update, given by the LFPs and introduce a new summary metric
        % which accounts for the spike variability within that longer time
        % windows containing a few original spike windows
        step = round(ff.usersettings.EEGtimewin / ff.usersettings.spktimewin); 
        sta = find(SPKi-step>0);
        nIDxs = SPKi(sta:length(SPKi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% % This implementation takes into account all the spike windows that occur between
% % consecutive slower LFP updates. This might wind up generating too many
% features, so for now we will switch to the average of the spike windows
% %         newF = zeros(length(nIDxs),(size(HanDataS,2)-1)*step);
% %         for sttt=1:length(nIDxs)
% %         pp = HanDataS([nIDxs(sttt)-step+1:nIDxs(sttt)],1:size(HanDataS,2)-1)';
% %         newF(sttt,:) = pp(:);
% %         end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
newF = zeros(length(nIDxs),size(HanDataS,2)-1);
    for sttt=1:length(nIDxs)
         pp = HanDataS([nIDxs(sttt)-step+1:nIDxs(sttt)],1:size(HanDataS,2)-1)';
          newF(sttt,:) = mean(pp);
    end
        HanDataM = [HanDataL(LFPi(sta:length(LFPi)),(1:size(HanDataL,2)-1)),...
           HanDataS(nIDxs,1:size(HanDataS,2)-1),newF,HanDataL(LFPi(sta:length(LFPi)),size(HanDataL,2))];
      
        ff = Find_Optimal_Classifiers_onlinefeats(ff,1,HanDataM,ElbDataM);         
      
        ff.usersettings.save_neural_spikes = 1;
        ff.usersettings.save_freq = 1;
        ff.usersettings.Classifiers.isCombineFeatures = 1;
        else
        
        ff = Find_Optimal_Classifiers_onlinefeats(ff,1,HanDataL,ElbDataL);
        ff.usersettings.save_freq = 1;

     end
    
else
    % if(ff.usersettings.save_neural_spikes == 1)
    if(~isempty(HanDataS) || ~isempty(ElbDataS))
       ff = Find_Optimal_Classifiers_onlinefeats(ff,1,HanDataS,ElbDataS);
        ff.usersettings.save_neural_spikes = 1;

    end
    
end
% ff = Find_Optimal_Classifiers_onlinefeats(ff,1,ElbData);

delete(progress.msgBox);

catch ME
    disp(ME.message)
    delete(progress.msgBox)
end
    
hh=ff;


