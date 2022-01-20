% OFFLINE ANALYSIS
% This can be run to identify optimal features in data set. It takes
% advantage of the Cortimo GUI for feature selcetion and automatically
% exports data to Cortimo, for real-time use


function [hh] = TrainCortimoClassifier_1(ff,Mods_to_train)


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
    
%      [progress, is2D] = RawDataToEpochs_3(ff.usersettings.DiscreteContinuous, ... 
%          ff.usersettings.CondLabels.ContinuousVR, ...
%          ff.usersettings.CondLabels.Continuous2D);
%     
     [progress, is2D] = RawDataToEpochs_4(ff.usersettings.DiscreteContinuous, ... 
        ff.usersettings.CondLabels.ContinuousVR, ...
        ff.usersettings.CondLabels.Continuous2D,...
        ff.usersettings.CondLabels.HandLabels,ff.usersettings.CondLabels.ElbowLabels);
    
    
    isCont = 0;
    
    if strcmp(ff.usersettings.DiscreteContinuous,"Continuous")
        isCont = 1;
    end
    
    
    % STEP 2
    % Get parameters from usersettings, extract features and get them ready
    % for classification training
    waitbar(0.35,progress.msgBox)
    
    dt = 0;
    
    if(ff.usersettings.Classifiers.isCombineFeatures == 1 && ...
            ff.usersettings.save_neural_spikes == 1 && ...
            ff.usersettings.save_freq == 1)
        
        if(~isempty(ff.usersettings.CondLabels.HandLabels) && ...
               ~isempty(ff.usersettings.CondLabels.ElbowLabels) )
           
            EpochsToTheseFeaturesALL_2(str2num(ff.usersettings.ch_to_spi),...
            ff.usersettings.spktimewin/1000,str2num(ff.usersettings.ch_to_freq),...
            ff.usersettings.EEGFFTWin, ...
            (1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000)),...
            ff.usersettings.Fstep,ff.usersettings.Frange, isCont, ff.usersettings.CondLabels.HandLabels, ...
            ff.usersettings.CondLabels.ElbowLabels);
        
        else
            if(~isempty(ff.usersettings.CondLabels.HandLabels))
                EpochsToTheseFeaturesALL_2(str2num(ff.usersettings.ch_to_spi),...
                ff.usersettings.spktimewin/1000,str2num(ff.usersettings.ch_to_freq),...
                ff.usersettings.EEGFFTWin, ...
                (1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000)),...
                ff.usersettings.Fstep,ff.usersettings.Frange, isCont, ff.usersettings.CondLabels.HandLabels, ...
                [0,1,2]);
            else
                 if(~isempty(ff.usersettings.CondLabels.ElbowLabels))
                      EpochsToTheseFeaturesALL_2(str2num(ff.usersettings.ch_to_spi),...
                       ff.usersettings.spktimewin/1000,str2num(ff.usersettings.ch_to_freq),...
                       ff.usersettings.EEGFFTWin, ...
                       (1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000)),...
                       ff.usersettings.Fstep,ff.usersettings.Frange, isCont, [0,1,2], ...
                       ff.usersettings.CondLabels.ElbowLabels);
                 else
                     
                      EpochsToTheseFeaturesALL_2(str2num(ff.usersettings.ch_to_spi),...
                       ff.usersettings.spktimewin/1000,str2num(ff.usersettings.ch_to_freq),...
                       ff.usersettings.EEGFFTWin, ...
                       (1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000)),...
                       ff.usersettings.Fstep,ff.usersettings.Frange, isCont, [0,1,2], ...
                       [0,1,2]);
                 end
            end
        end
        dt = ff.usersettings.EEGtimewin / 1000;
    else
        if(ff.usersettings.save_freq == 1)
            if(~isempty(ff.usersettings.CondLabels.HandLabels) && ...
                  ~isempty(ff.usersettings.CondLabels.HandLabels)  )
              
            EpochsToTheseFeaturesLFP(str2num(ff.usersettings.ch_to_freq),...
                ff.usersettings.EEGFFTWin,...
                1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000),...
                ff.usersettings.Fstep,ff.usersettings.Frange, isCont, ff.usersettings.CondLabels.HandLabels, ...
                ff.usersettings.CondLabels.ElbowLabels);
            else
                if(~isempty(ff.usersettings.CondLabels.HandLabels))
                EpochsToTheseFeaturesLFP(str2num(ff.usersettings.ch_to_freq),...
                ff.usersettings.EEGFFTWin,...
                1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000),...
                ff.usersettings.Fstep,ff.usersettings.Frange, isCont, ff.usersettings.CondLabels.HandLabels, ...
                [0,1,2]);
            
                else
                    if(~isempty(ff.usersettings.CondLabels.ElbowLabels))
                    EpochsToTheseFeaturesLFP(str2num(ff.usersettings.ch_to_freq),...
                    ff.usersettings.EEGFFTWin,...
                    1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000),...
                    ff.usersettings.Fstep,ff.usersettings.Frange, isCont, [0,1,2], ...
                    ff.usersettings.CondLabels.ElbowLabels);
                    else
                    EpochsToTheseFeaturesLFP(str2num(ff.usersettings.ch_to_freq),...
                    ff.usersettings.EEGFFTWin,...
                    1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000),...
                    ff.usersettings.Fstep,ff.usersettings.Frange, isCont, [0,1,2], ...
                    [0,1,2]);
                    end
                    
                end
            end
            
            dt = ff.usersettings.EEGtimewin / 1000;
        end
        
        if(ff.usersettings.save_neural_spikes == 1)
            if(~isempty(ff.usersettings.CondLabels.HandLabels) && ...
                  ~isempty(ff.usersettings.CondLabels.ElbowLabels)  )
              
            EpochsToTheseFeaturesSPK(str2num(ff.usersettings.ch_to_spi),...
                ff.usersettings.spktimewin /1000, isCont, ff.usersettings.CondLabels.HandLabels, ...
                ff.usersettings.CondLabels.ElbowLabels,ff.usersettings.isLeakyInt,...
                ff.usersettings.NormSpikes);
            else
                if(~isempty(ff.usersettings.CondLabels.HandLabels))
                 EpochsToTheseFeaturesSPK(str2num(ff.usersettings.ch_to_spi),...
                ff.usersettings.spktimewin /1000, isCont, ff.usersettings.CondLabels.HandLabels, ...
                [0,1,2],ff.usersettings.isLeakyInt,...
                ff.usersettings.NormSpikes);
                else
                    if(~isempty(ff.usersettings.CondLabels.ElbowLabels))
                       EpochsToTheseFeaturesSPK(str2num(ff.usersettings.ch_to_spi),...
                        ff.usersettings.spktimewin /1000, isCont, [0,1,2] , ...
                        ff.usersettings.CondLabels.ElbowLabels,ff.usersettings.isLeakyInt,...
                        ff.usersettings.NormSpikes);  
                        
                    else
                         EpochsToTheseFeaturesSPK(str2num(ff.usersettings.ch_to_spi),...
                        ff.usersettings.spktimewin /1000, isCont, [0,1,2] , ...
                        [0,1,2],ff.usersettings.isLeakyInt,...
                        ff.usersettings.NormSpikes);  
                    end
                    
                end
                
            end
            
            dt = ff.usersettings.spktimewin / 1000;
        end
        
    end
    
    
    % STEP 2
    % This might take some time depending on the time necessary to train the
    % classification models. Reducing the number of features will speed up this
    % step as well.
    waitbar(0.7,progress.msgBox)
    
    ff = Find_Optimal_Classifiers_3(ff,isCont,dt);
    
    delete(progress.msgBox);
    
catch ME
    disp(ME.message);
    if(exist('progress','var'))
    delete(progress.msgBox);
    end
end

hh=ff;

