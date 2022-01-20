function [progress, is2D] = RawDataToEpochs_Multiple_4(DiscCont, ContinuousLabelVR, ContinuousLabel2D, ...
    DiscLabelsHand, DiscLabelsElbow,filename,fileList)

% clearvars; close all; clc;
addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('Aux Funcs'));
% rmpath(['eeglab',filesep,'plugins']);
%% Open NS file containing Blackrock data files
% Stores continuous data
% [file,path] = uigetfile('*.*','Select raw data file');
% filename = fullfile(path,file);
% %%%%%%%%%%%%%%%%%%%%%%%
% if file==0
%     % user pressed cancel
%     return;
% end
% %%%%%%%%%%%%%%%%%%%%%%

%openNSx('report','read',filename,'sample', 'p:short', 's:75');
% openNSx('report','read',filename,'sample', 'p:short');

% Open corresponding NEV file to load the events (extracted spikes)
NEV = openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');
comments = NEV.Data.Comments;
events = double(NEV.Data.Spikes.TimeStamp);
electrode = double(NEV.Data.Spikes.Electrode);
Unit = NEV.Data.Spikes.Unit;

% If a continuous file has been loaded
if(contains(filename(end-3:end),'ns'))
    
    % Open NES file
    openNSx('report','read','uV',filename,'t:0:40','min','p:short');
    % openNSx('report','read','uV',filename,'p:short','s:2');
    
    % openNSx('report','read','uV',filename,'p:short');
    
    % Depending on the data format the structure names are different
    forma = filename(end);
    %rawData = NS5.Data;
    tp = eval(['NS',forma]);
    if iscell(tp.Data) == 1
        rawData = cat(2,tp.Data{:});
    else
        rawData = tp.Data;
    end
    
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
    
    clear tp;
    
    
%     %% Load the training file(s)
%     [fileList,pathT] = uigetfile(strcat(path,'*.bin'),'Select file with training data',...
%         'MultiSelect','on');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%
%     if fileList==0
%         % user pressed cancel
%         return;
%     end
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    if(iscell(fileList))
        LEN = length(fileList);
        
    else
        LEN = 1;
        tp{1} = fileList;
        fileList = tp;
    end
    
    
    
    progress.msgBox = waitbar(0.05,'PLEASE WAIT','Name','DATA PROCESSING',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    is2D = 0;
    
    % Concatenate all the training data files
    VRTrainingLabels = [];
    for file=1:LEN
        
        fileT = fileList{file};
        
        % filenameFeat = fullfile(pathT,fileT);
        filenameTrain = fullfile(fileT);
        
        % Open both files, the training data file and the EEG feature file
        % filenameTrain = strrep(filenameFeat,"FreqFeatures","ArmVR");
        try
            if contains(filenameTrain, "CenterOut")
                [VRTrainingLabels] = [VRTrainingLabels;loadCenterOutTrainingFile(filenameTrain)];
                is2D = 1;
            else
                [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile(filenameTrain)];
                is2D = 0;
            end
            % [EEGfeats] = loadEEGFeatureFile(filenameFeat);
        catch
            [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile_noDecoderOutput(filenameTrain)];
            is2D = 0;
            
        end
        
    end
    
    
    %  *********************************************************************  %
    %  THE FOLLOWING INFORMATION IS IMPORTANT TO INTERPRET THE BINARY DATA
    %  STREAM RECORDED FROM TRAINING/TESTING RT ACQUISITIONS USING CORTIMO BCI
    %  SYSTEM
    %
    %
    % The training binary file structure is as follows:
    % 15 doubles:
    % 1) Matlab acquisition cycle number (since last start/stop)
    % 2) time stamp from the Matlab runtime execution
    % 3) time stamp of the beginning of the last data sample read from the NSP
    % 4) Training/Testing BCI trial number if VR application is used
    % 5) & 6) Respectively wrist and elbow TARGET position if/when present
    % 7) & 8) Respectively wrist and elbow label based on real-time feedback
    % from the VR app. In other words, these indicate the "correction" command
    % that should be sent to the "ARM" (either MyoPro or VR) to hit the target
    % position
    % 9) & 10) respectively wrist and elbow target hit results. When there is a
    % goal oriented paradigm in progress, a 1 will indeicate whether the arm
    % segments are close to the target position and get there within the
    % specified trial duration
    % 11) & 12) respectively wrist and elbow positions of the VR RT feedback
    % 13) & 14) respectively wrist and elbow position of the MyoPro brace when
    % connected to the system
    
    %  *********************************************************************  %
    
    % Find start of the training file
    %traningStartID = find(round(VRTrainingLabels{:,3}(1),4) == round(rawtime,4));
    %trainingEndID = find(round(VRTrainingLabels{:,3}(end),4) == round(rawtime,4));
    
    % Synchronize training events and Matlab Cortimo cycles
    %TT =  SynchTrainingLabels(VRTrainingLabels);
    % This new synchronized time stamp array can be used to identify events in
    % the rawtime array and thus process the raw data points.
    % find(round(TT.NSPTime(3033),6) == round(rawtime,6))
    
    %% Group the data in epochs or trials or conditions
    %* There are two data sterams that need to be epoched and grouped by trials.
    %* 1) The raw voltage traces; 2) The extracted spikes from the NSP
    
    %
    % Approximate time stamps with Matlab system
    dd = find(diff(VRTrainingLabels.MatlabCycle)~=0);
    newTraining = VRTrainingLabels(dd,:);
    
    % Find start of the training file
    %traningStartID = find(round(rawtime,5)==round(newTraining{:,3}(1),5));
    %trainingEndID = find(round(rawtime,5)==round(newTraining{:,3}(end),5));
    
    fTime = round(rawtime,5);
    newTraining.NSPTime = round(newTraining.NSPTime,5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRAINING conditions can change based on patient specific training
    % paradigms. The deault training labels for function restoration will be
    % 0,1,2 respectively for hold, flex and extend. These might need to be
    % changed in case the training trials are modified.
    
    
    if (strcmp(DiscCont,"Discrete"))
        
        if(~isempty(DiscLabelsHand))
            %%%%%% HAND %%%%%%
            % REST correpsonds to WristTarget = -1 with no target
            % HOLD correpsonds to WristLab = 0;
            % FLEX correponds to WristLab = 1;
            % EXTEND correpsonds to WristLab = 2;
            conds = [0,1,2];
            HanEpochVolts = cell(4,1);
            HanEpochSpks = cell(4,1);
            HanEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                HanEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere*3) & events<(stophere*3));
                HanEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                HanEpochTS{1,co} = 3 .* [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                idcond = find(newTraining.WriLab == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                
                if(~isempty(A))
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    % The intervals between are to be used
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        HanEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere*3) & events<(stophere*3));
                        HanEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        HanEpochTS{cond+1,co} = 3 .* [starthere,stophere];
                        
                        co = co + 1;
                    end
                    
                end
                
            end
            
            waitbar(0.2,progress.msgBox);
            
            
            
            
            
        else
            
            DerivedLabels = [];
            % trialInds = find(newTraining.Trial > 0.1);
            
%             if strcmp(ContinuousLabelVR,"Target Position")
%                 DerivedLabels = [0;sign(diff(newTraining.WriTarget))];
%             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
%                 DerivedLabels = [0;sign(diff(newTraining.VRWriPos))];
%             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
%                 DerivedLabels = [0;sign(diff(newTraining.MyoProWriPos))];
%             end
%             
%             DerivedLabels(DerivedLabels == -1) = 2;
%             
             if strcmp(ContinuousLabelVR,"Target Position")
                 DerivedLabels = zeros(size(newTraining.WriTarget));
                 ddd = find(newTraining.WriTarget > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.WriTarget < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                 DerivedLabels = zeros(size(newTraining.VRWriPos));
                 ddd = find(newTraining.VRWriPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.VRWriPos < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                 DerivedLabels = zeros(size(newTraining.MyoProWriPos));
                 ddd = find(newTraining.MyoProWriPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.MyoProWriPos < 45);
                 DerivedLabels(ddd) = 2;

             end
            
            
            conds = [0,1,2];
            HanEpochVolts = cell(4,1);
            HanEpochSpks = cell(4,1);
            HanEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0, 200);
            
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                HanEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere*3) & events<(stophere*3));
                HanEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                HanEpochTS{1,co} = 3 .* [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                % idcond = find(newTraining.WriLab == conds(cond));
                idcond = find(DerivedLabels == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                
                if(~isempty(A))
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    % The intervals between are to be used
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        HanEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere*3) & events<(stophere*3));
                        HanEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        HanEpochTS{cond+1,co} = 3 .* [starthere,stophere];
                        
                        co = co + 1;
                    end
                    
                end
                
            end
            
            
            
        end
        
        
        
        if(~isempty(DiscLabelsElbow))
            %%%%%% ELBOW %%%%%%
            % REST correpsonds to Target = -1 with no target
            % HOLD correpsonds to WristLab = 0;
            % FLEX correponds to WristLab = 1;
            % EXTEND correpsonds to WristLab = 2;
            conds = [0,1,2];
            ElbEpochVolts = cell(4,1);
            ElbEpochSpks = cell(4,1);
            ElbEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                ElbEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere*3) & events<(stophere*3));
                ElbEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                ElbEpochTS{1,co} = 3 .* [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                idcond = find(newTraining.ElbLab == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                if(~isempty(A))
                    
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    
                    % The intervals between the indices in timeTouse are to be used
                    % The above steps detect training trials that are longer than 1 second. If
                    % a trial is shorter than that we can safely leave it out.
                    
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        ElbEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere*3) & events<(stophere*3));
                        ElbEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        ElbEpochTS{cond+1,co} = 3 .* [starthere,stophere];
                        co = co + 1;
                    end
                    
                end
                
            end
            
            
        else
            
%             DerivedLabels = [];
%             % trialInds = find(newTraining.Trial > 0.1);
%             
%             if strcmp(ContinuousLabelVR,"Target Position")
%                 DerivedLabels = [0;sign(diff(newTraining.ElbTarget))];
%             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
%                 DerivedLabels = [0;sign(diff(newTraining.VRElbPos))];
%             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
%                 DerivedLabels = [0;sign(diff(newTraining.MyoProElbPos))];
%             end
%             
%             
%             DerivedLabels(DerivedLabels == -1) = 2;
%             
            if strcmp(ContinuousLabelVR,"Target Position")
                 DerivedLabels = zeros(size(newTraining.ElbTarget));
                 ddd = find(newTraining.ElbTarget > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.ElbTarget < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                 DerivedLabels = zeros(size(newTraining.VRElbPos));
                 ddd = find(newTraining.WriElbPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.VRElbPos < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                 DerivedLabels = zeros(size(newTraining.MyoProElbPos));
                 ddd = find(newTraining.MyoProElbPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.MyoProElbPos < 45);
                 DerivedLabels(ddd) = 2;

             end
             
            
            waitbar(0.25,progress.msgBox);
            
            %%%%%% ELBOW %%%%%%
            % REST correpsonds to Target = -1 with no target
            % HOLD correpsonds to WristLab = 0;
            % FLEX correponds to WristLab = 1;
            % EXTEND correpsonds to WristLab = 2;
            conds = [0,1,2];
            ElbEpochVolts = cell(4,1);
            ElbEpochSpks = cell(4,1);
            ElbEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0,200);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                ElbEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere*3) & events<(stophere*3));
                ElbEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                ElbEpochTS{1,co} = 3 .* [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                % idcond = find(newTraining.ElbLab == conds(cond));
                idcond = find(DerivedLabels == conds(cond));
                
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                if(~isempty(A))
                    
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    
                    % The intervals between the indices in timeTouse are to be used
                    % The above steps detect training trials that are longer than 1 second. If
                    % a trial is shorter than that we can safely leave it out.
                    
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        ElbEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere*3) & events<(stophere*3));
                        ElbEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        ElbEpochTS{cond+1,co} = 3 .* [starthere,stophere];
                        co = co + 1;
                    end
                    
                end
                
            end
            
            
            
        end
        
        
        
        
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    else
        
        
        %%%%%% CONTINUOUS EPOCHS %%%%%%
        %%
        
        LabelsRaw = [];
        
        if is2D == 1
            trialInds = find(newTraining.SignalType ~= 0);
            
            if strcmp(ContinuousLabel2D,"Ball Location")
                LabelsRaw = [newTraining.BallPosX, newTraining.BallPosY];
            elseif strcmp(ContinuousLabel2D,"Target Location")
                LabelsRaw = [newTraining.TargetPosX, newTraining.TargetPosY];
            elseif strcmp(ContinuousLabel2D,"Ball Direction")
                LabelsRaw = [newTraining.BallPosX, newTraining.BallPosY];
            end
            
        else
            trialInds = find(newTraining.Trial > 0.1);
            
            if strcmp(ContinuousLabelVR,"Target Position")
                LabelsRaw = [newTraining.WriTarget, newTraining.ElbTarget];
            elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                LabelsRaw = [newTraining.VRWriPos, newTraining.VRElbPos];
            elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                LabelsRaw = [newTraining.MyoProWriPos, newTraining.MyoProElbPos];
            end
        end
        
        trialEndingInds = find(diff(trialInds) > 25);
        
        if isempty(trialEndingInds)
            trialIndCell{1,1} = trialInds(1:end);
        else
            trialIndCell = cell(length(trialEndingInds)+1, 1);
            
            % Build a cell of arrays where each array contains the indicies of each trial
            trialIndCell{1,1} = trialInds(1:trialEndingInds(1));
            for i=2:length(trialEndingInds)
                trialIndCell{i,1} = trialInds(trialEndingInds(i-1)+1:trialEndingInds(i));
            end
            trialIndCell{length(trialEndingInds)+1,1} = trialInds(trialEndingInds(end)+1:end);
        end
        

        
        ContLabels = cell(length(trialEndingInds)+1, 1);
        ContVolts = cell(length(trialEndingInds)+1, 1);
        ContSpks = cell(length(trialEndingInds)+1, 1);
        ContTS = cell(length(trialEndingInds)+1, 1);
        
        for i=1:length(trialIndCell)
            
            ContLabels{i,1} = [LabelsRaw(trialIndCell{i,1},1), LabelsRaw(trialIndCell{i,1},2)];
            
            [A,~] = intersect(fTime,newTraining.NSPTime(trialIndCell{i,1}));
            [~,rawtimeInds] = find(fTime == A);
            ContVolts{i,1} = rawData(:,rawtimeInds(1):rawtimeInds(end));
            
            r = 30000/Fs;
            idSp = find(events>=(rawtimeInds(1)*r) & events<(rawtimeInds(end)*r));
            ContSpks{i,1} = [double(electrode(idSp));double(events(idSp))];
            ContTS{i,1} = r.*[rawtimeInds(1),rawtimeInds(end)];
        end
        
        if ~isempty(trialEndingInds)
            
            %Get trial Inds INCLUDING dropped frames
            AllInclusiveTrialInds = [];
            for i=1:length(trialIndCell)
                AllInclusiveTrialInds = [AllInclusiveTrialInds;(trialIndCell{i,1}(1):trialIndCell{i,1}(end))'];
            end
            % Store indicies associated with rest (everything BUT trial indicies) in array.
            restInds = (1:height(newTraining))';
            restInds(AllInclusiveTrialInds) = [];
            restEndingInds = find(diff(restInds) > 25);
            restIndCell = cell(length(restEndingInds), 1);
            
            % The first start period begins three seconds before the beginning of
            % the first trial
            firstRestStart = find(round(newTraining.NSPTime,1) == round(newTraining.NSPTime(restEndingInds(1)) - 3.0,1), 1);
            
            restIndCell{1,1} = restInds(firstRestStart:restEndingInds(1));
            for i=2:length(restIndCell)
                
                % This checks if the rest time is longer than 4 seconds, indicating
                % the beginning of a new series of trials. In this case it only
                % takes takes the final 3 seconds of rest
                if (newTraining.NSPTime(restInds(restEndingInds(i))) - newTraining.NSPTime(restInds(restEndingInds(i-1)+1)) > 4 ...
                        && is2D == 1)
                    firstRestStart = find(round(newTraining.NSPTime,2) == round(newTraining.NSPTime(restEndingInds(i)),2) - 3, 1);
                    restIndCell{i,1} = restInds(firstRestStart:restEndingInds(i));
                else
                    % Otherwise use the entire rest period
                    restIndCell{i,1} = restInds(restEndingInds(i-1)+1:restEndingInds(i));
                end
            end
            
            RestContVolts = cell(length(restEndingInds), 1);
            RestContSpks = cell(length(restEndingInds), 1);
            RestContTS = cell(length(restEndingInds), 1);
            
            for i=1:length(trialIndCell)
                [A,~] = intersect(fTime,newTraining.NSPTime(restIndCell{i,1}));
                [~,rawtimeInds] = find(fTime == A);
                RestContVolts{i,1} = rawData(:,rawtimeInds(1):rawtimeInds(end));
                
                r = 30000/Fs;
                idSp = find(events>=(rawtimeInds(1)*r) & events<(rawtimeInds(end)*r));
                RestContSpks{i,1} = [double(electrode(idSp));double(events(idSp))];
                RestContTS{i,1} = r.*[rawtimeInds(1),rawtimeInds(end)];
            end
        end
    end
    % HanEpochVolts and ElbEpochVolts now contain the raw data grouped according to the
    % different trial conditions. These signals can now be used to derive
    % features and susequently feed those features into the machine learning
    % training system
    
    % HanEpochSpks and ElbEpochSpks now contain the spike extracted from the
    % NSP and epoched follwoing the same criterion used for the raw voltages.
    % These spike time stamps need to be broked in specific time bins to
    % extract spike features in next processing step.
    
    %-------------------------------------------------------------------------%
    % Export the Epoched data
    [~,nm,~] = fileparts(filename);
    dest = 'Epochs';
    if(exist(dest,'dir') ~= 7)
        mkdir(dest);
    end
    
    % Export all data to files to be used in subsequent analysis
    if strcmp(DiscCont,"Discrete")
        save([dest,filesep,nm,'-HanEpochVolts','.mat'],'HanEpochVolts','-v7.3');
        save([dest,filesep,nm,'-ElbEpochVolts','.mat'],'ElbEpochVolts','-v7.3');
        save([dest,filesep,nm,'-HanEpochSpks','.mat'],'HanEpochSpks','-v7.3');
        save([dest,filesep,nm,'-ElbEpochSpks','.mat'],'ElbEpochSpks','-v7.3');
        save([dest,filesep,nm,'-HanEpochTS','.mat'],'HanEpochTS','-v7.3');
        save([dest,filesep,nm,'-ElbEpochTS','.mat'],'ElbEpochTS','-v7.3');
    else
        save([dest,filesep,nm,'-ContLabels','.mat'],'ContLabels','-v7.3');
        save([dest,filesep,nm,'-ContVolts','.mat'],'ContVolts','-v7.3');
        save([dest,filesep,nm,'-ContSpks','.mat'],'ContSpks','-v7.3');
        save([dest,filesep,nm,'-ContTS','.mat'],'ContTS','-v7.3');
        if ~isempty(trialEndingInds)
            save([dest,filesep,nm,'-RestContVolts','.mat'],'RestContVolts','-v7.3');
            save([dest,filesep,nm,'-RestContSpks','.mat'],'RestContSpks','-v7.3');
            save([dest,filesep,nm,'-RestContTS','.mat'],'RestContTS','-v7.3');
        end
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else % No continuous files
    
    
    
    
    %% Load the training file(s)
%     [fileList,pathT] = uigetfile(strcat(path,'*.bin'),'Select file with training data',...
%         'MultiSelect','on');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%
%     if fileList==0
%         % user pressed cancel
%         return;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%
    
    
    if(iscell(fileList))
        LEN = length(fileList);
        
    else
        LEN = 1;
        tp{1} = fileList;
        fileList = tp;
    end
    
    
    
    progress.msgBox = waitbar(0.05,'PLEASE WAIT','Name','DATA PROCESSING',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    is2D = 0;
    
    % Concatenate all the training data files
    VRTrainingLabels = [];
    for file=1:LEN
        
        fileT = fileList{file};
        
        % filenameFeat = fullfile(pathT,fileT);
        filenameTrain = fullfile(fileT);
        
        % Open both files, the training data file and the EEG feature file
        % filenameTrain = strrep(filenameFeat,"FreqFeatures","ArmVR");
        try
            if contains(filenameTrain, "CenterOut")
                [VRTrainingLabels] = [VRTrainingLabels;loadCenterOutTrainingFile(filenameTrain)];
                is2D = 1;
            else
                [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile(filenameTrain)];
                is2D = 0;
            end
            % [EEGfeats] = loadEEGFeatureFile(filenameFeat);
        catch
            [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile_noDecoderOutput(filenameTrain)];
            is2D = 0;
            
        end
        
    end
    

%  *********************************************************************  %
%  THE FOLLOWING INFORMATION IS IMPORTANT TO INTERPRET THE BINARY DATA
%  STREAM RECORDED FROM TRAINING/TESTING RT ACQUISITIONS USING CORTIMO BCI
%  SYSTEM
%
%
% The training binary file structure is as follows:
% 15 doubles:
% 1) Matlab acquisition cycle number (since last start/stop)
% 2) time stamp from the Matlab runtime execution
% 3) time stamp of the beginning of the last data sample read from the NSP
% 4) Training/Testing BCI trial number if VR application is used
% 5) & 6) Respectively wrist and elbow TARGET position if/when present
% 7) & 8) Respectively wrist and elbow label based on real-time feedback
% from the VR app. In other words, these indicate the "correction" command
% that should be sent to the "ARM" (either MyoPro or VR) to hit the target
% position
% 9) & 10) respectively wrist and elbow target hit results. When there is a
% goal oriented paradigm in progress, a 1 will indeicate whether the arm
% segments are close to the target position and get there within the
% specified trial duration
% 11) & 12) respectively wrist and elbow positions of the VR RT feedback
% 13) & 14) respectively wrist and elbow position of the MyoPro brace when
% connected to the system

%  *********************************************************************  %

% Find start of the training file
% traningStartID = find(round(VRTrainingLabels{:,3}(1),4) == round(rawtime,4));
% trainingEndID = find(round(VRTrainingLabels{:,3}(end),4) == round(rawtime,4));
% 
% % Synchronize training events and Matlab Cortimo cycles
% TT =  SynchTrainingLabels(VRTrainingLabels);
% This new synchronized time stamp array can be used to identify events in
% the rawtime array and thus process the raw data points.
% find(round(TT.NSPTime(3033),6) == round(rawtime,6))

%% Group the data in epochs or trials or conditions
%* There are two data sterams that need to be epoched and grouped by trials.
%* 1) The raw voltage traces; 2) The extracted spikes from the NSP

%
% Approximate time stamps with Matlab system
dd = find(diff(VRTrainingLabels.MatlabCycle)~=0);
newTraining = VRTrainingLabels(dd,:);

% Find start of the training file
%traningStartID = find(round(rawtime,5)==round(newTraining{:,3}(1),5));
%trainingEndID = find(round(rawtime,5)==round(newTraining{:,3}(end),5));

%fTime = round(rawtime,5);
ev = events(1):events(length(events));
fTime = round(ev/30e3,5);

newTraining.NSPTime = round(newTraining.NSPTime,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAINING conditions can change based on patient specific training
% paradigms. The deault training labels for function restoration will be
% 0,1,2 respectively for hold, flex and extend. These might need to be
% changed in case the training trials are modified.


    if (strcmp(DiscCont,"Discrete"))
    
    if(~isempty(DiscLabelsHand))
    %%%%%% HAND %%%%%%
    % REST correpsonds to WristTarget = -1 with no target
    % HOLD correpsonds to WristLab = 0;
    % FLEX correponds to WristLab = 1;
    % EXTEND correpsonds to WristLab = 2;
    conds = [0,1,2];
   % HanEpochVolts = cell(4,1);
    HanEpochSpks = cell(4,1);
    HanEpochTS = cell(4,1);
    
    % idREST = find(newTraining.WriTarget == -1);
    idREST = find(newTraining.Trial == 0);
    [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
    
    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
    % The intervals between are to be used
    co = 1;
    for interval=1:length(timeTouse)-1
        starthere = find(fTime == A(timeTouse(interval)+1));    
        stophere = find(fTime == A(timeTouse(interval+1)));
        % HanEpochVolts{1,co} = rawData(:,starthere:stophere);
        % Sample numbers represented here by starthere and stophere, need to be
        % adjusted to accoutn for the fact that the spike events are stored as
        % sample numbers @30KHz not @10KHz (this only applies to the continuous
        % signals).
        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
        idSp = find(events>=(starthere) & events<(stophere));
        HanEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
        HanEpochTS{1,co} =  [starthere,stophere];
        co = co + 1;
    end
  
   
            for cond=1:length(conds)
                idcond = find(newTraining.WriLab == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                
                if(~isempty(A))
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    % The intervals between are to be used
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        %  HanEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere) & events<(stophere));
                        HanEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        HanEpochTS{cond+1,co} =  [starthere,stophere];
                        
                        co = co + 1;
                    end
                    
                end
                
            end
            
            waitbar(0.2,progress.msgBox);
            
            
            
    else
            
            DerivedLabels = [];
            % trialInds = find(newTraining.Trial > 0.1);
%             
%             if strcmp(ContinuousLabelVR,"Target Position")
%                 DerivedLabels = [0;sign(diff(newTraining.WriTarget))];
%             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
%                 DerivedLabels = [0;sign(diff(newTraining.VRWriPos))];
%             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
%                 DerivedLabels = [0;sign(diff(newTraining.MyoProWriPos))];
%             end
%             
%             DerivedLabels(DerivedLabels == -1) = 2;
%             
                      
          if strcmp(ContinuousLabelVR,"Target Position")
                 DerivedLabels = zeros(size(newTraining.WriTarget));
                 ddd = find(newTraining.WriTarget > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.WriTarget < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                 DerivedLabels = zeros(size(newTraining.VRWriPos));
                 ddd = find(newTraining.VRWriPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.VRWriPos < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                 DerivedLabels = zeros(size(newTraining.MyoProWriPos));
                 ddd = find(newTraining.MyoProWriPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.MyoProWriPos < 45);
                 DerivedLabels(ddd) = 2;

          end
             
            
            conds = [0,1,2];
            % HanEpochVolts = cell(4,1);
            HanEpochSpks = cell(4,1);
            HanEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0,200);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                %    HanEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere) & events<(stophere));
                HanEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                HanEpochTS{1,co} = [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                % idcond = find(newTraining.WriLab == conds(cond));
                idcond = find(DerivedLabels == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                
                if(~isempty(A))
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    % The intervals between are to be used
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        %         HanEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere) & events<(stophere));
                        HanEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        HanEpochTS{cond+1,co} = [starthere,stophere];
                        
                        co = co + 1;
                    end
                    
                end
                
            end
            
    end
      
        waitbar(0.25,progress.msgBox);
        
        
        
        if(~isempty(DiscLabelsElbow))
            %%%%%% ELBOW %%%%%%
            % REST correpsonds to Target = -1 with no target
            % HOLD correpsonds to WristLab = 0;
            % FLEX correponds to WristLab = 1;
            % EXTEND correpsonds to WristLab = 2;
            conds = [0,1,2];
            % ElbEpochVolts = cell(4,1);
            ElbEpochSpks = cell(4,1);
            ElbEpochTS = cell(4,1);
            
            idREST = find(newTraining.ElbTarget == -1);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                %    ElbEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere) & events<(stophere));
                ElbEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                ElbEpochTS{1,co} = [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                idcond = find(newTraining.ElbLab == conds(cond));
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                if(~isempty(A))
                    
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    
                    % The intervals between the indices in timeTouse are to be used
                    % The above steps detect training trials that are longer than 1 second. If
                    % a trial is shorter than that we can safely leave it out.
                    
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        %    ElbEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere) & events<(stophere));
                        ElbEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        ElbEpochTS{cond+1,co} = [starthere,stophere];
                        co = co + 1;
                    end
                    
                end
                
            end
            
            
        else
            
            DerivedLabels = [];
            % trialInds = find(newTraining.Trial > 0.1);
            
%             if strcmp(ContinuousLabelVR,"Target Position")
%                 DerivedLabels = [0;sign(diff(newTraining.ElbTarget))];
%             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
%                 DerivedLabels = [0;sign(diff(newTraining.VRElbPos))];
%             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
%                 DerivedLabels = [0;sign(diff(newTraining.MyoProElbPos))];
%             end
%             
%             DerivedLabels(DerivedLabels == -1) = 2;
    if strcmp(ContinuousLabelVR,"Target Position")
                 DerivedLabels = zeros(size(newTraining.ElbTarget));
                 ddd = find(newTraining.ElbTarget > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.ElbTarget < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                 DerivedLabels = zeros(size(newTraining.VRElbPos));
                 ddd = find(newTraining.VRElbPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.VRElbPos < 45);
                 DerivedLabels(ddd) = 2;

             elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                 DerivedLabels = zeros(size(newTraining.MyoProElbPos));
                 ddd = find(newTraining.MyoProElbPos > 45);
                 DerivedLabels(ddd) = 1;
                 ddd = find(newTraining.MyoProElbPos < 45);
                 DerivedLabels(ddd) = 2;
                 
    end

            
            waitbar(0.25,progress.msgBox);
            
            %%%%%% ELBOW %%%%%%
            % REST correpsonds to Target = -1 with no target
            % HOLD correpsonds to WristLab = 0;
            % FLEX correponds to WristLab = 1;
            % EXTEND correpsonds to WristLab = 2;
            conds = [0,1,2];
            %  ElbEpochVolts = cell(4,1);
            ElbEpochSpks = cell(4,1);
            ElbEpochTS = cell(4,1);
            
            idREST = find(newTraining.Trial == 0,200);
            [A,~] = intersect(fTime,newTraining.NSPTime(idREST));
            
            timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
            % The intervals in between are to be used
            co = 1;
            for interval=1:length(timeTouse)-1
                starthere = find(fTime == A(timeTouse(interval)+1));
                stophere = find(fTime == A(timeTouse(interval+1)));
                %   ElbEpochVolts{1,co} = rawData(:,starthere:stophere);
                % Sample numbers represented here by starthere and stophere, need to be
                % adjusted to accoutn for the fact that the spike events are stored as
                % sample numbers @30KHz not @10KHz (this only applies to the continuous
                % signals).
                %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                idSp = find(events>=(starthere) & events<(stophere));
                ElbEpochSpks{1,co} = [double(electrode(idSp));double(events(idSp))];
                ElbEpochTS{1,co} = [starthere,stophere];
                co = co + 1;
            end
            
            
            for cond=1:length(conds)
                % idcond = find(newTraining.ElbLab == conds(cond));
                idcond = find(DerivedLabels == conds(cond));
                
                [A,~] = intersect(fTime,newTraining.NSPTime(idcond));
                if(~isempty(A))
                    
                    timeTouse = [0;find(diff(A)>0.2);length(A)]; % if this is true, it means that there was an active trial going on
                    
                    % The intervals between the indices in timeTouse are to be used
                    % The above steps detect training trials that are longer than 1 second. If
                    % a trial is shorter than that we can safely leave it out.
                    
                    co = 1;
                    for interval=1:length(timeTouse)-1
                        starthere = find(fTime == A(timeTouse(interval)+1));
                        stophere = find(fTime == A(timeTouse(interval+1)));
                        %     ElbEpochVolts{cond+1,co} = rawData(:,starthere:stophere);
                        % Sample numbers represented here by starthere and stophere, need to be
                        % adjusted to accoutn for the fact that the spike events are stored as
                        % sample numbers @30KHz not @10KHz (this only applies to the continuous
                        % signals).
                        %%%%%%%%%%% idSp = find(events>=starthere & events<stophere);
                        idSp = find(events>=(starthere) & events<(stophere));
                        ElbEpochSpks{cond+1,co} = [double(electrode(idSp));double(events(idSp))];
                        ElbEpochTS{cond+1,co} = [starthere,stophere];
                        co = co + 1;
                    end
                    
                end
                
            end
            
            
        end
        
        
        
    else
        
        
        %%%%%% CONTINUOUS EPOCHS %%%%%%
        %%
        
        LabelsRaw = [];
        Fs = 30000;
        
        if is2D == 1
            trialInds = find(newTraining.SignalType ~= 0);
            
            if strcmp(ContinuousLabel2D,"Ball Location")
                LabelsRaw = [newTraining.BallPosX, newTraining.BallPosY];
            elseif strcmp(ContinuousLabel2D,"Target Location")
                LabelsRaw = [newTraining.TargetPosX, newTraining.TargetPosY];
            elseif strcmp(ContinuousLabel2D,"Ball Direction")
                LabelsRaw = [newTraining.BallPosX, newTraining.BallPosY];
            end
            
        else
            trialInds = find(newTraining.Trial > 0.1);
            
            if strcmp(ContinuousLabelVR,"Target Position")
                LabelsRaw = [newTraining.WriTarget, newTraining.ElbTarget];
            elseif strcmp(ContinuousLabelVR,"VR Arm Position")
                LabelsRaw = [newTraining.VRWriPos, newTraining.VRElbPos];
            elseif strcmp(ContinuousLabelVR,"MyoPro Position")
                LabelsRaw = [newTraining.MyoProWriPos, newTraining.MyoProElbPos];
            end
        end
        
        trialEndingInds = find(diff(trialInds) > 25);
        trialIndCell = cell(length(trialEndingInds)+1, 1);
        
        % Build a cell of arrays where each array contains the indicies of each trial
        trialIndCell{1,1} = trialInds(1:trialEndingInds(1));
        for i=2:length(trialEndingInds)
            trialIndCell{i,1} = trialInds(trialEndingInds(i-1)+1:trialEndingInds(i));
        end
        trialIndCell{length(trialEndingInds)+1,1} = trialInds(trialEndingInds(end)+1:end);
        
        ContLabels = cell(length(trialEndingInds)+1, 1);
        % ContVolts = cell(length(trialEndingInds)+1, 1);
        ContSpks = cell(length(trialEndingInds)+1, 1);
        ContTS = cell(length(trialEndingInds)+1, 1);
        
        for i=1:length(trialIndCell)
            
            ContLabels{i,1} = [LabelsRaw(trialIndCell{i,1},1), LabelsRaw(trialIndCell{i,1},2)];
            
            [A,~] = intersect(fTime,newTraining.NSPTime(trialIndCell{i,1}));
            [~,rawtimeInds] = find(fTime == A);
            %   ContVolts{i,1} = rawData(:,rawtimeInds(1):rawtimeInds(end));
            
            r = 30000/Fs;
            idSp = find(events>=(rawtimeInds(1)*r) & events<(rawtimeInds(end)*r));
            ContSpks{i,1} = [double(electrode(idSp));double(events(idSp))];
            ContTS{i,1} = r.*[rawtimeInds(1),rawtimeInds(end)];
        end
        
        %Get trial Inds INCLUDING dropped frames
        AllInclusiveTrialInds = [];
        for i=1:length(trialIndCell)
            AllInclusiveTrialInds = [AllInclusiveTrialInds;(trialIndCell{i,1}(1):trialIndCell{i,1}(end))'];
        end
        % Store indicies associated with rest (everything BUT trial indicies) in array.
        restInds = (1:height(newTraining))';
        restInds(AllInclusiveTrialInds) = [];
        restEndingInds = find(diff(restInds) > 25);
        restIndCell = cell(length(restEndingInds), 1);
        
        % The first start period begins three seconds before the beginning of
        % the first trial
        firstRestStart = find(round(newTraining.NSPTime,1) == round(newTraining.NSPTime(restEndingInds(1)) - 3.0,1), 1);
        
        restIndCell{1,1} = restInds(firstRestStart:restEndingInds(1));
        for i=2:length(restIndCell)
            
            % This checks if the rest time is longer than 4 seconds, indicating
            % the beginning of a new series of trials. In this case it only
            % takes takes the final 3 seconds of rest
            if (newTraining.NSPTime(restInds(restEndingInds(i))) - newTraining.NSPTime(restInds(restEndingInds(i-1)+1)) > 4 ...
                    && is2D == 1)
                firstRestStart = find(round(newTraining.NSPTime,2) == round(newTraining.NSPTime(restEndingInds(i)),2) - 3, 1);
                restIndCell{i,1} = restInds(firstRestStart:restEndingInds(i));
            else
                % Otherwise use the entire rest period
                restIndCell{i,1} = restInds(restEndingInds(i-1)+1:restEndingInds(i));
            end
        end
        
        % RestContVolts = cell(length(restEndingInds), 1);
        RestContSpks = cell(length(restEndingInds), 1);
        RestContTS = cell(length(restEndingInds), 1);
        
        for i=1:length(trialIndCell)
            [A,~] = intersect(fTime,newTraining.NSPTime(restIndCell{i,1}));
            [~,rawtimeInds] = find(fTime == A);
            %  RestContVolts{i,1} = rawData(:,rawtimeInds(1):rawtimeInds(end));
            
            r = 30000/Fs;
            idSp = find(events>=(rawtimeInds(1)*r) & events<(rawtimeInds(end)*r));
            RestContSpks{i,1} = [double(electrode(idSp));double(events(idSp))];
            RestContTS{i,1} = r.*[rawtimeInds(1),rawtimeInds(end)];
        end
        
    end
    % HanEpochVolts and ElbEpochVolts now contain the raw data grouped according to the
    % different trial conditions. These signals can now be used to derive
    % features and susequently feed those features into the machine learning
    % training system
    
    % HanEpochSpks and ElbEpochSpks now contain the spike extracted from the
    % NSP and epoched follwoing the same criterion used for the raw voltages.
    % These spike time stamps need to be broked in specific time bins to
    % extract spike features in next processing step.
    
    %-------------------------------------------------------------------------%
    % Export the Epoched data
    [~,nm,~] = fileparts(filename);
    dest = [pwd,filesep,'Epochs'];
    if(exist(dest,'dir') ~= 7)
        mkdir(dest);
    end
    
    
    
    % Export all data to files to be used in subsequent analysis
    if strcmp(DiscCont,"Discrete")
        %save([dest,filesep,nm,'-HanEpochVolts','.mat'],'HanEpochVolts','-v7.3');
        %save([dest,filesep,nm,'-ElbEpochVolts','.mat'],'ElbEpochVolts','-v7.3');
        save([dest,filesep,nm,'-HanEpochSpks','.mat'],'HanEpochSpks','-v7.3');
        save([dest,filesep,nm,'-ElbEpochSpks','.mat'],'ElbEpochSpks','-v7.3');
        save([dest,filesep,nm,'-HanEpochTS','.mat'],'HanEpochTS','-v7.3');
        save([dest,filesep,nm,'-ElbEpochTS','.mat'],'ElbEpochTS','-v7.3');
    else
        save([dest,filesep,nm,'-ContLabels','.mat'],'ContLabels','-v7.3');
        %save([dest,filesep,nm,'-ContVolts','.mat'],'ContVolts','-v7.3');
        save([dest,filesep,nm,'-ContSpks','.mat'],'ContSpks','-v7.3');
        save([dest,filesep,nm,'-ContTS','.mat'],'ContTS','-v7.3');
        %save([dest,filesep,nm,'-RestContVolts','.mat'],'RestContVolts','-v7.3');
        save([dest,filesep,nm,'-RestContSpks','.mat'],'RestContSpks','-v7.3');
        save([dest,filesep,nm,'-RestContTS','.mat'],'RestContTS','-v7.3');
    end
    
    
    
end


end


%% FEATURE EXTRACTION NEXT


