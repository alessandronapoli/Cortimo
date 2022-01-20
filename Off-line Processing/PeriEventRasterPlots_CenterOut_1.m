% Build raster plots for trial sequences.

clearvars; close all; clc;

% STEP 1
% This will grab the raw Blackrock files and Cortimo training files,
% synchronize them and generate new files that are the epoched version of
% the continuous raw data ready to be used for feature extraction
% addpath(genpath('C:\Users\aless\OneDrive - Thomas Jefferson University and its Affiliates\Jefferson\Cortimo Project\MatlabCode\Off-line Analysis'));
addpath(genpath([pwd,filesep,'NPMK']));

% Stores continuous data
[file,path] = uigetfile('*.*','Select raw data file');
filename = fullfile(path,file);

  %%%%%%%%%%%%%%%%%%%%%%%
    if filename==0
        % user pressed cancel
        return;
    end
    
    
% Open corresponding NEV file to load the events (extracted spikes)
NEV = openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');
comments = NEV.Data.Comments;
events = double(NEV.Data.Spikes.TimeStamp);
electrode = double(NEV.Data.Spikes.Electrode);
Unit = NEV.Data.Spikes.Unit;


  %% Load the training file(s)
    [fileList,pathT] = uigetfile(strcat(path,'*.bin'),'Select file with training data',...
        'MultiSelect','on');
    
    %%%%%%%%%%%%%%%%%%%%%%%
    if (~iscell(fileList) & fileList==0)
        % user pressed cancel
        return;
    end
    
    
     if(iscell(fileList))
        LEN = length(fileList);
        
    else
        LEN = 1;
        tp{1} = fileList;
        fileList = tp;
     end
    
     
     
      % Concatenate all the training data files
    VRTrainingLabels = [];
    EMGMatrix = [];
    
    for file=1:LEN
        
        fileT = fileList{file};
        
        % filenameFeat = fullfile(pathT,fileT);
        filenameTrain = fullfile(pathT,fileT);
        
        % Open both files, the training data file and the EEG feature file
        % filenameTrain = strrep(filenameFeat,"FreqFeatures","ArmVR");
        try
            if contains(filenameTrain, "CenterOut")
                [VRTrainingLabels] = [VRTrainingLabels;loadCenterOutTrainingFile(filenameTrain)];
                is2D = 1;
            else
                [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile(filenameTrain)];
                is2D = 0;
                
                fileE = strrep(fileT,'ArmVR','MyoEMG');
                if(isfile([pathT,filesep,fileE]))
                    s = dir([pathT,filesep,fileE]);
                    if(s.bytes>0)
                       [EMGMatrix] = [EMGMatrix;load_MyoEMG([pathT,filesep,fileE])];
                       isEMG = true;
                       EMGMatrix.NSPTime = round(EMGMatrix.NSPTime,5);

                    else
                        isEMG = false;
                    end
                end
                
            end
            % [EEGfeats] = loadEEGFeatureFile(filenameFeat);
        catch
            [VRTrainingLabels] = [VRTrainingLabels;loadVrTrainingFile_noDecoderOutput(filenameTrain)];
            is2D = 0;
            
        end
        
    end
    
    
    
    dd = find(diff(VRTrainingLabels.MatlabCycle)~=0);
    newTraining = VRTrainingLabels(dd,:);
    newTraining.NSPTime = round(newTraining.NSPTime,5);
    
    %%%%%% SPIKES   %%%%%%
% Spike parameters
spkChs = 30:100;
spkTWs = [0.05];

OnlyUseTheseHandLabels = [-1,1,2];
OnlyUseTheseElbowLabels = [-1,1,2];

Fs = 30e3;


%% Re-Organize the spikes for raster plot
    
TrainingTimeStartID = find(events./Fs >= newTraining.NSPTime(1),1);
TrainingTimeEndID = find(events./Fs > newTraining.NSPTime(length(newTraining.NSPTime)),1);

timeVector = (events(TrainingTimeStartID):events(TrainingTimeEndID))./Fs;
step = 200*30e3;
chunks = 1:step:length(timeVector);
bin = 1/Fs;

ff1 = figure(1);

storage = cell(1,length(chunks));
spkbinning = 0.1; % in sec
left = 0;

for cc=1:length(chunks)
    
    if(cc<length(chunks))
    thistimeVector = timeVector(chunks(cc):chunks(cc+1)-1);
    else
    thistimeVector = timeVector(chunks(cc):length(timeVector));
    end

spikeMatrix = zeros(length(thistimeVector),length(spkChs));


    for ch=spkChs(1):spkChs(end)
    
        ChId = find(electrode == ch);
        
        [C,ia,ib] = intersect(events(ChId)./Fs,thistimeVector);
        spikeMatrix(ib,ch) = 1;
        
    end
    
    toPlot = spikeMatrix;
    % toPlot(isnan(toPlot))=0;
    spikes = logical(toPlot)';
    [~,~] = plotSpikeRaster(spikes,'PlotType','vertline',...
    'TimePerBin',bin,'rasterWindowOffset',thistimeVector(1));

    
    %idd = thistimeVector(1):spkbinning:thistimeVector(end); 
    idd = 1:spkbinning*Fs:size(spikeMatrix,1);
    spC = zeros(length(idd),size(spikeMatrix,2));
    
    for w=1:length(idd)
        
        if(w<length(idd))
            if(w==1)
              spC(w,:) = sum(spikeMatrix(idd(w):idd(w+1),:)) + left;
              left = 0;
            else

               spC(w,:) = sum(spikeMatrix(idd(w):idd(w+1),:));
            end
        else
            left = sum(spikeMatrix(idd(w):end,:));
        end
        
    end
    
    storage{1,cc} = spC;
    storage{2,cc} = thistimeVector(1) + idd ./Fs;
end



% f1 = figure(1);
%  %Find its handle
% aH = gca;
% fH = ancestor(aH,'fig');
% fH(2) = figure(2); %Figure you want to copy the stuff to
% %Copy axes or see my previous code to copy lines, store its handle
% aH(2) = copyobj(aH,fH(2));

%%
figure(1);
hold on
% plot(newTraining.NSPTime-newTraining.NSPTime(1),newTraining.MyoProWriPos,...
%      'LineWidth',3,'color','magenta');
 plot(newTraining.NSPTime,100.*newTraining.BallPosX ./ norm(newTraining.BallPosX),...
     'LineWidth',3,'color','magenta');
 plot(newTraining.NSPTime,100.*newTraining.BallPosY ./ norm(newTraining.BallPosY),...
     'LineWidth',3,'color','cyan');

 view(0,-90);
 hold off;
title('Norm Cursor Position [x,y]  |||   spike counts');
xlabel('Time [s]');
ylabel('Channels    |||    Cursor Positions [%]');
%legend('','MyoPro Pos');


%% Find correlation between signals and spike count
spC = [];
TT = [];

for c=1:size(storage,2)
    
    spC = [spC; storage{1,c}];
    TT = [TT;storage{2,c}'];
    
end


% Center Out Cursor X Position
tssig = timeseries(newTraining.BallPosX,newTraining.NSPTime);
%Ttime = newTraining.NSPTime;
tsSync = resample(tssig,TT);
CursorPosXSync = tsSync.data;


% Center Out Cursor Y Position
tssig = timeseries(newTraining.BallPosY,newTraining.NSPTime);
tsSync = resample(tssig,TT);
CursorPosYSync = tsSync.data;


CMatCursorPos = [];

for ch=1:size(spC,2)
   
   % Hand 1)Extension 2) Flexion
    [R,P] = corrcoef(spC(:,ch),CursorPosXSync,'Rows','Complete');
    CMatCursorPos(1,ch) = R(2);
    [R,P] = corrcoef(spC(:,ch),CursorPosYSync,'Rows','Complete');
    CMatCursorPos(2,ch) = R(2);
    
end


[RCurX,BestCurX] = sort(CMatCursorPos(1,:),'descend');
nna = find(~isnan(RCurX)); 
RCurX = RCurX(nna);
BestCurX = BestCurX(nna);

[RCurY,BestCurY] = sort(CMatCursorPos(2,:),'descend');
nna = find(~isnan(RCurY)); 
RCurY = RCurY(nna);
BestCurY = BestCurY(nna);


%% Select a few correlated channels

% Best Cursor X
nCh = 5;
sigs = [];
ti = [];
for c=1:size(storage,2)

sigs = [sigs;storage{1,c}(:,BestCurX(1:nCh))];
ti = [ti,storage{2,c}];
end
figure;
hold on;
plot(ti,sigs,'LineWidth',2);
plot(TT,100.*CursorPosXSync./norm(CursorPosXSync),'Color','Yellow','LineWidth',3,'LineStyle','--');
tt = cellstr(num2str([BestCurX(1:nCh)',round(RCurX(1:nCh),3)']));
tt{nCh+1} = 'Cursor X Position';
legend(tt);
title('Spike Count in time bins with Cursor data');
ylabel('Spike Count  ||  Norm Cursor Position [%]');
xlabel('Time [s]')
hold off;


% Best Cursor Y
nCh = 5;
sigs = [];
ti = [];

for c=1:size(storage,2)

sigs = [sigs;storage{1,c}(:,BestCurY(1:nCh))];
ti = [ti,storage{2,c}];
end
figure;
hold on;
plot(ti,sigs,'LineWidth',2);
plot(TT,100.*CursorPosYSync./norm(CursorPosYSync),'Color','Green','LineWidth',3,'LineStyle','--');
tt = cellstr(num2str([BestCurY(1:nCh)',round(RCurY(1:nCh),3)']));
tt{nCh+1} = 'Cursor Y Position';
legend(tt);
title('Spike Count in time bins with Cursor data');
ylabel('Spike Count  ||  Norm Cursor Position [%]');
xlabel('Time [s]')
hold off;