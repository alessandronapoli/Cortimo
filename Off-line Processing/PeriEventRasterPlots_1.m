

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
    
%     EEMGMatrix = EMGMatrix(dd,:);

DiscCont = "Discrete"; % Set to Discrete for discrete classifiers

% Parameters for creating new training labels from position data 
% "Target Position", "VR Arm Position", "MyoPro Position"
ContinuousLabelVR =  "VR Arm Position";
% ContinuousLabelVR =  "Target Position";

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

%% Add the training label information

%%%%%%%%%%%%%%% This works for repeated trials when VR Labels are used %%%%%%%%%%%%%  
idd = find(diff(newTraining.Trial)~=0) + 1;
if(~isempty(idd))
TrialON = sign(diff([0;newTraining.Trial]));
typeofTrial = [newTraining.WriLab(idd+1),newTraining.ElbLab(idd+1)];
condLabelsHand = {'Hand Hold';'Hand Flex'; 'Hand Extend'};
condLabelsElbow = {'Elbow Hold';'Elbow Flex'; 'Elbow Extend'};

for i=1:length(idd)
if(TrialON(idd(i)) == 1)
    
%     vline(newTraining.NSPTime(idd(i))-newTraining.NSPTime(1),'g',...
%     ['START ',condLabelsHand{typeofTrial(i)},condLabelsElbow{typeofTrial(i)}]);
% vline(newTraining.NSPTime(idd(i))-newTraining.NSPTime(1),'g',[]);
vline(newTraining.NSPTime(idd(i)),'g',[]);

%     text(newTraining.NSPTime(idd(i))-newTraining.NSPTime(1),2,{'START ';condLabelsHand{typeofTrial(i)+1},...
%         ;condLabelsElbow{typeofTrial(i)+1}},'color','green');
text(newTraining.NSPTime(idd(i)),-2,{'START ';condLabelsHand{typeofTrial(i)+1},...
        ;condLabelsElbow{typeofTrial(i)+1}},'color','green');
    
else
    if(TrialON(idd(i)) == -1)
%     vline(newTraining.NSPTime(idd(i))-newTraining.NSPTime(1),'r','');
%      text(newTraining.NSPTime(idd(i))-newTraining.NSPTime(1),-2,'STOP');
     vline(newTraining.NSPTime(idd(i)),'r','');
     text(newTraining.NSPTime(idd(i)),-2,'STOP','Color','r');
    end
    
end
end


end



%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% THIS IS TO BE USED WHEN VR TARGET IS USED %%%%%%%%%%%%%
% 
% if(isempty(idd))
%     
%  
% 
% f1 = figure(1);
%  %Find its handle
% aH = gca;
% fH = ancestor(aH,'fig');
% fH(2) = figure(2); %Figure you want to copy the stuff to
% %Copy axes or see my previous code to copy lines, store its handle
% aH(2) = copyobj(aH,fH(2));
% 
% figure(1);
% hold on
% plot(newTraining.NSPTime-newTraining.NSPTime(1),newTraining.WriTarget,...
%      'LineWidth',3);
%  view(0,-90);
%  hold off;
% title('Hand Position   |||   spike counts');
% xlabel('Time [s]');
% ylabel('Channels    |||    Degrees of flexion');
% 
% figure(2);
% hold on
% plot(newTraining.NSPTime-newTraining.NSPTime(1),newTraining.ElbTarget,...
%      'LineWidth',3);
% view(0,-90);
% hold off;
% title('Elbow Position   |||   spike counts');
% xlabel('Time [s]');
% ylabel('Channels   |||   Degrees of flexion');
% 
%     
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% THIS IS TO BE USED WHEN MYOPRO POSITION IS USED %%%%%%%%%%%%%

if(isempty(idd) && isEMG == 0)

f1 = figure(1);
 %Find its handle
aH = gca;
fH = ancestor(aH,'fig');
fH(2) = figure(2); %Figure you want to copy the stuff to
%Copy axes or see my previous code to copy lines, store its handle
aH(2) = copyobj(aH,fH(2));

figure(1);
hold on
% plot(newTraining.NSPTime-newTraining.NSPTime(1),newTraining.MyoProWriPos,...
%      'LineWidth',3,'color','magenta');
 plot(newTraining.NSPTime,newTraining.MyoProWriPos,...
     'LineWidth',3,'color','magenta');

 view(0,-90);
 hold off;
title('Hand Position   |||   spike counts');
xlabel('Time [s]');
ylabel('Channels    |||    Degrees of flexion');
%legend('','MyoPro Pos');

figure(2);
hold on
% plot(newTraining.NSPTime-newTraining.NSPTime(1),newTraining.MyoProElbPos,...
%      'LineWidth',3,'color','magenta');
plot(newTraining.NSPTime,newTraining.MyoProElbPos,...
     'LineWidth',3,'color','magenta');

view(0,-90);
hold off;
title('Elbow Position   |||   spike counts');
xlabel('Time [s]');
ylabel('Channels   |||   Degrees of flexion');
%legend('','MyoPro Pos')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% THIS IS TO BE USED WHEN MYOPRO EMG SIGS ARE USED  %%%%%%%%%%%%%
%%
if(isEMG == true)
    
f1 = figure(1);
 %Find its handle
aH = gca;
fH = ancestor(aH,'fig');
fH(2) = figure(2); %Figure you want to copy the stuff to
%Copy axes or see my previous code to copy lines, store its handle
aH(2) = copyobj(aH,fH(2));

figure(1);
    hold on;
%     plot(EMGMatrix.NSPTime-newTraining.NSPTime(1),EMGMatrix.HandExtend,...
%         'linewidth',3,'color','yellow');
%     plot(EMGMatrix.NSPTime-newTraining.NSPTime(1),EMGMatrix.HandFlex,...
%         'linewidth',3,'color','green');
    plot(newTraining.NSPTime,newTraining.MyoProWriPos,...
        'LineWidth',3,'color','magenta');    
    plot(EMGMatrix.NSPTime,EMGMatrix.HandExtend,...
        'linewidth',3,'color','yellow');
    plot(EMGMatrix.NSPTime,EMGMatrix.HandFlex,...
        'linewidth',3,'color','green');
    plot(newTraining.NSPTime,newTraining.MyoProElbPos,...
        'LineWidth',3,'Color','blue','LineStyle','-.');
    
 
    hold off;
    legend('','MyoPro Pos','Ext','Flex','MyoPro Elb')
    title('Hand Position   |||   spike counts   |||   EMG');
    xlabel('Time [s]');
    ylabel('Channels   |||   Degrees of flexion  |||   EMG');
    view(0,-90);

    
%     %Find its handle
%     aH = gca;
%     fH = ancestor(aH,'fig');
%     fH(3) = figure(3); %Figure you want to copy the stuff to
%     %Copy axes or see my previous code to copy lines, store its handle
%     aH(3) = copyobj(aH,fH(3));
    figure(2);
    hold on;
%     plot(EMGMatrix.NSPTime-newTraining.NSPTime(1),EMGMatrix.ElbowExtend,...
%         'linewidth',3,'color','red');
%     plot(EMGMatrix.NSPTime-newTraining.NSPTime(1),EMGMatrix.ElbowFlex,...
%         'linewidth',3,'color','blue');
    plot(newTraining.NSPTime,newTraining.MyoProElbPos,...
        'LineWidth',3,'color','magenta');
    plot(EMGMatrix.NSPTime,EMGMatrix.ElbowExtend,...
        'linewidth',3,'color','red');
    plot(EMGMatrix.NSPTime,EMGMatrix.ElbowFlex,...
        'linewidth',3,'color','blue');
    hold off;
    legend('','MyoPro Pos','Triceps','Biceps')
    title('Elbow Position   |||   spike counts   |||   EMG');
    xlabel('Time [s]');
    ylabel('Channels   |||   Degrees of flexion  |||   EMG');
    view(0,-90);

end


%% Plot spike count
% figure(33)
% hold on;
% for c=1:length(storage)
% bar(storage{2,c},storage{1,c}(:,33),'EdgeColor','blue','FaceColor','blue')
% end
% hold off;


%% Find correlation between signals and spike count
spC = [];
TT = [];

for c=1:length(storage)
    
    spC = [spC; storage{1,c}];
    TT = [TT;storage{2,c}'];
    
end

% CMat = [];
% for ch=1:size(spC,2)
%    
%     CMat(:,ch) = xcorr(spC(:,ch),newTraining.MyoProWriPos);
% end

% MyoPro Hand position
tssig = timeseries(newTraining.MyoProWriPos,newTraining.NSPTime);
%Ttime = newTraining.NSPTime;
tsSync = resample(tssig,TT);
MyoProWriPosSync = tsSync.data;


% MyoPro Elbow position
tssig = timeseries(newTraining.MyoProElbPos,newTraining.NSPTime);
tsSync = resample(tssig,TT);
MyoProElbPosSync = tsSync.data;

if(isEMG == 1)
% EMG Hand Extension
tssig = timeseries(EMGMatrix.HandExtend,EMGMatrix.NSPTime);
tsSync = resample(tssig,TT);
EMGHandExtSync = tsSync.data;


% EMG Hand Flexion
tssig = timeseries(EMGMatrix.HandFlex,EMGMatrix.NSPTime);
tsSync = resample(tssig,TT);
EMGHandFlexSync = tsSync.data;


% EMG Elbow Extension
tssig = timeseries(EMGMatrix.ElbowExtend,EMGMatrix.NSPTime);
tsSync = resample(tssig,TT);
EMGElbowExtSync = tsSync.data;


% EMG Elbow Flexion
tssig = timeseries(EMGMatrix.ElbowFlex,EMGMatrix.NSPTime);
tsSync = resample(tssig,TT);
EMGElbowFlexSync = tsSync.data;

else
   EMGHandExtSync = [];
   EMGHandFlexSync = [];
   EMGElbowExtSync = [];
   EMGElbowFlexSync = [];
end

CMatWriPos = [];
CMatElbPos = [];
CMatHandEMG = [];
CMatElbEMG = [];


for ch=1:size(spC,2)
   
    % CMat(:,ch) = xcorr(spC(:,ch),MyoProWriPosSync);
    [R,P] = corrcoef(spC(:,ch),MyoProWriPosSync,'Rows','Complete'); 
    CMatWriPos(:,ch) = R(2);
    [R,P] = corrcoef(spC(:,ch),MyoProElbPosSync,'Rows','Complete'); 
    CMatElbPos(:,ch) = R(2);
    
    if(isEMG == 1)
    % Hand 1)Extension 2) Flexion
    [R,P] = corrcoef(spC(:,ch),EMGHandExtSync,'Rows','Complete');
    CMatHandEMG(1,ch) = R(2);
    [R,P] = corrcoef(spC(:,ch),EMGHandFlexSync,'Rows','Complete');
    CMatHandEMG(2,ch) = R(2);
    
    % Elbow 1)Extension 2) Flexion
    [R,P] = corrcoef(spC(:,ch),EMGElbowExtSync,'Rows','Complete');
    CMatElbEMG(1,ch) = R(2);
    [R,P] = corrcoef(spC(:,ch),EMGElbowFlexSync,'Rows','Complete');
    CMatElbEMG(2,ch) = R(2);
    end

end

% The follwoing variables will give you the most correlated channels
% relative to several control signals 
[RHanPos,BestHanPos] = sort(CMatWriPos,'descend');
nna = find(~isnan(RHanPos)); 
RHanPos = RHanPos(nna);
BestHanPos = BestHanPos(nna);

[RElbPos,BestElbPos] = sort(CMatElbPos,'descend');
nna = find(~isnan(RElbPos)); 
RElbPos = RElbPos(nna);
BestElbPos = BestElbPos(nna);

if(isEMG == 1)
[RHanExt,BestHanExt] = sort(CMatHandEMG(1,:),'descend');
nna = find(~isnan(RHanExt)); 
RHanExt = RHanExt(nna);
BestHanExt = BestHanExt(nna);

[RHanFlex,BestHanFlex] = sort(CMatHandEMG(2,:),'descend');
nna = find(~isnan(RHanFlex)); 
RHanFlex = RHanFlex(nna);
BestHanFlex = BestHanFlex(nna);


[RElbExt,BestElbExt] = sort(CMatElbEMG(1,:),'descend');
nna = find(~isnan(RElbExt)); 
RElbExt = RElbExt(nna);
BestElbExt = BestElbExt(nna);

[RElbFlex,BestElbFlex] = sort(CMatElbEMG(2,:),'descend');
nna = find(~isnan(RElbFlex)); 
RElbFlex = RElbFlex(nna);
BestElbFlex = BestElbFlex(nna);

end



%% Select a few correlated channels


if(isEMG == 1)
% Best Hand Ext
nCh = 5;
sigs = [];
ti = [];
for c=1:length(storage)

sigs = [sigs;storage{1,c}(:,BestHanExt(1:nCh))];
ti = [ti,storage{2,c}];
end
figure;
hold on;
plot(ti,sigs,'LineWidth',2);
plot(TT,EMGHandExtSync,'Color','Yellow','LineWidth',3,'LineStyle','--');
tt = cellstr(num2str([BestHanExt(1:nCh)',round(RHanExt(1:nCh),3)']));
tt{nCh+1} = 'Hand Extensors';
legend(tt);
title('Spike Count in time bins with Hand Extension data');
ylabel('Spike Count  ||  EMG');
xlabel('Time [s]')
hold off;

% Best Hand Flex
nCh = 5;
sigs = [];
ti = [];

for c=1:length(storage)

sigs = [sigs;storage{1,c}(:,BestHanFlex(1:nCh))];
ti = [ti,storage{2,c}];
end
figure;
hold on;
plot(ti,sigs,'LineWidth',2);
plot(TT,EMGHandFlexSync,'Color','Green','LineWidth',3,'LineStyle','--');
tt = cellstr(num2str([BestHanExt(1:nCh)',round(RHanFlex(1:nCh),3)']));
tt{nCh+1} = 'Hand Flexors';
legend(tt);
title('Spike Count in time bins with Hand Flexion data');
ylabel('Spike Count  ||  EMG');
xlabel('Time [s]')
hold off;

end


% Best Hand Pos
nCh = 5;
sigs = [];
ti = [];

for c=1:length(storage)

sigs = [sigs;storage{1,c}(:,BestHanPos(1:nCh))];
ti = [ti,storage{2,c}];
end
figure;
hold on;
plot(ti,sigs,'LineWidth',2);
plot(TT,MyoProWriPosSync,'Color','Magenta','LineWidth',3,'LineStyle','--');
tt = cellstr(num2str([BestHanPos(1:nCh)',round(RHanPos(1:nCh),3)']));
tt{nCh+1} = 'Hand Pos';
tt{nCh+2} = 'Elbow Pos';
plot(TT,MyoProElbPosSync,'Color','Blue','LineWidth',2,'LineStyle','-.');
legend(tt);
title('Spike Count in time bins with Hand data');
ylabel('Spike Count  ||  Degrees of Flexion');
xlabel('Time [s]')
hold off;



