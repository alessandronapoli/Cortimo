

% Build raster plots for trial sequences.

clearvars; close all; clc;

% STEP 1
% This will grab the raw Blackrock files and Cortimo training files,
% synchronize them and generate new files that are the epoched version of
% the continuous raw data ready to be used for feature extraction
%addpath(genpath('C:\Users\aless\OneDrive - Thomas Jefferson University and its Affiliates\Jefferson\Cortimo Project\MatlabCode\Off-line Analysis'));
addpath(genpath([pwd,filesep,'NPMK']));


% Stores continuous data
[file,path] = uigetfile('*.*','Select raw data file');
filename = fullfile(path,file);

  %%%%%%%%%%%%%%%%%%%%%%%
    if (filename == 0)
        % user pressed cancel
        return;
    end
    
    
% Open corresponding NEV file to load the events (extracted spikes)
NEV = openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');
comments = NEV.Data.Comments;
events = double(NEV.Data.Spikes.TimeStamp);
electrode = double(NEV.Data.Spikes.Electrode);
Unit = NEV.Data.Spikes.Unit;

Fs = 30e3;

%% Re-Organize the spikes for raster plot
    
TrainingTimeStartID = find(events >= comments.TimeStamp(1),1) - 100;
TrainingTimeEndID = find(events >= comments.TimeStamp(length(comments.TimeStamp)),1);

timeVector = (events(TrainingTimeStartID):events(TrainingTimeEndID))./Fs;
step = 200*30e3;
chunks = 1:step:length(timeVector);

bin = 1/Fs; % Spike bin in seconds

spkChs = 1:120;

ff1 = figure(1);
storage = cell(1,length(chunks));
spkbinning = 0.05; % in sec
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


%% Add the manual annotations to the plot
hold on;

for com=1:length(comments.TimeStampSec)

    vline(comments.TimeStampSec(com),'r','');
    text(comments.TimeStampSec(com),-2,comments.Text(com,:),...
        'color','r');

end
hold off;



%% Find time stamps for specific events and look for cross-correlations
coms = cellstr(comments.Text);
Estring = find(contains(coms,'elbow'));
Estring2 = find(contains(coms(Estring),'fle'));
% Grab the identified strings
theseStrings = Estring(Estring2);

spC = [];
TT = [];

for c=1:length(storage)
    
    spC = [spC; storage{1,c}];
    TT = [TT;storage{2,c}'];
    
end

select = coms(theseStrings);
selectTimes = comments.TimeStampSec(theseStrings);

CorThis = [];
for tt=1:length(TT)
    if(tt<length(TT))
        idd = find(selectTimes<=TT(tt+1) & selectTimes>TT(tt));
        if(~isempty(idd))
        CorThis(tt) = 100;
        else
        CorThis(tt) = 0;
        end
    else
        idd = find(selectTimes);
         if(~isempty(idd))
        CorThis(tt) = 100;
        else
        CorThis(tt) = 0;
        end
    end
    
 
end


CMat = [];
for ch=1:size(spC,2)
 % CMat(:,ch) = xcorr(spC(:,ch),MyoProWriPosSync);
    [R,P] = corrcoef(spC(:,ch),CorThis,'Rows','Complete'); 
    CMat(:,ch) = R(2);
end


