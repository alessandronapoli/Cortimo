%% Open NEV file
% % Stores spike count
% openNEV('report','read','uV','nomat','nosave');
% 
% events = NEV.Data.Spikes.TimeStamp;
% electrode = NEV.Data.Spikes.Electrode;
% Unit = NEV.Data.Spikes.Unit;
clearvars; close all; clc;
addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('AuxFunc'));
% rmpath(['eeglab',filesep,'plugins']);
%% Open corresponding NS file
% Stores continuous data
[file,path] = uigetfile('*.*','Select raw data file');
filename = fullfile(path,file);
%openNSx('report','read',filename,'sample', 'p:short', 's:75');
% openNSx('report','read',filename,'sample', 'p:short');
openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');
comments = NEV.Data.Comments;

openNSx('report','read','uV',filename,'t:0:40','min','p:short');
% openNSx('report','read','uV',filename,'p:short','s:2');

% openNSx('report','read','uV',filename,'p:short');

% Depending on the data format the structure names are different
forma = filename(end);
%rawData = NS5.Data;
tp = eval(['NS',forma]);
rawData = tp.Data;
rawtime = (0:length(rawData)-1)*(1/tp.MetaTags.SamplingFreq); % in seconds
Fs = double(tp.MetaTags.SamplingFreq);

% %% Load the training file
% [fileT,pathT] = uigetfile('*.bin','Select file with training data');
% filenameT = fullfile(pathT,fileT);
% 
% % training_data = fread(fid,inf,'double');
% training_data = loadVrTrainingFile(filenameT);


% 
% stream_labels = {'AcqCycle','TimeStamp','Trial','WriTarget',...
%     'ElbTarget','WriLab','ElbLab',...
%     'WriHit','ElbHit','VRWriPos','VRElbPos','VRMyomoWriPos','VRMyomoElbPos',...
%     'MyoProWriPos','MyoProElbPos'};
% 

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
% traningStartID = find(training_data{1,3} == rawtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data have been loaded. The following analysis is composed of severeal
% sections.
% Section 1: Find whether there are any difference between electrodes on
% the skull defect VS. the homologous contralateral electrodes

%----------------------- SECTION 1 -------------------------------------- %
% Analyze data time windows from homologous electrodes
ChLabels = {'FP1','FP2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T3','T4','T5','T6','Z-Gr','Fz','Cz','Pz','A1','A2','C2',...
    'C6','FC4','CP4','CP2','CP6','FC2','FC6'};
Ch = 1:30;

% Group the elctrodes in hemicraniectomy and skull
hemi = {'F4','C4','P4','T4'};
skull = {'F3','C3','P3','T3'};

additionalHemi = {'C2','C6','FC4','CP4','CP2','CP6','FC2','FC6'};


%% Break down the data into 500ms time windows, then average and display
% power spectra
Timedur = 1200;
% st = 700 * Fs; % Use this for file with baseline
st = 600 * Fs; % Use this for file with training trials
sig = rawData(:,st:st-1+Timedur*Fs)';

% Filter
Hd1 = Butt_bandpass(0.5,1000,Fs);
sig = filter(Hd1,sig);
% Downsample to 2000 sps
Fsdown = 2000;
factor = Fs/Fsdown;
sig = downsample(sig,factor);


%%%%%%%%%%%%%%%%%%%%%%%%% Notch Filter  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Notch filter for 60Hz
 d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fsdown);
 sig = filtfilt(d,sig); 
 
 % detrend the data
 sig = detrend(sig);


%%%%%%%%%%%%%%% USE CHANNEL ZERO TO RE-FERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%% RE-REFERENCE %%%%%%%%%%%%%%%%%%
% Signals were recorded with A1 as only reference
% Let's use the average between A1 and A2
[~,~,a2] = intersect({'A1','A2'},ChLabels);
sig = sig - mean(sig(:,a2),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% CAR filter %%%%%%%%%%%%%%%%%%
% Use Common Average Reference
% sig = sig - mean(sig,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Quick plot
% figure
% ll = 60*Fsdown;
% plot((1:ll)/Fsdown,epoched(1:ll,7))
% hold on
% plot((1:ll)/Fsdown,epoched(1:ll,8))
% title('Homologous Electrode Amplitudes')
% xlabel('uV')
% ylabel('Time [s]')
% legend('P3','P4');
% 
%%
% Get spectra of time windows
w = 0.5 * Fsdown;
% Select overlap between windows
overlap = 0 ;

% count = 1;
stepby = w*(1-overlap);
epoched = zeros(w,size(sig,2),round(size(sig,1)/stepby)-1);
here = w;
stepin = w*overlap;
for window=1:round((size(sig,1)/stepby))-1
    if(window == 1)
        epoched(:,:,window) = sig(1:w,:);
    else
        ti = here-stepin+1; 
        epoched(:,:,window) = sig(ti:ti+w-1,:);
        here = ti+w-1;      
        
    end
end

% for i=1:stepby:size(sig,1)
% epoched(:,:,count) = sig(i:i+w-1,:);
% count = count + 1;
% end

Fstep = 2;
Fr = 0:Fstep:Fsdown/2-1/Fsdown/2;
MeanSpectra = zeros(length(Fr),size(sig,2),size(epoched,3));

for trial=1:size(epoched,3)
    
    MeanSpectra(:,:,trial) = periodogram(epoched(:,:,trial),[],Fr,Fsdown);

end

% Average the Power Spectral Density
MeanSpectra = mean(MeanSpectra,3);
figure
surf(1:size(MeanSpectra,2),Fr,10*log10(MeanSpectra))
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(Ch);
xticklabels(ChLabels);
colormap(jet);
caxis([-20 30]); % Change limits
axis('tight')
% ylim([0,500])


%% Look for homologous electrode differences

% idHemi = find(strcmp(ChLabels,hemi));
[~,~,idHemi] = intersect(hemi,ChLabels);
[~,~,idSkull] = intersect(skull,ChLabels);
idHemi = sort(idHemi);
idSkull = sort(idSkull);

% idSkull = find(strcmp(ChLabels,skull)); 
[~,~,addHemi] = intersect(additionalHemi,ChLabels);

% Plot the means of the two groups
figure
% x = mean(MeanSpectra(:,idHemi),2);
x = mean(MeanSpectra(:,[idHemi;addHemi]),2);
h = plot(Fr,10*log10(x));
plotFonts(h);
hold on
x = mean(MeanSpectra(:,idSkull),2);
h = plot(Fr,10*log10(x));
title('Mean Power Spectra')
legend('Hemicraniectomy','Skull')
ylabel('PSD dB/Hz')
xlabel('Hz')
plotFonts(h);


%% Calculate the Correlation between electrodes according tho their distance

% Calculate the Pearson correlation coefficients between electrodes
[R,pvals] = corrcoef(sig);
% Try to quantify how "connected" an electrode is in terms of how
% correlated it is to the other electrodes on average
vec = zeros(size(sig,2),1);
for cha=1:size(sig,2)
    thisv = R(cha,:);
    thisv(cha) = [];
    vec(cha) = mean(thisv);
end

% Topographical plot of the average one-to-all correlation coefficients
% between electrodes. It shows on average how correlated each electrode is
% compared to all the others.
figure
handle = topoplot(vec,'DHC-locs.ced','electrodes','on',...
    'electrodes','labels','shrink','on','gridscale',200,...
    'plotchans',[1:16,18:20,23:30]);
caxis([0,1]);
title('Mean Pearson Correlation Coeffiecients')
plotFonts(handle);

t = (0:length(sig)-1)/Fsdown;
figure
chtoPlot = find(strcmp(ChLabels,'C6'));
plot(t,sig(:,chtoPlot));

% Get the RMS values for all channels to quantify level
rm = rms(sig);
figure
handle = topoplot(rm,'DHC-locs.ced','electrodes','on',...
    'electrodes','labels','shrink','on','gridscale',200,...
    'plotchans',[1:16,18:20,23:30]);
% caxis([0,1]);
title('RMS values')
plotFonts(handle);

% rmpath(genpath('eeglab'));


% % Average trials
% MeanWin = mean(epoched,3);
% 
% % Calculate Spectra
% Fstep = 3;
% Fr = 0:Fstep:Fsdown/2-1/Fsdown/2;
% MeanSpectra = periodogram(MeanWin,[],Fr,Fsdown);
% figure
% surf(1:size(MeanSpectra,2),Fr,10*log10(MeanSpectra))
% view(90,-90)
% ylabel('Frequency Hz')
% xlabel('Channels')
% xticks(Ch);
% xticklabels(ChLabels);
% colormap(jet);
% caxis([-20 30]); % Change limits
% axis('tight')
% ylim([0,500])







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Perform here the trial analysis
% 
% % id = 2:6:length(training_data);
% % take advantage of the data labels
% %events = find(strcmp(stream_labels,'TimeStamp'));
% 
% % id = ts:length(stream_labels):length(training_data);
% % training_time = training_data(id);
% % freq = diff(training_time)*1000;
% 
% events = training_data{:,'NSPTime'};
% 
% condnum = 3;
% TimeEvents = cell(1,condnum);
% epoch = cell(1,condnum);
% for cond=1:condnum
%     matri = [];
%     idd = find(training_data{:,'WriLab'} == 0);
%     win = find(diff(idd)>1);
%     win = [1;win];
%     for i=2:length(win)
%         [~,idRaw,~]= intersect(rawtime,training_data{win,'WriLab'});
%         
%     end
%     [~,~,TimeEvents{cond}] = intersect(training_data{idd,'NSPTime'},rawtime);
%     
% end
% 
% % Find the flexion/extension "labels" that can be used to train a
% % classifier. Based on the 
% ts = find(strcmp(stream_labels,'WriLab'));
% id = ts:length(stream_labels):length(training_data);
% training_labels = [training_data(id),training_data(id+1)];
% 
% % Add the motion data from the VR and MyoPro to the training data
% ts = find(strcmp(stream_labels,'VRWriPos'));
% id = ts:length(stream_labels):length(training_data);
% VRpos = [training_data(id),training_data(id+1)]; % Wri Elb
% ts = find(strcmp(stream_labels,'MyoProWriPos'));
% id = ts:length(stream_labels):length(training_data);
% Myopos = [training_data(id),training_data(id+1)]; % Wri Elb
% 
% 
% 
% % Synchronize raw continuous acquisition and training labels
% start_training = training_time(1); % First training sample in training data
% idraw = find(rawtime==start_training); % First training sample in raw data
% idspk = find(events>start_training*(Fs),1); % First training sample in event data (spikes) 
% % It is also necessary to find the end of the training information
% stop_training = training_time(end); % Last training sample in training data
% last_idraw = find(rawtime==stop_training); % Last training sample in raw data
% last_idspk = find(events<=stop_training*(Fs) & ...
%     events>start_training*(Fs)); % Last training sample in event data
% last_idspk = last_idspk(end);
% 
% % Raw training data
% training_raw = rawData(:,idraw:last_idraw);
% % Event data to synch
% training_events = events(idspk:last_idspk);
% training_electrode = electrode(idspk:last_idspk);
% training_unit = Unit(idspk:last_idspk);
% 
% 
% %%  Now the raw continuos data and the events have been synchronized
% %   Implement signal processing that is similar to the online BCI feature
% %   extraction
% 
% feats = myoffline_processing_2(Fs,3000,training_raw,training_events,...
%     training_time,training_electrode,training_unit,training_labels,...
%     VRpos,Myopos);
% 
