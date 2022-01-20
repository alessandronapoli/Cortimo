
%close all; clearvars; clc;
addpath('Aux Funcs');
addpath('NPMK');

%% Open corresponding NS file
% Stores continuous data
[file,path] = uigetfile('*.*','Select raw data file');
filename = fullfile(path,file);
%openNSx('report','read',filename,'sample', 'p:short', 's:75');
openNSx('report','read',filename,'sample', 'p:short');

% Depending on the data format the structure names are different
forma = filename(end);
%rawData = NS5.Data;
tp = eval(['NS',forma]);

% Raw data matrix
rawData = tp.Data';
% Time vector corresponding to raw file data
rawtime = (0:length(rawData)-1)*(1/tp.MetaTags.SamplingFreq); % in seconds

% Load corresponding spike file with annotations
% Open corresponding NEV file to load the events (extracted spikes)
NEV = openNEV('report',[filename(1:end-3), 'nev'],'nomat','nosave');

comments = NEV.Data.Comments; % This variable stores user annotations
events = double(NEV.Data.Spikes.TimeStamp); % This variable stores spike time stamps expressed in number of samples
electrode = double(NEV.Data.Spikes.Electrode); % This variable contains electrode number corresponding to the time stamp
Unit = NEV.Data.Spikes.Unit; % This is for spike sorting. Not implemented in the existing files


%% Filter
% If you want to bandpass filter the data
Fs = tp.MetaTags.SamplingFreq;
TiWin = 300 .* Fs;
Hifreq = 6000;

Hd = Butt_bandpass(300,Hifreq,Fs);

rawDatachunk = double(rawData(1:TiWin,:));
rawTimechunk = rawtime(1:TiWin);

rawDatafilt = zeros(size(rawDatachunk));

for ch=1:size(rawDatachunk,2)
    rawDatafilt(:,ch) = filter(Hd,rawDatachunk(:,ch));
end

%% Reduce sampling frequency
% After bandpass filtering you can drop the sampling frequency. Useful for
% further processing
fac = floor(Fs / (Hifreq*2));
rawDatafilt = downsample(rawDatafilt,fac);
rawTimefilt = downsample(rawTimechunk,fac);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER STUFF

% %% FFT
% y = zeros(size(rawDatafilt));
% for ch=1:size(rawDatafilt,2)
%     y(:,ch) = fft(rawDatafilt(:,ch));
% end
% 
% n = length(rawTimefilt);
% f = (0:n-1)*(Fs/n);     % frequency range
% Px = abs(y).^2/n;       % power of the DFT
% 
% 
% %% Plots
% % Time Plot
% figure(1)
% chToPlot = 1;
% plot(rawTimefilt,rawDatafilt(:,chToPlot))
% 
% % Frequency Plot
% figure(2)
% % plot(f(1:floor(n/2)),db(Px(1:floor(n/2),chToPlot)))
% plot(f(1:floor(n/2)),Px(1:floor(n/2),chToPlot))
% xlabel('Frequency')
% ylabel('Power')
% xlim([0 Hifreq])
% 
