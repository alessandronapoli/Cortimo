%% Load the EEG acquisition file

clearvars; clc;
[file,path] = uigetfile('*.edf','Select edf file to open');
fullname = [path,filesep,file];
addpath('AuxFunc');
addpath('NPMK');

fd = read_edf_file(fullname);

signal_to_plot = 9;
nspr = fd.nSamplesPerRec(signal_to_plot);
total_samples = sum(fd.nSamplesPerRec);
fs = nspr/fd.duration;
t_err_correction = 535.4 - 516.2;
ns = fd.ns;
% t = (1:length(s3))/fs + t_err_correction; % not sure what this error
% correction is for
% t = t / 60; % If you want time in minutes


%% 
EEGMatrix = zeros(fd.nDataRec*fs,ns);

for ch=1:ns
    signal_to_plot = ch;
    total_samples = sum(fd.nSamplesPerRec);
    nspr = fd.nSamplesPerRec(signal_to_plot);
    s1 = reshape(fd.data,total_samples,[]);
    iStart = sum(fd.nSamplesPerRec(1:signal_to_plot-1)) + 1;
    iEnd = iStart + nspr - 1;
    s2 = s1(iStart:iEnd,:);
    s3 = reshape(s2,[],1)';
   

    EEGMatrix(:,ch) = s3;
    
end

 % t = (0:length(s3)-1)/fs; 
 chLabels = cellstr(fd.labels);
 
 
 %% Preprocess the data
 
 % Notch filter for 60Hz
 d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
 EEGMatrix = filtfilt(d,EEGMatrix); 
 
 % detrend the data
 EEGMatrix = detrend(EEGMatrix);
 
 % Select channels
 EEGMatrix = EEGMatrix(:,[1:19,23:32]);
 channels = chLabels([1:19,23:32]);
%% Try comparing these signals with Blackorock acquisitions

% Disregard the first t minutes
t = 30;
samps = t*60*fs;
EEGMatrix = EEGMatrix(samps:end,:);

% Re-reference the channels to A1. This is the initial setting for the
% Blackrock system
nRef1 = find(contains(channels,'A1'));
nRef2 = find(contains(channels,'A2')); % A2 is really nosiy for ID 0001

% EEGMatrixRR = EEGMatrix - mean(EEGMatrix(:,[nRef1,nRef2]));
EEGMatrixRR = EEGMatrix - EEGMatrix(:,nRef1);


% [MeanSpectra,H] = mySpecAnalysis_1(EEGMatrix,0.5,0.5,fs,2,1,channels);
[MeanSpectra,H,Fr] = mySpecAnalysis_1(EEGMatrixRR,0.5,0.5,fs,2,1,channels);
 

% Blackrock data

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

sig = rawData';
% Analyze data time windows from homologous electrodes
BRCLabels = {'FP1','FP2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T3','T4','T5','T6','Z-Gr','Fz','Cz','Pz','A1','A2','C2',...
    'C6','FC4','CP4','CP2','CP6','FC2','FC6'};
Ch = 1:30;

% Filter
Hd1 = Butt_bandpass(0.5,250,Fs);
sig = filter(Hd1,sig);
% Downsample to 500 sps
Fsdown = 500;
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

 %%%%%%%%%%%%%%% RE-REFERENCE %%%%%%%%%%%%%%%%%%
% Signals were recorded with A1 as only reference
% Let's use the average between A1 and A2
% [~,~,a2] = intersect({'A1','A2'},ChLabels);
thisID = find(contains(BRCLabels,'A1'));
sig = sig - sig(:,thisID);

% 
[MeanSpectraBR,H,Fr] = mySpecAnalysis_1(sig,0.5,0.5,Fsdown,2,1,BRCLabels);



figure
plot(sig(:,1))
hold on
plot(EEGMatrixRR(1:fs*2000,1))

% Get spectrum of two equivalent channels
ch = 16;
SPNK = periodogram(EEGMatrix(:,ch),[],Fr,fs);
SPBR = periodogram(sig(:,ch),[],Fr,Fsdown);
figure
plot(Fr,10*log10(SPNK))
hold on
plot(Fr,10*log10(SPBR))

%% Look for homologous electrode differences
% Group the elctrodes in hemicraniectomy and skull
hemi = {'F4','C4','P4','T4'};
skull = {'F3','C3','P3','T3'};

additionalHemi = {'C2','C6','FC4','CP4','CP2','CP6','FC2','FC6'};



idHemi = find(contains(channels,hemi));
idSkull = find(contains(channels,skull));
addHemi = find(contains(channels,additionalHemi));

% [~,~,idHemi] = intersect(hemi,channels);
% [~,~,idSkull] = intersect(skull,channels);
idHemi = sort(idHemi);
idSkull = sort(idSkull);

% idSkull = find(strcmp(ChLabels,skull)); 
% [~,~,addHemi] = intersect(additionalHemi,channels);

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


%%
% Plot channel Cp4
chtoPlot = find(contains(channels,'CP4'));
figure
plot((0:length(EEGMatrixRR)-1)./fs,EEGMatrixRR(:,chtoPlot));



 %% Apply filter bank

% From "Hemicraniectomy: A New Model for Human Electrophysiology with High 
% Spatio-Temporal Resolution"
% 75 pass-bands in 2Hz increments from 0-150 Hz using Gaussian shaped
% filterring windows. Gaussian window standard deviation was chosen as 10%
% of its center frequency. Then the analytic amplitude, which is the
% absolute value of the Hilbert transform of the signal was used to create
% a grand average time-frequency feature (ERPs in the paper).

% Approach 2 with just a regular Butterworth filter bank with 2Hz frequency
% bins.

% Build Gaussian windows and Butterworth filters

wid = 2;
highFreqlim = 100;
len = length(EEGMatrix);
gsWin = zeros(len,highFreqlim/wid);
ButCoef = cell(2,highFreqlim/wid);
cc = 1;
Fbands = 2:wid:highFreqlim;
for frr=0:wid:highFreqlim-wid
    cFreq = frr+(wid/2);
   gsWin(:,cc) = myGaussian(1:len,cFreq,0.1*cFreq);
   if(frr==0)
       fnorm = frr+wid / (fs/2); % Low-pass at frr+wid
   else
       fnorm = [frr,frr+wid] / (fs/2);
   end
   [ButCoef{1,cc},ButCoef{2,cc}] = butter(4,fnorm);
    cc = cc+1;
end


%% Perform the filtering

% FFT of the data
FFTMatrix = fft(EEGMatrix);
%GaussFiltered = zeros(size(FFTMatrix));
%ButterFiltered = zeros(size(EEGMatrix));
GaussFiltered = cell(1,cc-1);
ButterFiltered = cell(1,cc-1);

cc = 1;
for frr=0:wid:highFreqlim-wid
    
    % Gaussian Window Filter
    % In the frequency domain
    for column=1:size(FFTMatrix,2)
    GaussFiltered{cc}(:,column) = gsWin(:,cc) .* FFTMatrix(:,column);
     % Butterworth Filter
    ButterFiltered{cc}(:,column) = filter(ButCoef{1,cc},ButCoef{2,cc},...
        EEGMatrix(:,column));
    end
    
   
    cc = cc+1;
end

channelToPlot = 1;
FF = (0:length(t)-1)*(fs/length(t));
plot(FF(1:length(FF)/2),abs(FFTMatrix(1:length(FFTMatrix)/2,channelToPlot)))
hold on
thisfft = fft(ButterFiltered{5});
plot(FF(1:length(FF)/2),abs(thisfft(1:length(thisfft)/2,channelToPlot)))
plot(FF(1:length(FF)/2),abs(GaussFiltered{5}(1:length(GaussFiltered{5})/2,channelToPlot)))


%% Get the average power in frequency bands and display
PP = zeros(length(Fbands),ch);
for band=1:length(Fbands)

    % Dividing by wid we convert this into the power spectral density
    PP(band,:) = mean(ButterFiltered{band}.^2) / wid;

end

figure
surf(db(PP(1:end,1:23))')
xlabel('Frequency Hz')
ylabel('Channel')
view([0,90])
% caxis([-10,10]);
set(gca,'YTick',(1:23));
% yticklabels(chLabels);
set(gca,'XTick',1:50)
set(gca,'XTicklabel',Fbands)
title(['Subject ',fd.loc_pt_info(1:8)]);
% axis tight;
%set(gca,'ydir','reverse')

figure
stem(Fbands(3:end),PP(3:end,1:23),'*')

% Split the channels in left vs. right

figure
subplot(1,2,1)
stem(Fbands(3:end),PP(3:end,1:2:16),'*')
legend(chLabels(1:2:16))
subplot(1,2,2)
stem(Fbands(3:end),PP(3:end,2:2:16),'*')
legend(chLabels(2:2:16))