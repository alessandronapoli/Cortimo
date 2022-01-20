close all; clearvars; clc;

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
rawData = tp.Data';
rawtime = (0:length(rawData)-1)*(1/tp.MetaTags.SamplingFreq); % in seconds


%% Open NEV file
% Stores spike count
nevf = [filename(1:end-3),'nev'];
% If you want to get spikes
% openNEV('report','read',nevf,'uV','nomat','nosave');

% If you want to get the comments
openNEV('report','read',nevf,'nomat','nosave');

%events = NEV.Data.Spikes.TimeStamp;
%electrode = NEV.Data.Spikes.Electrode;
%Unit = NEV.Data.Spikes.Unit;
comments = NEV.Data.Comments.Text;
commentsTimeStart = NEV.Data.Comments.TimeStampStartedSec;
commentsTimeOK = NEV.Data.Comments.TimeStampSec;
%% Filters
Fc1 = 1;
Fc2 = 100;
Fs = tp.MetaTags.SamplingFreq;
Hd1 = Butt_bandpass_8(Fc1,Fc2,Fs);
% Filter signals
sfi = filter(Hd1,double(rawData));
% sfi2 = filtfilt(Hd1.sosMatrix,Hd1.scaleValues,double(rawData));\

Hd1 = Butt_bandpass_8(8,15,Fs);
alpha = filter(Hd1,sfi);


%% Downsample the signals
% Set new Fs at 500Hz




%% Window the data
events = commentsTimeOK;
comments = cellstr(comments);
% labels = [0,1,0,1];

% Make sure to plot the right channels here!!!
% For the preliminary occiptal rhythm acquisition, only Bank B was used and
% only the first 3 channels from Bank B, which are 33, 34, 35

% Pick the channel to plot
chs = {'C3','L-EOG','O1'};
c = find(strcmp(chs,'O1'));
% c = find(strcmp(chs,'C3'));
ch = c + 32;

figure
plot(rawtime,rawData(:,ch),'black')
hold on
lims = [min(rawData(:,ch)),max(rawData(:,ch))];
X = [events;events];
Y = ones(size(X)).*double(lims');
line(X,Y)
xlabel('Time [s]')
ylabel('Amplitude [uV]')


% Group the data
sig = cell(length(comments),1);
for gg=1:length(comments)
    if gg<length(comments)
sig{gg} = sfi(find(rawtime>=events(gg),1):find(rawtime>=events(gg+1),1)-1,ch);
    else
sig{gg} = sfi(find(rawtime>=events(gg),1):length(rawtime),ch);
    end
end

%% Spectra using fft
SIG = cell(length(comments),2);
for gg=1:length(comments)
    t = fft(sig{gg});
    SIG{gg,2} = t(1:length(sig{gg})/2+1);
    SIG{gg,1} = 0:Fs/length(sig{gg}):Fs/2;
end

uplim = 50;

figure
len = length(events);
for fig=1:len
    subplot(len,1,fig)
    plot(SIG{fig,1},abs(SIG{fig,2}))
    xlim([0 uplim])
    ylabel('Absolute value of FFT')
    title('FFT of signals')
    xlabel('Frequency [Hz]');
    legend(comments{fig})

end


pmax = 0;
for fig=1:len
    subplot(len,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
end

for fig=1:len
    subplot(len,1,fig)
    ylim([0,pmax])
end



% %% Try spectral analysis with periodogram
% 
% % Spectra
% [B,w1] = periodogram(baseline,rectwin(length(baseline)),length(baseline));
% %f1 = (0:length(B/2)-1).*(Fs/2)/length(B/2);
% 
% [C1,w2] = periodogram(closed1,rectwin(length(closed1)),length(closed1));
% %f2 = (0:length(C1/2)-1).*(Fs/2)/length(C1/2);
% 
% [O1,w3] = periodogram(open1,rectwin(length(open1)),length(open1));
% %f3 = (0:length(O1/2)-1).*(Fs/2)/length(O1/2);
% 
% [C2,w4] = periodogram(closed2,rectwin(length(closed2)),length(closed2));
% %f4 = (0:length(C2/2)-1).*(Fs/2)/length(C2/2);
% 
% [O2,w5] = periodogram(open2,rectwin(length(open2)),length(open2));
% %f5 = (0:length(O2/2)-1).*(Fs/2)/length(O2/2);
% 
% 
% uplim = 50;
% 
% figure
% subplot(5,1,1)
% %plot(f1(1:length(f1)/2),abs(B(1:length(f1)/2,:)))
% plot((w1*Fs/2)/pi,10*log10(B))
% legend('Baseline')
% xlim([0 uplim])
% ylabel('Power Density [dB/Hz]')
% title('Periodogram of signals')
% subplot(5,1,2)
% %plot(f2(1:length(f2)/2),abs(C1(1:length(f2)/2,:)))
% plot((w2*Fs/2)/pi,10*log10(C1))
% legend('Closed 1')
% xlim([0 uplim])
% ylabel('Power Density [dB/Hz]')
% subplot(5,1,3)
% %plot(f3(1:length(f3)/2),abs(O1(1:length(f3)/2,:)))
% plot((w3*Fs/2)/pi,10*log10(O1))
% legend('Open 1')
% xlim([0 uplim])
% ylabel('Power Density [dB/Hz]')
% subplot(5,1,4)
% %plot(f4(1:length(f4)/2),abs(C2(1:length(f4)/2,:)))
% plot((w4*Fs/2)/pi,10*log10(C2))
% legend('Closed 2')
% xlim([0 uplim])
% ylabel('Power Density [dB/Hz]')
% subplot(5,1,5)
% %plot(f5(1:length(f5)/2),abs(O2(1:length(f5)/2,:)))
% plot((w5*Fs/2)/pi,10*log10(O2))
% legend('Open 2')
% xlim([0 uplim])
% xlabel('Frequency [Hz]');
% ylabel('Power Density [dB/Hz]')
% 
% pmax = 0;
% for fig=1:5
%     subplot(5,1,fig)
%     a = axis;
%     pmax = max(pmax,a(4));
% end
% 
% for fig=1:5
%     subplot(5,1,fig)
%     ylim([0,pmax])
% end


%% Try extracting the alhpa band and then plotting

% Group the data
Bsig = cell(length(comments),1);
for gg=1:length(comments)
    if gg<length(comments)
Bsig{gg} = alpha(find(rawtime>=events(gg),1):find(rawtime>=events(gg+1),1)-1,ch);
    else
Bsig{gg} = alpha(find(rawtime>=events(gg),1):length(rawtime),ch);
    end
end

%% Spectra using fft
BSIG = cell(length(comments),2);
for gg=1:length(comments)
    t = fft(Bsig{gg});
    BSIG{gg,2} = t(1:length(Bsig{gg})/2+1);
    BSIG{gg,1} = 0:Fs/length(Bsig{gg}):Fs/2;
end


uplim = 50;

figure
len = length(events);
for fig=1:len
    subplot(len,1,fig)
    plot(BSIG{fig,1},abs(BSIG{fig,2}))
    xlim([0 uplim])
    ylabel('Absolute value of FFT')
    title('FFT of signals')
    xlabel('Frequency [Hz]');
    legend(comments{fig})

end


pmax = 0;
for fig=1:len
    subplot(len,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
end

for fig=1:len
    subplot(len,1,fig)
    ylim([0,pmax])
end


% Let's see what happens in time
figure
sig1 = Bsig{1};
sig2 = Bsig{2};
subplot(2,1,1)
plot(sig1(1:30e4))
subplot(2,1,2)
plot(sig2(1:30e4))
xlabel('# of samples');
ylabel('Amplitude [uV]')
title('Amplitudes in the alpha band')

pmax = 0;
pmin = 0;
for fig=1:2
    subplot(2,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
    pmin = min(pmin,a(3));
    legend(comments{fig})
end

for fig=1:2
    subplot(2,1,fig)
    ylim([pmin,pmax])
    pmin = min(pmin,a(3));
end

figure
sig1 = sig{1};
sig2 = sig{2};
time = 5;
ll = Fs*time;
subplot(2,1,1)
plot((1:ll)/Fs,sig1(1:ll))
subplot(2,1,2)
plot((1:ll)/Fs,sig2(1:ll))
xlabel('Time [s]');
ylabel('Amplitude [uV]')
title('Amplitudes with all frequencies')
pmax = 0;
pmin = 0;
for fig=1:2
    subplot(2,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
    pmin = min(pmin,a(3));
    legend(comments{fig})
end

for fig=1:2
    subplot(2,1,fig)
    ylim([pmin,pmax])
    pmin = min(pmin,a(3));
end



%% Run the analysis to simulate real-time acquisition
newFs = 1000;
n = Fs/newFs;
Hd1 = Butt_bandpass(0.5,400,Fs);
tm = filter(Hd1,rawData);
newData = downsample(tm,n);
newtime = downsample(rawtime,n);

% Desing a filter to use in RT
HdRT = Butt_bandpass(1,100,newFs);
HdRT2 = Butt_bandpass(8,15,newFs);

% Simulate now the RT execution
step = 0.2 * newFs; % Grab data every step size

% Filter the entire sequence
% filtered = filter(HdRT2,newData);
filtered = filter(HdRT,newData);
alpha = filter(HdRT2,filtered);

%% This for loop simulates RT data polling
Buffer = zeros(size(newData));
win_time = zeros(round(length(newData)/step),1);


for cycle=1:step:length(newData)-step
    
    %if(cycle+step-1 <= length(newData))
    % seldata = newData(cycle:(cycle+step-1),:);
     seldata = filter(HdRT,newData(cycle:(cycle+step-1),:));
    % Buffer = cycleBuffermulti(Buffer, seldata);
     Buffer = cycleBuffermulti(Buffer, filter(HdRT2,seldata));
    % Buffer = cycleBuffermulti(Buffer, filter(HdRT,seldata));

     win_time(cycle) = (cycle+step-1)*newFs;
   % else
   % Buffer = cycleBuffermulti(Buffer, filter(HdRT,newData(cycle:length(newData),:)));
   % en  
   
   if(cycle+2*step >=length(newData))
      % Buffer =  cycleBuffermulti(Buffer, newData(cycle+step:length(newData),:)); 
       % seldata = newData(cycle+step:length(newData),:);
      seldata = filter(HdRT,newData(cycle+step:length(newData),:));
      Buffer =  cycleBuffermulti(Buffer, filter(HdRT2,seldata));
     % Buffer =  cycleBuffermulti(Buffer, filter(HdRT,seldata));
     win_time(cycle+1) = length(newData)*newFs;
   end
    
   % Extract the alpha power content in each time window with time
   % information
  
   
   
   
end

figure
% plot(newtime,filtered(:,ch),'--*')
plot(newtime,alpha(:,ch),'--*')
hold on
plot(newtime,Buffer(:,ch),'--o')


figure
% [f,ft] = fftToPlot(filtered(:,ch),newFs); 
[f,ft] = fftToPlot(alpha(:,ch),newFs); 
% plot(f,20*log10(ft),'--*')
plot(f,ft,'--*')
hold on
[f,ft] = fftToPlot(Buffer(:,ch),newFs); 
% plot(f,20*log10(ft),'--o')
plot(f,ft,'--o')
xlim([0,100])
