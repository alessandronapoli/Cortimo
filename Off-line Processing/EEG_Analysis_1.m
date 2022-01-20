close all; clearvars; clc;

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
%% Quick filter
Fc1 = 0.5;
Fc2 = 150;
Fs = tp.MetaTags.SamplingFreq;
Hd1 = Butt_bandpass_8(Fc1,Fc2,Fs);
% Filter signals
sfi = filter(Hd1,double(rawData));
% sfi2 = filtfilt(Hd1.sosMatrix,Hd1.scaleValues,double(rawData));


Hd1 = Butt_bandpass_8(8,15,Fs);
alpha = filter(Hd1,sfi);

%% Window the data
events = commentsTimeOK;
comments = cellstr(comments);
% labels = [0,1,0,1];

% Make sure to plot the right channels here!!!
% For the preliminary occiptal rhythm acquisition, only Bank B was used and
% only the first 3 channels from Bank B, which are 33, 34, 35

% Pick the channel to plot
chs = {'C3','L-EOG','O1'};
% c = find(strcmp(chs,'O1'));
c = find(strcmp(chs,'C3'));
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
baseline = sfi(1:find(rawtime>=events(1),1)-1,ch);
closed1 = sfi(find(rawtime>=events(1),1):find(rawtime>=events(2),1)-1,ch);
open1 = sfi(find(rawtime>=events(2),1):find(rawtime>=events(3),1)-1,ch);
closed2 = sfi(find(rawtime>=events(1),1):find(rawtime>=events(4),1)-1,ch);
open2 = sfi(find(rawtime>=events(4),1):length(rawtime),ch);


%% Spectra using fft

Bt = fft(baseline);
B = Bt(1:length(baseline)/2+1);
% If you want to get the power spectral density estimate
%Bt(2:end-1) = 2*Bt(2:end-1);
%psdest = 1/length(baseline*Fs)*abs(Bt).^2;
f1 = 0:Fs/length(baseline):Fs/2;

C1t = fft(closed1);
C1 = C1t(1:length(closed1)/2+1);
f2 = 0:Fs/length(closed1):Fs/2;

O1t = fft(open1);
O1 = O1t(1:length(open1)/2+1);
f3 = 0:Fs/length(open1):Fs/2;

C2t = fft(closed2);
C2 = C2t(1:length(closed2)/2+1);
f4 = 0:Fs/length(closed2):Fs/2;

O2t = fft(open2);
O2 = O2t(1:length(open2)/2+1);
f5 = 0:Fs/length(open2):Fs/2;


uplim = 50;
llabels = {'Baseline','Eyes closed','Eyes open'};
%llabels = {'Baseline','Moving feet','Moving hand'};

figure
subplot(5,1,1)
plot(f1,abs(B))
legend(llabels{1})
xlim([0 uplim])
ylabel('Absolute value of FFT')
title('FFT of signals')
subplot(5,1,2)
plot(f2,abs(C1))
legend(llabels{2})
xlim([0 uplim])
subplot(5,1,3)
plot(f3,abs(O1))
legend(llabels{3})
xlim([0 uplim])
ylabel('Absolute value of FFT')
a = axis;
subplot(5,1,4)
plot(f4,abs(C2))
legend('Closed 2')
xlim([0 uplim])
ylabel('Absolute value of FFT')
subplot(5,1,5)
plot(f5,abs(O2))
legend('Open 2')
xlim([0 uplim])
xlabel('Frequency [Hz]');
ylabel('Absolute value of FFT')

pmax = 0;
for fig=1:5
    subplot(5,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
end

for fig=1:5
    subplot(5,1,fig)
    ylim([0,pmax])
end

%% Try spectral analysis with periodogram


% Spectra
[B,w1] = periodogram(baseline,rectwin(length(baseline)),length(baseline));
%f1 = (0:length(B/2)-1).*(Fs/2)/length(B/2);

[C1,w2] = periodogram(closed1,rectwin(length(closed1)),length(closed1));
%f2 = (0:length(C1/2)-1).*(Fs/2)/length(C1/2);

[O1,w3] = periodogram(open1,rectwin(length(open1)),length(open1));
%f3 = (0:length(O1/2)-1).*(Fs/2)/length(O1/2);

[C2,w4] = periodogram(closed2,rectwin(length(closed2)),length(closed2));
%f4 = (0:length(C2/2)-1).*(Fs/2)/length(C2/2);

[O2,w5] = periodogram(open2,rectwin(length(open2)),length(open2));
%f5 = (0:length(O2/2)-1).*(Fs/2)/length(O2/2);


uplim = 50;

figure
subplot(5,1,1)
%plot(f1(1:length(f1)/2),abs(B(1:length(f1)/2,:)))
plot((w1*Fs/2)/pi,10*log10(B))
legend('Baseline')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
title('Periodogram of signals')
subplot(5,1,2)
%plot(f2(1:length(f2)/2),abs(C1(1:length(f2)/2,:)))
plot((w2*Fs/2)/pi,10*log10(C1))
legend('Closed 1')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,3)
%plot(f3(1:length(f3)/2),abs(O1(1:length(f3)/2,:)))
plot((w3*Fs/2)/pi,10*log10(O1))
legend('Open 1')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,4)
%plot(f4(1:length(f4)/2),abs(C2(1:length(f4)/2,:)))
plot((w4*Fs/2)/pi,10*log10(C2))
legend('Closed 2')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,5)
%plot(f5(1:length(f5)/2),abs(O2(1:length(f5)/2,:)))
plot((w5*Fs/2)/pi,10*log10(O2))
legend('Open 2')
xlim([0 uplim])
xlabel('Frequency [Hz]');
ylabel('Power Density [dB/Hz]')

pmax = 0;
for fig=1:5
    subplot(5,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
end

for fig=1:5
    subplot(5,1,fig)
    ylim([0,pmax])
end


%% Try extracting the alhpa band and then plotting

% Group the data
b = alpha(1:find(rawtime>=events(1),1)-1,ch);
c1 = alpha(find(rawtime>=events(1),1):find(rawtime>=events(2),1)-1,ch);
o1 = alpha(find(rawtime>=events(2),1):find(rawtime>=events(3),1)-1,ch);
c2 = alpha(find(rawtime>=events(1),3):find(rawtime>=events(4),1)-1,ch);
o2 = alpha(find(rawtime>=events(4),1):length(rawtime),ch);


% Spectra
[B,w1] = periodogram(b,rectwin(length(b)),length(b));
%f1 = (0:length(B/2)-1).*(Fs/2)/length(B/2);

[C1,w2] = periodogram(c1,rectwin(length(c1)),length(c1));
%f2 = (0:length(C1/2)-1).*(Fs/2)/length(C1/2);

[O1,w3] = periodogram(o1,rectwin(length(o1)),length(o1));
%f3 = (0:length(O1/2)-1).*(Fs/2)/length(O1/2);

[C2,w4] = periodogram(c2,rectwin(length(c2)),length(c2));
%f4 = (0:length(C2/2)-1).*(Fs/2)/length(C2/2);

[O2,w5] = periodogram(o2,rectwin(length(o2)),length(o2));
%f5 = (0:length(O2/2)-1).*(Fs/2)/length(O2/2);


uplim = 50;

figure
subplot(5,1,1)
%plot(f1(1:length(f1)/2),abs(B(1:length(f1)/2,:)))
plot((w1*Fs/2)/pi,10*log10(B))
legend('Baseline')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
title('Periodogram of alpha filtered signals')
subplot(5,1,2)
%plot(f2(1:length(f2)/2),abs(C1(1:length(f2)/2,:)))
plot((w2*Fs/2)/pi,10*log10(C1))
legend('Closed 1')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,3)
%plot(f3(1:length(f3)/2),abs(O1(1:length(f3)/2,:)))
plot((w3*Fs/2)/pi,10*log10(O1))
legend('Open 1')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,4)
%plot(f4(1:length(f4)/2),abs(C2(1:length(f4)/2,:)))
plot((w4*Fs/2)/pi,10*log10(C2))
legend('Closed 2')
xlim([0 uplim])
ylabel('Power Density [dB/Hz]')
subplot(5,1,5)
%plot(f5(1:length(f5)/2),abs(O2(1:length(f5)/2,:)))
plot((w5*Fs/2)/pi,10*log10(O2))
legend('Open 2')
xlim([0 uplim])
xlabel('Frequency [Hz]');
ylabel('Power Density [dB/Hz]')

pmax = 0;
for fig=1:5
    subplot(5,1,fig)
    a = axis;
    pmax = max(pmax,a(4));
end

for fig=1:5
    subplot(5,1,fig)
    ylim([0,pmax])
end


% Let's see what happens in time
figure
subplot(2,1,1)
plot(b(1:30e4))
subplot(2,1,2)
plot(c1(1:30e4))
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
end

for fig=1:2
    subplot(2,1,fig)
    ylim([pmin,pmax])
    pmin = min(pmin,a(3));
end

figure
time = 5;
ll = Fs*time;
subplot(2,1,1)
plot((1:ll)/Fs,baseline(1:ll))
subplot(2,1,2)
plot((1:ll)/Fs,closed1(1:ll))
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
end

for fig=1:2
    subplot(2,1,fig)
    ylim([pmin,pmax])
    pmin = min(pmin,a(3));
end




