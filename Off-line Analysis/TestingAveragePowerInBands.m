
Fs = 2000;
t = (0:1/Fs:900)';


fff = [21,45,340,600];

x1 = 5.*randn(size(t,1),length(fff)) + 3.* cos(2*pi.*fff.*t);

xs = sum(x1,2);

FS = 10;
fr = 0:FS:1000;

wind = 0.5 * Fs;
PD = zeros(length(fr),size(xs,2),round(length(xs)/wind));
PW = zeros(length(fr),size(xs,2),round(length(xs)/wind));

trial = 1;
for w=1:wind:length(xs)-wind
    
[PD(:,:,trial),f] = periodogram(xs(w:w+wind-1,:),[],fr,Fs);
[PW(:,:,trial),ff] = pwelch(xs(w:w+wind-1,:),[],[],fr,Fs);
PBB(:,trial) = bandpower(PD(:,:,trial),fr,'psd');

trial = trial + 1;


end

figure
plot(f,db(PD(:,:,2)))
hold on
plot(ff,db(PW(:,:,2)))
plot(ff,db(PBB(:,2)))



