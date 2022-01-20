function [MeanSpectra,H,Fr] = mySpecAnalysis_1(sig,TimeWin,overlap,Fs,...
    FreqBinWidth,isPlot,ChLabels)

% Get spectra of time windows
w = TimeWin * Fs;
% Select overlap between windows
% overlap = 0.5;

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


Fstep = FreqBinWidth;
Fr = 0:Fstep:Fs/2-1/Fs/2;
MeanSpectra = zeros(length(Fr),size(sig,2),size(epoched,3));

for trial=1:size(epoched,3)
    
    MeanSpectra(:,:,trial) = periodogram(epoched(:,:,trial),[],Fr,Fs);

end

% Average the Power Spectral Density
MeanSpectra = mean(MeanSpectra,3);

if (isPlot == 1)
    figure
H = surf(1:size(MeanSpectra,2),Fr,10*log10(MeanSpectra));
view(90,-90)
ylabel('Frequency Hz')
xlabel('Channels')
xticks(1:size(sig,2));
xticklabels(ChLabels);
colormap(jet);
% caxis([-20 30]); % Change limits
axis('tight')
end
    
end