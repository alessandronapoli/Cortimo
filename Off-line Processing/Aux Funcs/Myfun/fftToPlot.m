
%% Frequency analysis

function [Freq, absSig] = fftToPlot(sig,Fs)

t = fft(sig);
absSig = abs(t(1:length(sig)/2+1));
Freq = 0:Fs/length(sig):Fs/2;

end
    