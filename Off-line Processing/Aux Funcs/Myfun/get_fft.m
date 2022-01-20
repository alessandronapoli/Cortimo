% This function will calculate a real-time FFT
function [F,P] = get_fft(signal,method,Fs,HighFreq,varargin)

if(~isempty(HighFreq))
n = length(signal) / HighFreq;
else
    n = 1;
end

switch method
    
    case 1
        
        X = fft(signal);
        F = (0:round(length(X)/2)-1) * Fs/(2*length(X));
        P = X(1:round(length(X)/2));
    case 2
        if(isempty(varargin))
        freq = 0:((length(signal)/n)/2-1);
        else
            freq = (0:cell2mat(varargin):(length(signal)/n)/2-1);
        end
        
        signal = double(signal);
        [P,F] = periodogram(signal,[],freq,Fs);
        
end
        
        

end

