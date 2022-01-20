function [thresholds] = calcthresholds(traindata,labels,channel,freqRange,freqbwidth,Fs)
% [thresholds] = calcthresholds(channel, frequency, bandwidth, samplingrate)
%
% calculate decision thresholds for a given subject, session, channel, 
% center frequency, bandwidth, and samplingrate

% traindata is a 3D matrix whose dimensions are ch x freq x trial
% labels is a 1D matrix whose dimension matches the length(trial). In other
% words labels contains the class assigned to a specific trial and it is
% used to train the classification thresholds

% Identify the number of classes in the training data set
classes = unique(labels);

avgspectrum = cell(1,length(classes));
ffvec = 0:freqbwidth:Fs/2-1/Fs/2;
% avgspectrum1=[];
% avgspectrum2=[];
% avgspectrum3=[];
% avgspectrum4=[];
for cur_trial=1:size(traindata,3)
   
    trialdata = traindata(:,:,cur_trial);
   % trialdata=signal(trialidx, channel);
   % truncate the trial (just a precaution ... should always be of equal length)
   % trialdata=trialdata(1:triallength);
   % calculate the power spectral density
   % [Pxx,F,Pxxc,F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P,DFLAG)
   % [Pxx, Pxxc, F] = PSD(trialdata, triallength, samplingrate, triallength, 0, 0.95, 'linear');
   [Pxx, F] = periodogram(trialdata,[],ffvec,Fs);
   cur_targetcode = labels(cur_trial);
   
   avgspectrum{cur_targetcode} = cat(2,avgspectrum{cur_targetcode},Pxx);
   
%    if (cur_targetcode == 1)
%       avgspectrum1=cat(2, avgspectrum1, Pxx);
%    end
%    if (cur_targetcode == 2)
%       avgspectrum2=cat(2, avgspectrum2, Pxx);
%    end
%    if (cur_targetcode == 3)
%       avgspectrum3=cat(2, avgspectrum3, Pxx);
%    end
%    if (cur_targetcode == 4)
%       avgspectrum4=cat(2, avgspectrum4, Pxx);
%    end
end

res = [];
meanFeat = [];
for cond=1:length(classes)
    
% calculate average spectra for each condition
res(cond)=mean(avgspectrum{cond}, 2);
res2=mean(avgspectrum2, 2);
res3=mean(avgspectrum3, 2);
res4=mean(avgspectrum4, 2);

% calculate feature means
idxs=find((F >= TopFreq) & (F < frequency+bandwidth/2));
meanFeat(cond)=mean(res1(idxs));
mean2=mean(res2(idxs));
mean3=mean(res3(idxs));
mean4=mean(res4(idxs));

% set decision thresholds between feature means
% of course this could be done much fancier, but this is just a demo
thresholds(1)=mean([mean1, mean2]);
thresholds(2)=mean([mean2, mean3]);
thresholds(3)=mean([mean3, mean4]);

end