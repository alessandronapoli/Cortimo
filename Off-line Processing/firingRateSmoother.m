
function [smoothed] = firingRateSmoother(vector,binwidth,tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [smoothed] = firingRateSmoother(spikeVector,binwidth,timeConstant)
%   Written by Mijail Serruya, 2/21/2005
% tau is time constant
% vector is spike vector, vector of spike counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tau = 0.02; %0.01 % 0.05/0.02! good at showing; if making discrete thresholder, exceed running average?
up = tau*(-1);
a =  0; b = 4;
sprintf('Bindwidth is %d msec, time constant is %d msec .',round(binwidth*1000),round(tau*1000));
n = 0;
nsk=1;
k = 1; m = 1;
while k == 1
    now = vector(m);
    if now > 0
        k = 2;
    end
    m = m+1;
end
first = vector(m); last = vector(size(vector,1));
first = (round(first*100))/100; last = (round(last*100))/100;
elapsed = last - first;
numspikes = length(vector);
meanrate = round(numspikes/elapsed);
meanisi = round(mean(diff(vector))*100);
bins = first:binwidth:last;
binned = hist(vector,bins);
binsizebinnedsize = [size(bins,2) size(binned,2)];
nxbin=  [min(binned) max(binned)];
k = 0;
allk = [];
for ii=1:size(binned,2)
    k = (k+binned(ii))*(1 - exp(binwidth/up));%/tau
    allk = [allk k];
end
leaky = allk;
normalize = leaky./(max(leaky));
smoothed = normalize;
nsk=nsk+1;
sizenormalize = length(normalize);
nindex = bins;
nindex = nindex(:);
vector = [ones(size(vector)).*3 vector(:)];
subplot(2,1,1), plot(nindex,normalize,'r-',vector,ones(size(vector)).*3,'.k')
%axis([0 last a b])
fs = 1;
filter = ones(fs,1)./fs;
subplot(2,1,2), bar(conv(filter,normalize))
title(sprintf('Binwidth %d msec, time constant %d msec',binwidth*1000,tau*1000))
