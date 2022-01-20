function [feats,chs] = EpochsToLFPFeatures(file,Fs,timeWin,overlap,chs,freqBin,freqRange,condLab,isCont)

addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('Aux Funcs'));
% [~,nm,~] = fileparts(file);
% nm = extractAfter(nm,'-');
tp = load(file);
tpp = fieldnames(tp);
Epochs = tp.(tpp{1});

% condLab = -1:size(Epochs,1)-2;
freq = freqRange(1):freqBin:freqRange(2)-1;
Fsdown = 2000;
factor = Fs/Fsdown;

HdMain = Butt_bandpass(0.5,1000,Fs);
HdMain.PersistentMemory = true;

if(isempty(chs))
    chs = 1:size(Epochs{1,1},1);
end

% Windows in number of samples
w = timeWin * Fsdown;
% Calculate the step size
stepby = w*(1-overlap);

% Set a transitory time window in which there might be a transition
% between intent or trials or to account for subject's delays in responding
% to visual/auditory instructions or delays in motion execution.
% transTime = 0.5;

%% FEATURE EXTRACTION
feats = [];
rshfeats = [];

if isCont == 0
    
    for i=1:length(condLab)
        
        if condLab(i) == -1
            cond = 1;
        elseif condLab(i) == 0
            cond = 2;
        elseif condLab(i) == 1
            cond = 3;
        elseif condLab(i) == 2
            cond = 4;
        end
        
        
        TPfeats = zeros(length(freq),length(chs),1);
        
        if(~isempty(Epochs{cond,1}))
            for trial=1:size(Epochs,2)
                reset(HdMain);
                if(~isempty(Epochs{cond,trial}))
                    data = Epochs{cond,trial}';
                    data = filter(HdMain,data(:,chs));
                    
                    data = downsample(data,factor);
                    here = w;
                    stepin = w*overlap;
                    % for window=1:round((size(data,1)/stepby))-1
                    if(round((size(data,1)-w-1)/stepby)>1)
                        for window=1:round((size(data,1)-w-1)/stepby)
                            if(window == 1)
                                TPfeats(:,:,window) = periodogram(data(1:w,:),[],freq,Fsdown);
                            else
                                ti = here-stepin+1;
                                TPfeats(:,:,window) = periodogram(data(ti:ti+w-1,:),[],freq,Fsdown);
                                here = ti+w-1;
                                
                            end
                        end
                        
                        
                        
                        % Features have been extracted and now can be concatenated to compute a
                        % single feature vector for each time point
                        rshfeats = reshape(TPfeats,size(TPfeats,3),size(TPfeats,1)*size(TPfeats,2));
                        % Finally add the labels to the trials
                        rshfeats = [rshfeats, ones(size(rshfeats,1),1).*condLab(i)];
                        
                        feats = [feats; rshfeats];
                        
                        TPfeats = zeros(length(freq),length(chs),1);
                    end
                end
                
            end
            
        end
        
        
    end
    
else
    for trial=1:size(Epochs,1)
        reset(HdMain);
        if(~isempty(Epochs{trial,1}))
            
            
            
            data = Epochs{trial,1}'; 
            data = filter(HdMain,data(:,chs));
            
            data = downsample(data,factor); 
            
            here = w; 
            stepin = w*overlap;
            % for window=1:round((size(data,1)/stepby))-1
            if(round((size(data,1)-w-1)/stepby)>1)
                for window=1:round((size(data,1)-w-1)/stepby)
                    if(window == 1)
                        TPfeats(:,:,window) = periodogram(data(1:w,:),[],freq,Fsdown);
                    else
                        ti = here-stepin+1; 
                        TPfeats(:,:,window) = periodogram(data(ti:ti+w-1,:),[],freq,Fsdown);
                        here = ti+w-1;
                        
                    end
                end
                
                
                % Features have been extracted and now can be concatenated to compute a 
                % single feature vector for each time point 
                rshfeats = reshape(TPfeats,size(TPfeats,3),size(TPfeats,1)*size(TPfeats,2));
                
                
                % Resample continuous labels to match length of LFP feats
                labels = zeros(window,2);
                
                
                reSampFactor = length(condLab.ContLabels{trial,1})/window;
                xq = 1:reSampFactor:length(condLab.ContLabels{trial,1});
                
                if length(xq) < window
                    xq(window) = length(condLab.ContLabels{trial,1});
                end
                
                % resampling using 1-D interpolation
                foo = interp1(condLab.ContLabels{trial,1}(:,1), xq, 'spline');
                labels(:,1) = foo(1:window);
                foo = interp1(condLab.ContLabels{trial,1}(:,2), xq, 'spline');
                labels(:,2) = foo(1:window);
                
                
                
                % Finally add the labels to the trials
                rshfeats = [rshfeats, labels];
                
                feats = [feats; rshfeats];
                
                TPfeats = zeros(length(freq),length(chs),1);
            end
        end
        
    end
    
end

end

