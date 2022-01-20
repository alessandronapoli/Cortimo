function Combfeats = EpochsToCMBFeatures(file,Fs,timeWin,overlap,chs,freqBin,freqRange,condLab,Spkchs,SpksTW,isCont)

addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('Aux Funcs'));
% [~,nm,~] = fileparts(file);
% nm = extractAfter(nm,'-');
tp = load(file);
tpp = fieldnames(tp);
EpochsLFP = tp.(tpp{1});

nfile = strrep(file,'Volts','TS');
tppp = load(nfile);
ttpp = fieldnames(tppp);
Nsamples = tppp.(ttpp{1});

sfile = strrep(nfile,'TS','Spks');
tppp = load(sfile);
ttpp = fieldnames(tppp);
EpochsSPK = tppp.(ttpp{1});


% condLab = -1:size(Epochs,1)-2;
freq = freqRange(1):freqBin:freqRange(2)-1;
Fsdown = 2000;
factor = Fs/Fsdown;

HdMain = Butt_bandpass(100,1000,Fs);
HdMain.PersistentMemory = true;

if(isempty(chs))
    chs = 1:size(EpochsLFP{1,1},1);
end


if isempty(Spkchs)
    if(~isempty(EpochsSPK{1,1}))
        Spkchs = unique(EpochsSPK{1,1}(1,:));
    else
        if(~isempty(EpochsSPK{2,1}))
            Spkchs = unique(EpochsSPK{2,1}(1,:));
        else
            if(~isempty(EpochsSPK{3,1}))
                Spkchs = unique(EpochsSPK{2,1}(1,:));
            end
        end
    end
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
Combfeats = [];
Fss = 30e3;
Spkfeats = [];
binsinthepast = round((stepby/Fsdown) / SpksTW);
tfe = [];
avtfe = [];

TPfeats = zeros(length(freq),length(chs),1);

stepin = w*overlap;

% Express bins in original sampling frequency at 30KHz
sampStep = SpksTW * Fss;
ffact = Fss / Fsdown;


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
        
        rshfeats = [];
        
        if(~isempty(EpochsLFP{cond,1}))
            
            
            for trial=1:size(EpochsLFP,2)
                reset(HdMain);
                
                if(~isempty(EpochsLFP{cond,trial}) && ~isempty(EpochsSPK{cond,trial}))
                    data = EpochsLFP{cond,trial}';
                    data = filter(HdMain,data(:,chs));
                    
                    spdata = EpochsSPK{cond,trial};
                    data = downsample(data,factor);
                    here = w;
                    
                    if(round((size(data,1)-w-1)/stepby)>1)
                        
                        TPfeats = zeros(length(freq),length(chs),1);
                        
                        % for window=1:round((size(data,1)/stepby))-1
                        for window=1:round((size(data,1)-w-1)/stepby)
                            
                            if(window == 1)
                                TPfeats(:,:,window) = periodogram(data(1:w,:),[],freq,Fsdown);
                                
                                timeBins = Nsamples{cond,trial}(1):sampStep:Nsamples{cond,trial}(1)+w*ffact;
                                Spkf = zeros(size(timeBins,2)-1,size(Spkchs,2));
                                
                                for step =1:length(timeBins)-1
                                    idxs = find((spdata(2,:) >= timeBins(step)) & ...
                                        (spdata(2,:) < timeBins(step+1)-1));
                                    for c=1:size(Spkchs,2)
                                        ch = Spkchs(c);
                                        Spkf(step,c) = length(find(spdata(1,idxs) == ch));
                                    end
                                end
                                % tfe = [totalfeats;[tpFeats,ones(size(tpFeats,1),1).*condLab(cond)]];
                                tfe = [tfe;Spkf(size(Spkf,1),:)];
                                if(size(Spkf,1)-(binsinthepast-1)>0)
                                    avtfe = [avtfe;mean(Spkf(size(Spkf,1)-(binsinthepast-1):size(Spkf,1),:))];
                                else
                                    avtfe = [avtfe;mean(Spkf(1:size(Spkf,1),:))];
                                end
                            else
                                ti = here-stepin+1;
                                TPfeats(:,:,window) = periodogram(data(ti:ti+w-1,:),[],freq,Fsdown);
                                
                                timeBins = Nsamples{cond,trial}(1)+(ti*ffact):sampStep:Nsamples{cond,trial}(1)+(ti*ffact)+...
                                    (w-1)*ffact;
                                Spkf = zeros(size(timeBins,2)-1,size(Spkchs,2));
                                
                                for step =1:length(timeBins)-1
                                    idxs = find((spdata(2,:) >= timeBins(step)) & ...
                                        (spdata(2,:) < timeBins(step+1)-1));
                                    for c=1:size(Spkchs,2)
                                        ch = Spkchs(c);
                                        Spkf(step,c) = length(find(spdata(1,idxs) == ch));
                                    end
                                end
                                % tfe = [totalfeats;[tpFeats,ones(size(tpFeats,1),1).*condLab(cond)]];
                                tfe = [tfe;Spkf(size(Spkf,1),:)];
                                if(size(Spkf,1)-(binsinthepast-1)>0)
                                    avtfe = [avtfe;mean(Spkf(size(Spkf,1)-(binsinthepast-1):size(Spkf,1),:))];
                                else
                                    avtfe = [avtfe;mean(Spkf(1:size(Spkf,1),:))];
                                end
                                
                                here = ti+w-1;
                            end
                            
                            
                        end
                        
                        % Features have been extracted and now can be concatenated to compute a
                        % single feature vector for each time point
                        rshfeats = reshape(TPfeats,size(TPfeats,3),size(TPfeats,1)*size(TPfeats,2));
                        % Finally add the labels to the trials
                        Combfeats = [Combfeats; [rshfeats,tfe,avtfe,ones(size(rshfeats,1),1).*condLab(i)]];
                        tfe = [];
                        avtfe = [];
                    end
                    
                end
                
            end
            
        end
        
    end
    
else
    
    for trial=1:size(EpochsLFP,1)
        reset(HdMain);
        
        if(~isempty(EpochsLFP{trial,1}) && ~isempty(EpochsSPK{trial,1}))
            data = EpochsLFP{trial,1}';
            data = filter(HdMain,data(:,chs));
            
            spdata = EpochsSPK{trial,1};
            data = downsample(data,factor);
            here = w;
            
            if(round((size(data,1)-w-1)/stepby)>1)
                
                TPfeats = zeros(length(chs),1);
                
                % for window=1:round((size(data,1)/stepby))-1
                for window=1:round((size(data,1)-w-1)/stepby)
                    
                    if(window == 1)
                        %TPfeats(:,:,window) = periodogram(data(1:w,:),[],freq,Fsdown);
                        TPfeats(:,window) = sum(data(1:w,:).^2)';
                        timeBins = Nsamples{trial,1}(1):sampStep:Nsamples{trial,1}(1)+w*ffact;
                        Spkf = zeros(size(timeBins,2)-1,size(Spkchs,2));
                        
                        for step =1:length(timeBins)-1
                            idxs = find((spdata(2,:) >= timeBins(step)) & ...
                                (spdata(2,:) < timeBins(step+1)-1));
                            for c=1:size(Spkchs,2)
                                ch = Spkchs(c);
                                Spkf(step,c) = length(find(spdata(1,idxs) == ch));
                            end
                        end
                        % tfe = [totalfeats;[tpFeats,ones(size(tpFeats,1),1).*condLab(cond)]];
                        tfe = [tfe;Spkf(size(Spkf,1),:)];
                        if(size(Spkf,1)-(binsinthepast-1)>0)
                            avtfe = [avtfe;mean(Spkf(size(Spkf,1)-(binsinthepast-1):size(Spkf,1),:))];
                        else
                            avtfe = [avtfe;mean(Spkf(1:size(Spkf,1),:))];
                        end
                    else
                        ti = here-stepin+1;
                        % TPfeats(:,:,window) = periodogram(data(ti:ti+w-1,:),[],freq,Fsdown);
                        TPfeats(:,window) = sum(data(ti:ti+w-1,:).^2)';
                        timeBins = Nsamples{trial,1}(1)+(ti*ffact):sampStep:Nsamples{trial,1}(1)+(ti*ffact)+...
                            (w-1)*ffact;
                        Spkf = zeros(size(timeBins,2)-1,size(Spkchs,2));
                        
                        for step =1:length(timeBins)-1
                            idxs = find((spdata(2,:) >= timeBins(step)) & ...
                                (spdata(2,:) < timeBins(step+1)-1));
                            for c=1:size(Spkchs,2)
                                ch = Spkchs(c);
                                Spkf(step,c) = length(find(spdata(1,idxs) == ch));
                            end
                        end
                        % tfe = [totalfeats;[tpFeats,ones(size(tpFeats,1),1).*condLab(cond)]];
                        tfe = [tfe;Spkf(size(Spkf,1),:)];
                        if(size(Spkf,1)-(binsinthepast-1)>0)
                            avtfe = [avtfe;mean(Spkf(size(Spkf,1)-(binsinthepast-1):size(Spkf,1),:))];
                        else
                            avtfe = [avtfe;mean(Spkf(1:size(Spkf,1),:))];
                        end
                        
                        here = ti+w-1;
                    end
                    
                    
                end
                
                % Features have been extracted and now can be concatenated to compute a
                % single feature vector for each time point
                rshfeats = reshape(TPfeats,size(TPfeats,3),size(TPfeats,1)*size(TPfeats,2));
                
                % Upsample continuous labels to match length of LFP feats
                labels = zeros(window,2);
                
                reSampFactor = floor(length(condLab.ContLabels{trial,1})/window);
                
                foo = decimate(condLab.ContLabels{trial,1}(:,1), reSampFactor);
                labels(:,1) = foo(1:window);
                foo = decimate(condLab.ContLabels{trial,1}(:,2), reSampFactor);
                labels(:,2) = foo(1:window);
                
                % Finally add the labels to the trials
                Combfeats = [Combfeats; [rshfeats,tfe,avtfe,labels]];
                tfe = [];
                avtfe = [];
            end
            
        end
        
    end
    
end


end