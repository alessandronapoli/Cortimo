function [spkfeats,theseChs] = EpochsToSPKFeatures(file,theseChs,TWin,condLab,isCont,isLeaky,isNormalize)

addpath(genpath('NPMK'));
addpath(genpath(['eeglab',filesep,'functions']));
% addpath('eeglab');
addpath(genpath('Aux Funcs'));

tp = load(file);
tpp = fieldnames(tp);
Epochs = tp.(tpp{1});

nfile = strrep(file,'Spks','TS');
tppp = load(nfile);
ttpp = fieldnames(tppp);
Nsamples = tppp.(ttpp{1});


if isempty(theseChs)
    if(~isempty(Epochs{1,1}))
        theseChs = unique(Epochs{1,1}(1,:));
    else
        if(~isempty(Epochs{2,1}))
            theseChs = unique(Epochs{2,1}(1,:));
        else
            if(~isempty(Epochs{3,1}))
                theseChs = unique(Epochs{2,1}(1,:));
            end
        end
    end
end

Fs = 30e3;
sampleWin = TWin * Fs;
totalfeats = [];

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
        %
        % if(~isempty(Epochs{cond,1}))
            
            for trial=1:size(Epochs,2)
                if(~isempty(Epochs{cond,trial}))
                    
                    % Build the vector to be used to epoch the spikes in time bins
                    % timeBins = Epochs{cond,trial}(2,1):sampleWin:Epochs{cond,trial}(2,size(Epochs{cond,trial},2));
                    timeBins = Nsamples{cond,trial}(1):sampleWin:Nsamples{cond,trial}(2);
                    tpFeats = zeros(size(timeBins,2)-1,size(theseChs,2));
                    
                    for step =1:length(timeBins)-1
                        idxs = find((Epochs{cond,trial}(2,:) >= timeBins(step)) & ...
                            (Epochs{cond,trial}(2,:) < timeBins(step+1)-1));
                        for c=1:size(theseChs,2)
                            ch = theseChs(c);
                            tpFeats(step,c) = length(find(Epochs{cond,trial}(1,idxs) == ch));
                        end
                    end
                    totalfeats = [totalfeats;[tpFeats,ones(size(tpFeats,1),1).*condLab(i)]];
                end
            end
            
            
        % end
        
    end
    
else
    
    for trial=1:size(Epochs,1)
        if(~isempty(Epochs{trial,1}))
            
            % Build the vector to be used to epoch the spikes in time bins
            % timeBins = Epochs{cond,trial}(2,1):sampleWin:Epochs{cond,trial}(2,size(Epochs{cond,trial},2));
            timeBins = Nsamples{trial,1}(1):sampleWin:Nsamples{trial,1}(2);
            tpFeats = zeros(size(timeBins,2)-1,size(theseChs,2));
            labels = zeros(size(timeBins,2)-1,size(condLab.ContLabels{trial,1},2));
            
            
            
            reSampFactor = length(condLab.ContLabels{trial,1}) / (length(timeBins) - 1);
            xq = 1:reSampFactor:length(condLab.ContLabels{trial,1});
            
            if length(xq) < (length(timeBins) - 1)
                xq((length(timeBins) - 1)) = length(condLab.ContLabels{trial,1});
            end
            
            % resampling using 1-D interpolation
            foo = interp1(condLab.ContLabels{trial,1}(:,1), xq, 'spline');
            labels(:,1) = foo(1:(length(timeBins) - 1));
            foo = interp1(condLab.ContLabels{trial,1}(:,2), xq, 'spline');
            labels(:,2) = foo(1:(length(timeBins) - 1));
            
            for step =1:length(timeBins)-1
                idxs = find((Epochs{trial,1}(2,:) >= timeBins(step)) & ...
                    (Epochs{trial,1}(2,:) < timeBins(step+1)-1));
                
                for c=1:size(theseChs,2)
                    ch = theseChs(c);
                    tpFeats(step,c) = length(find(Epochs{trial,1}(1,idxs) == ch));
                end
            end
            totalfeats = [totalfeats;[tpFeats,labels]];
        end
    end
end

spkfeats = totalfeats;

if(isLeaky)
    if(isCont)
        sub = 2;
    else
        sub = 1;
    end
    
    newspkfeats = spkfeats;
    for ll=2:size(spkfeats,1)
    newspkfeats(ll,1:size(spkfeats,2)-sub) = calcMyLeakyIntegrator(spkfeats(ll-1,1:size(spkfeats,2)-sub),...
        spkfeats(ll,1:size(spkfeats,2)-sub),-1,[]);
    end
    
    spkfeats = newspkfeats;
end


if(isNormalize)
     if(isCont)
        sub = 2;
    else
        sub = 1;
     end
    
         newspkfeats = spkfeats;
         step = 10;
for ll=2:step:size(spkfeats,1)-step
    bs = CalcBaseline(spkfeats(ll-1:ll-1+step-1,1:size(spkfeats,2)-sub),2);
   newspkfeats(ll:ll+step-1,1:size(spkfeats,2)-sub) = spkfeats(ll:ll+step-1,1:size(spkfeats,2)-sub) ./ ...
       bs(:,1);
end

end