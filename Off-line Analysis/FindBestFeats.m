% clear vars; close all; clc;
function FindBestFeats(chs,TrainingLabels)

if(isempty(chs))
   chs = 'all';
end

filelist = dir('Features');
filelist = filelist(3:end);

destdirectory = 'Reduced Features';

 if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
 end
 
% TrainingLabels = [-1,0,1,2];

     
for fi=1:length(filelist)
    
    Session = [filelist(fi).folder,filesep,filelist(fi).name];
    %[~,nm,~] = fileparts(file);

    newfilelist = dir(Session);
    newfilelist = newfilelist(3:end);
    
    for featset=1:length(newfilelist)
        thisFeat = [newfilelist(featset).folder,filesep,...
            newfilelist(featset).name];
        load(thisFeat);
        
        
        
 if(contains(Session,'Volts'))
         
         
                 
    CVP = cvpartition(feats(:,size(feats,2)),'kfold',5);
   % BestFeats = cell(1,length(TrainingLabels)-1);
   % BestChs = cell(1,length(TrainingLabels)-1);
   % BestFreqs = cell(1,length(TrainingLabels)-1);
    Best = cell(4,length(TrainingLabels));
    
     for cond=1:length(TrainingLabels)-1
         h = zeros(CVP.NumTestSets,size(feats,2)-1);
         p = zeros(CVP.NumTestSets,size(feats,2)-1);
         
        for test=1:CVP.NumTestSets
    dataTrain = feats(CVP.training(test),1:size(feats,2)-1);
    dataLabs = feats(CVP.training(test),size(feats,2));
            
           dataTrain1 = dataTrain(dataLabs==0,:);
           dataTrain2 = dataTrain(dataLabs==cond+1,:);
           [h(test,:),p(test,:),ci,stat] = ttest2(dataTrain1,dataTrain2,'Vartype','unequal');
             
        end
        
        if(~isempty(dataTrain1) && ~isempty(dataTrain2))
            
pav = nanmean(p);
ecdf(pav);
xlabel('P value');
ylabel('CDF value');
% Sort the features by mean p-value
[vals,featureIdxSortbyP] = sort(pav,2); % sort the features
idxx = find(vals<0.001); % These are the most discrimintive features
BestFeats = featureIdxSortbyP(1:length(idxx));

% Calculate the number of frequency bins
nnn1 = extractBetween(newfilelist(featset).name,'[','-');

nnn2 = extractBetween(newfilelist(featset).name,'[',']');
nnn2 = nnn2{1};
nnn2 = extractBetween(nnn2,'-',length(nnn2));
nnn1 = str2double(nnn1{1});
nnn2 = str2double(nnn2{1});
nnnstep = extractBetween(newfilelist(featset).name,'-','-');
nnnstep = str2double(nnnstep{2});

bins = length(nnn1:nnnstep:nnn2-1);
freqs = nnn1:nnnstep:nnn2-1;
% theseChs = ceil(BestFeats{cond}/bins);
% BestChs{cond} = unique(theseChs);

ccs = (size(feats,2)-1) / bins;

[ffr,BestChs] = ind2sub([bins,ccs],BestFeats);
BestFreqs = freqs(ffr);
     
    
    Best(:,1) = {'BestFeatID';'BestChs';'BestFreq';'OrigChs'};
    Best{1,cond+1} = BestFeats;
    Best{2,cond+1} = BestChs;
    Best{3,cond+1} = BestFreqs;
    Best{4,cond+1} = chs;
    thisdir = strrep(thisFeat,'Features','Reduced Features');
        end
     end
     
    [td,~,~] = fileparts(thisdir);
     if(~exist(td,'dir'))
         mkdir(td);
     end
    save(thisdir,'Best');
    
         
         
  end
    
    
    end
    
            
    
end



end
            
            
            
            