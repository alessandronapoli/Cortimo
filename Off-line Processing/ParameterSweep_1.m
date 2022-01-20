clearvars; close all; clc;


filelist = dir('Epochs');
filelist = filelist(3:end);

Fs = 10000;
chs = [1:16,18:20,23:30];

% Parameter 1
TimeWins = [0.3,0.5,0.8,1,2];
% Parameter 2
overlaps = [0,0.25,0.5,0.75];
% Parameter 3
freqBins = [3,5,10,20,50];
% Parameter 4
% freqRanges = {[0,1000],[0,800],[0,500],[0,250],[0,100],...
%     [2,1000],[2,800],[2,500],[2,250],[2,100]};
freqRanges = {[0,800],[0,500],[0,250],[0,100],...
    [2,800],[2,500],[2,250],[2,100]};


% Spike parameters
spkChs = [];
spkTWs = [0.05,0.1,0.2];


% Training labels to be used to train classifiers
TrainingLabels = [0,0,1,2];

    

for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    
    
% This section is for processing LFPs
    if(contains(file,'Volts'))
        for par1=1:length(TimeWins)
            TW = TimeWins(par1);

            for par2=1:length(overlaps)
                OL = overlaps(par2);
                
                for par3=1:length(freqBins)
                    FB = freqBins(par3);
                    
                    for par4=1:length(freqRanges)
                        FR = freqRanges{par4};
            
            feats = EpochsToLFPFeatures(file,Fs,TW,OL,chs,FB,FR,TrainingLabels);
            theseFeats = [num2str(TW),'-',num2str(OL),'-',num2str(FB),'-',...
                '[',num2str(FR(1)),'-',num2str(FR(2)),']'];
            destdirectory = ['Features',filesep,nm];
            if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
            end
      
            save(['Features',filesep,nm,filesep,...
                'LFP-',theseFeats,'.mat'],'feats');
          
                    end
                end
            end
        end
        
  
    
    else
        
        % This section is for processing spike data
        if(contains(file,'Spks'))
            for par1=1:length(spkTWs)
                TW = spkTWs(par1);
                
                % for par2=1:length(spkChs)
                 feats = EpochsToSPKFeatures(file,spkChs,TW,TrainingLabels);
                 if(~isempty(feats))
            theseFeats = [num2str(TW),'-',num2str(size(feats,2)-1)];
                    
              destdirectory = ['Features',filesep,nm];
            if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
            end
      
            save(['Features',filesep,nm,filesep,...
                'SPK-',theseFeats,'.mat'],'feats');
          
                % end
                 end
            end
        end
            
    end
end