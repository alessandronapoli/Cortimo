%clearvars; close all; clc;

function ParameterSweep_2(chsLFP,TimeWins,overlaps,freqBins,freqRanges,spkChs,spkTWs,...
    isCont,labelsforTrainingHand,labelsforTrainingElbow)


filelist = dir('Epochs');
filelist = filelist(3:end);

Fs = 10000;
    

for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
   
    if(contains(file,'Han'))
    % Training labels to be used to train classifiers
    TrainingLabels = labelsforTrainingHand;
    else
        if(contains(file,'Elb'))
        % Training labels to be used to train classifiers
        TrainingLabels = labelsforTrainingElbow;
        end
    end
    
        
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
            
            feats = EpochsToLFPFeatures(file,Fs,TW,OL,chsLFP,FB,FR,TrainingLabels,isCont);
            theseFeats = [num2str(TW),'-',num2str(OL),'-',num2str(FB),'-',...
                '[',num2str(FR(1)),'-',num2str(FR(2)),']'];
            
            
             while contains(theseFeats,'  ')
                     theseFeats = strrep(theseFeats,'  ',' ');
             end
               theseFeats = strrep(theseFeats,' ','_');
               
            destdirectory = ['Features',filesep,nm];
            if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
            end
      
                  % File names too long cannot be used
        if(length(theseFeats)>70)
           theseFeats = [theseFeats(1:70),'---'] ;
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
                 feats = EpochsToSPKFeatures(file,spkChs,TW,TrainingLabels,isCont);
                 if(~isempty(feats))
            theseFeats = [num2str(TW),'-',num2str(size(feats,2)-1)];
                    
              destdirectory = ['Features',filesep,nm];
            if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
            end
            
               while contains(theseFeats,'  ')
                     theseFeats = strrep(theseFeats,'  ',' ');
               end
               theseFeats = strrep(theseFeats,' ','_');
        
       % File names too long cannot be used
        if(length(theseFeats)>70)
           theseFeats = [theseFeats(1:70),'---'] ;
        end
        
            save(['Features',filesep,nm,filesep,...
                'SPK-',theseFeats,'.mat'],'feats');
          
                % end
                 end
            end
        end
            
    end
end


end