%clearvars; close all; clc;

function ParameterSweepAllFeats_1(chsLFP,TimeWins,overlaps,freqBins,freqRanges,spkChs,spkTWs,isCont)


filelist = dir('Epochs');
filelist = filelist(3:end);

Fs = 10000;


% Training labels to be used to train classifiers
TrainingLabels = [-1,0,1,2];

    

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
            
            [feats,chss] = EpochsToLFPFeatures(file,Fs,TW,OL,chsLFP,FB,FR,TrainingLabels,isCont);
           theseFeats = [num2str(TW),'-',num2str(OL),'-',...
                '[',num2str(FR(1)),'-',num2str(FB),'-',...
                num2str(FR(2)),']','-',num2str(chss)];
           % theseFeats = strrep(theseFeats,' ','_');
            
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
                 [feats,spk] = EpochsToSPKFeatures(file,spkChs,TW,TrainingLabels,isCont);
                 if(~isempty(feats))
            theseFeats = [num2str(TW),'-',num2str(spk)];
            
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
    
%%%%%%%%%%%%% THIS CYCLES THROUGH THE FILES IN EPOCHS %%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%  COMBINED FEATURE SWEEP  %%%%%%%%%%%%%%%%%%%%%
   
   for par1=1:length(spkTWs)
       TWspks = spkTWs(par1);
       
       for par2=1:length(TimeWins)
           TW = TimeWins(par2);
            
           for par3=1:length(overlaps)
                OL = overlaps(par3);
                
                for par4=1:length(freqBins)
                    FB = freqBins(par4);
                    
                    for par5=1:length(freqRanges)
                        FR = freqRanges{par5};   
                        %%%%% This computes and saves the features passed
                        %%%%% as arguments to the function if the
                        %%%%% corresponding Epoch files are present in the
                        %%%%% Epochs folder
                        EpochsToTheseFeaturesALL_2(spkChs,TWspks,chsLFP,TW,OL,FB,FR,isCont);
                        
                        
                    end
                    
                end
                
           end
           
       end
       
   end
   
    

    

end