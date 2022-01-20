function EpochsToTheseFeaturesSPK(theseChs,TWin,isCont,conditionLabelsHand, conditionLabelsElbow,...,
                                  isLeaky,isNormalize)

filelist = dir('Epochs');
filelist = filelist(3:end);


for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    
    if(contains(file,'Labels'))
        conditionLabels = load(file);
    end
    
    if(contains(file,'Spks'))
        
        if isCont == 0
            if(contains(file,'Han'))
                [feats,Ch] = EpochsToSPKFeatures(file,theseChs,TWin,conditionLabelsHand,isCont,...
                    isLeaky,isNormalize);
            elseif(contains(file,'Elb'))
                [feats,Ch] = EpochsToSPKFeatures(file,theseChs,TWin,conditionLabelsElbow,isCont,...
                    isLeaky,isNormalize);
            end
        else
            [feats,Ch] = EpochsToSPKFeatures(file,theseChs,TWin,conditionLabels,isCont,...
                isLeaky,isNormalize);
        end
        
        theseFeats = [num2str(TWin),'-',num2str(Ch)];
        
        while contains(theseFeats,'  ')
            theseFeats = strrep(theseFeats,'  ',' ');
        end
        
        theseFeats = strrep(theseFeats,' ','_');
        
        destdirectory = [pwd,filesep,'Features',filesep,nm];
        if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
        end
        
        % File names too long cannot be used
        if(length(theseFeats)>70)
           theseFeats = [theseFeats(1:70),'---'] ;
        end
        
        save(['Features',filesep,nm,filesep,...
            'SPK-',theseFeats,'.mat'],'feats','-v7.3');
        
        
    end
    
    
    
    
end