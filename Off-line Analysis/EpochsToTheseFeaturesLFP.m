function EpochsToTheseFeaturesLFP(chs,TW,OL,FB,FR,isCont,conditionLabelsHand, conditionLabelsElbow)


filelist = dir('Epochs');
filelist = filelist(3:end);

Fs = 10000;
% chs = [];

for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    
    
    
    if(contains(file,'Labels'))
        conditionLabels = load(file);
    end
    
    if(contains(file,'Volts'))
        if isCont == 0
            if(contains(file,'Han'))
                [feats,chss] = EpochsToLFPFeatures(file,Fs,TW,OL,chs,FB,FR,conditionLabelsHand,isCont);
            elseif(contains(file,'Elb'))
                [feats,chss] = EpochsToLFPFeatures(file,Fs,TW,OL,chs,FB,FR,conditionLabelsElbow,isCont);
            end
        else
            [feats,chss] = EpochsToLFPFeatures(file,Fs,TW,OL,chs,FB,FR,conditionLabels,isCont);
        end
        
        theseFeats = [num2str(TW),'-',num2str(OL),'-',...
            '[',num2str(FR(1)),'-',num2str(FB),'-',...
            num2str(FR(2)),']','-',num2str(chss)];
                
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
            'LFP-',theseFeats,'.mat'],'feats','-v7.3');
    end
    
end


