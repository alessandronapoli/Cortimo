function EpochsToTheseFeaturesALL_2(theseChsSPK,TWspks,theseChsLFP,TW,OL,FB,FR,isCont,conditionLabelsHand, conditionLabelsElbow)

filelist = dir('Epochs');
filelist = filelist(3:end);

Fs = 10000;




for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    
    % tp = find(contains({filelist.name},nm(1:end-4)));
    LFPfeats = [];
    SPKfeats = [];
    Combfeats = [];
    
    if(contains(file,'Labels'))
        conditionLabels = load(file);
    end
    
    if(contains(file,'Volts'))
        
        spkfile = strrep(file,'Volts','Spks');
        if(exist(spkfile,'file') == 2)
            
            
            % Compute just LFP features
            %    LFPfeats = EpochsToLFPFeatures(file,Fs,TW,OL,theseChsLFP,FB,FR,...
            %        conditionLabels);
            
            LFPfeatLabs = [num2str(TW),'-',num2str(OL),'-',...
                '[',num2str(FR(1)),'-',num2str(FB),'-',...
                num2str(FR(2)),']','-',num2str(theseChsLFP)];
            
        while contains(LFPfeatLabs,'  ')
            LFPfeatLabs = strrep(LFPfeatLabs,'  ',' ');
        end

        LFPfeatLabs = strrep(LFPfeatLabs,' ','_');
            
            % Compute just spike features
            %    SPKfeats = EpochsToSPKFeatures(spkfile,theseChsSPK,TWspks,conditionLabels);
            if(isempty(theseChsSPK))
                theseChs = 'ALL';
            else
                theseChs = theseChsSPK;
            end
            
            SPKfeatLabs = [num2str(TWspks),'-',num2str(theseChs)];
        
        while contains(SPKfeatLabs,'  ')
            SPKfeatLabs = strrep(SPKfeatLabs,'  ',' ');
        end
        
        SPKfeatLabs = strrep(SPKfeatLabs,' ','_');
                    
            % Compute combined feautures
            if isCont == 0
                if(contains(file,'Han'))
                    Combfeats = EpochsToCMBFeatures(file,Fs,TW,OL,theseChsLFP,FB,...
                        FR,conditionLabelsHand,theseChsSPK,TWspks,isCont);
                elseif(contains(file,'Elb'))
                    Combfeats = EpochsToCMBFeatures(file,Fs,TW,OL,theseChsLFP,FB,...
                        FR,conditionLabelsElbow,theseChsSPK,TWspks,isCont);
                end
            else
                Combfeats = EpochsToCMBFeatures(file,Fs,TW,OL,theseChsLFP,FB,...
                    FR,conditionLabels,theseChsSPK,TWspks,isCont);
            end
            
            
            destdirectory = [pwd,filesep,'Features',filesep,nm];
            destdirectory = strrep(destdirectory,'Volts','Cmb');
            
            if(exist(destdirectory,'dir') ~= 7)
                mkdir(destdirectory);   %create the directory
            end
            
            nnname = [LFPfeatLabs(1:16),'---',SPKfeatLabs(1:4),'.mat'];
          % File names too long cannot be used
        if(length(nnname)>70)
           nnname = [nnname(1:70),'---'] ;
        end
        
            %      save([destdirectory,filesep,...
            %                'Combined-',[LFPfeatLabs,'---',SPKfeatLabs],'.mat'],'Combfeats');
            save([destdirectory,filesep,...
                'Combined-',nnname],'Combfeats',...
                '-v7.3');
            
            %      save(['Features',filesep,nm,filesep,...
            %                'SPK-',SPKfeatLabs,'.mat'],'SPKfeats');
            %            save(['Features',filesep,nm,filesep,...
            %                'LFP-',LFPfeatLabs,'.mat'],'LFPfeats');
        end
        
    end
    
    
    
    
end





end



