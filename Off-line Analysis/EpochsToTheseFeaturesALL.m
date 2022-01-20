function EpochsToTheseFeaturesALL(theseChsSPK,TWspks,theseChsLFP,TW,OL,FB,FR)

filelist = dir('Epochs');
filelist = filelist(3:end);

conditionLabels = [0,0,1,2];
Fs = 10000;

for fi=1:2:length(filelist)

 file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);

    
tp = find(contains({filelist.name},nm(1:end-4)));
LFPfeats = [];
SPKfeats = [];

for thisfile = 1:length(tp)
nfile = [filelist(tp(thisfile)).folder,filesep,filelist(tp(thisfile)).name];
%filevol = [filelist(tp(2)).folder,filesep,filelist(tp(2)).name];


    if(contains(nfile,'Volts'))
    LFPfeats = EpochsToLFPFeatures(nfile,Fs,TW,OL,theseChsLFP,FB,FR,...
        conditionLabels);
  
    LFPfeatLabs = [num2str(TW),'-',num2str(OL),'-',num2str(FB),'-',...
                '[',num2str(FR(1)),'-',num2str(FR(2)),']'];            
%             destdirectory = ['Features',filesep,nm];
%             if(exist(destdirectory,'dir') ~= 7)
%             mkdir(destdirectory);   %create the directory
%             end
      
            
    end
    
    if(contains(nfile,'Spks'))

    SPKfeats = EpochsToSPKFeatures(nfile,theseChsSPK,TWspks,conditionLabels);
    
            SPKfeatLabs = [num2str(TWspks),'-',num2str(theseChsSPK)];
%            destdirectory = ['Features',filesep,nm];
%             if(exist(destdirectory,'dir') ~= 7)
%             mkdir(destdirectory);   %create the directory
%             end
%         
    end
end

% Combine desired features and export data


  save(['Features',filesep,nm,filesep,...
                'Combined-',theseFeats,'.mat'],'feats');
            
end

              

    

end



