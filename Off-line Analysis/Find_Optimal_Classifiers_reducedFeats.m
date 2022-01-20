%clear all; close all; clc;

function Find_Optimal_Classifiers_reducedFeats()
% [filename,pathname] = uigetfile('*mat','Select the classification file to be optimized');



% list_models = {'LDA','Coarse Decision Tree','Quadratic SVM','Linear SVM',...
%     'Medium Gaussian SVM','Medium KNN','Ensemble Bagged Trees'};

dirlist = dir('Reduced Features');
dirlist = dirlist(3:end);


for d=1:length(dirlist)
% uigetdir('Features','Select the folder with features to be classified');
%if(length(dirlist)>1)
filelist = dir(['Reduced Features',filesep,dirlist(d).name]);
%else
%filelist = dir(['Features',filesep,dirlist.name]);
%end
filelist = filelist(3:end);

for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    ffeatsID = load(file);
    
    file2 = strrep(file,'Reduced Features','Features');
    if(exist(file2,'file')==2)
    load(file2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Select the best features using the data loaded from the reduced
    % features folder files. In addition always take the last column in
    % feats because that represents the responses.
    if(~isempty(ffeatsID.Best{1,2}) || ~isempty(ffeatsID.Best{1,3}))
        sel_feats = feats(:,[ffeatsID.Best{1,2},ffeatsID.Best{1,3},size(feats,2)]);
    [Mod, Accuracy, list_models] = Classifiers(sel_feats,5,0);
    
    
    Mod(2,:) = num2cell(Accuracy');
    destdirectory = ['Classifiers',filesep,'Reduced Features',filesep,dirlist(d).name];
    if(exist(destdirectory,'dir') ~= 7)
    mkdir(destdirectory);   %create the directory
    end
    % DO NOT SAVE ALL THE MODELS FOR NOW
%     save([destdirectory,filesep,...
%         nm,'.mat'],'Mod');
  Accuracy = Accuracy';
  save([destdirectory,filesep,...
         nm,'Accuracy.mat'],'Accuracy');
    end
    
    end
    
end

maxAcc = 0;
BestMod = [];
BestFeats = [];

newlist = dir(destdirectory);
newlist = newlist(3:end);
for fi=1:length(newlist)
   file = [newlist(fi).folder,filesep,newlist(fi).name];
    [~,nm,~] = fileparts(file);
    Accuracy = [];
    load(file);
    
    [temp,mod] = max(Accuracy);
     if(temp > maxAcc)
    maxAcc = temp;
    BestFeats = nm;
    BestMod = list_models{mod};
    end
    
end

% Store the optimal conditions for this trial
 destdirectory = ['Results',filesep,'Reduced Features',filesep,dirlist(d).name];
    if(exist(destdirectory,'dir') ~= 7)
    mkdir(destdirectory);   %create the directory
    end
    opti(1,:) = {'BestFeats','BestClassifier','BestAccuracy'};
    opti(2,:) = {BestFeats,BestMod,maxAcc};
    save([destdirectory,filesep,...
         'Optima','.mat'],'opti')
     
 if(contains(destdirectory,'Elb'))
     
   saveCompactModel(Mod{1,mod},[destdirectory,filesep,'CompactElb.mat'])
  
 end
 
 
  if(contains(destdirectory,'Han'))
      
   saveCompactModel(Mod{1,mod},[destdirectory,filesep,'CompactHan.mat'])

  end
 
end



end