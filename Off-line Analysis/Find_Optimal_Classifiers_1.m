%clear all; close all; clc;

function Find_Optimal_Classifiers_1(isCont,thisClassifier,theseFeats)

% [filename,pathname] = uigetfile('*mat','Select the classification file to be optimized');

% list_models = {'LDA','Coarse Decision Tree','Quadratic SVM','Linear SVM',...
%     'Medium Gaussian SVM','Medium KNN','Ensemble Bagged Trees'};

if(theseFeats == -1)
    
    ffolder = 'Features';
    
else
    if(theseFeats == 1)
        ffolder = 'AverageFeatures';
    end
    
end

dirlist = dir([pwd,filesep,ffolder]);
dirlist = dirlist(3:end);


for d=1:length(dirlist)
% uigetdir('Features','Select the folder with features to be classified');
%if(length(dirlist)>1)
filelist = dir([ffolder,filesep,dirlist(d).name]);
%else
%filelist = dir(['Features',filesep,dirlist.name]);
%end
if(length(filelist) > 1)
filelist = filelist(3:end);
else
end

for fi=1:length(filelist)
    
    file = [filelist(fi).folder,filesep,filelist(fi).name];
    [~,nm,~] = fileparts(file);
    tp = load(file);
    tpp = fieldnames(tp);
    feats = tp.(tpp{1});
    
    if contains(nm,'SPK')
        newS = extractBetween(nm,'-','-');
        newS = newS{1};
    else
        if(contains(nm,'LFP'))
        newS = extractBetween(nm,'-','-');
        newS = newS{1};    
        end
        
    end
    
    if(contains(nm,'Mean'))
        LC = nm(length(nm));
        newS = extractBetween(nm,'t',LC);
        newS = [newS{1},LC];
    end
    

    dt = str2double(newS) / 1000; % pass the dt value in ms
    
    [Mod, Accuracy, list_models] = Classifiers(feats,5,isCont,dt,...
        thisClassifier);
  for i=1:length(Accuracy)
        Mod{2,i} = Accuracy(i);
  end
  
  destdirectory = [pwd,filesep,'Classifiers',filesep,dirlist(d).name];
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
 destdirectory = ['Results',filesep,dirlist(d).name];
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

