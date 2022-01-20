% clear all; close all; clc;

% [filename,pathname] = uigetfile('*mat','Select the classification file to be optimized');
function [out] = Find_Optimal_Classifiers_onlinefeats(ff,saveThisClassifier,hanfeats,elbfeats)


% list_models = {'LDA','Coarse Decision Tree','Quadratic SVM','Linear SVM',...
%     'Medium Gaussian SVM','Medium KNN','Ensemble Bagged Trees'};

FEATS{1} = hanfeats;
FEATS{2} = elbfeats;
joints = {'HAN','ELB'};

c = clock;
c = c(1:end-1);
dir = num2str(fix(c));
dir= dir(find(~isspace(dir)));

for joint=1:2
    
    feats = FEATS{joint};
    [Mod, Accuracy, list_models] = Classifiers(feats,5);
    Mod(2,:) = num2cell(Accuracy');
    destdirectory = ['Classifiers',filesep,dir];
    if(exist(destdirectory,'dir') ~= 7)
    mkdir(destdirectory);   %create the directory
    end
    % DO NOT SAVE ALL THE MODELS FOR NOW
%     save([destdirectory,filesep,...
%         nm,'.mat'],'Mod');
    
    [maxAcc,mod] = max(Accuracy);
    BestMod = list_models{mod};
    
    opti(1,:) = {'TheseFeats','BestClassifier','BestAccuracy'};
%     theseFeats = [joints{joint},'-LFP-',num2str(ff.usersettings.EEGFFTWin),...
%         '-',num2str(ff.usersettings.EEGtimewin/1000),'-',...
%         num2str(ff.usersettings.Fstep),'-',...
%         '[',num2str(ff.usersettings.Frange(1)),...
%         '-',num2str(ff.usersettings.Frange(2)),']'];

 theseFeats = [joints{joint},'-LFP-',num2str(ff.usersettings.EEGFFTWin),...
        '-',num2str(1-(ff.usersettings.EEGtimewin) ./ (ff.usersettings.EEGFFTWin*1000)),'-',...
        num2str(ff.usersettings.Fstep),'-',...
        '[',num2str(ff.usersettings.Frange(1)),...
        '-',num2str(ff.usersettings.Frange(2)),']'];
    
    opti(2,:) = {theseFeats,BestMod,maxAcc};
     
     save([destdirectory,filesep,'Results-',...
            theseFeats,'.mat'],'opti');
        
    
     % Save the model to export!!!
     if(joint==2)
             
         if(saveThisClassifier < 1 )
         ff.usersettings.Classifiers.Name.Elbow = opti{2,2};
         saveCompactModel(Mod{1,mod},[destdirectory,filesep,'best',filesep,...
            theseFeats,'-CompactElb.mat'])
         %ff.usersettings.Classifiers.MdlElbow = Mod{mod};
         ff.usersettings.Classifiers.MdlElbow = loadCompactModel([destdirectory,filesep,'best',filesep,...
            theseFeats,'-CompactElb.mat']);
         ff.usersettings.Classifiers.Feats.Elbow = opti{2,1};
         
          
         else
         ff.usersettings.Classifiers.Name.Elbow = list_models{saveThisClassifier};
         saveCompactModel(Mod{1,saveThisClassifier},[destdirectory,filesep,...
             theseFeats,...
             '-CompactElb.mat']);
         %ff.usersettings.Classifiers.MdlElbow = Mod{mod};
         ff.usersettings.Classifiers.MdlElbow = loadCompactModel([destdirectory,filesep,...
             theseFeats,...
             '-CompactElb.mat']);
         ff.usersettings.Classifiers.Feats.Elbow = theseFeats;
         end
        
     
     end
     
     
     if(joint==1)
         
         if(saveThisClassifier < 1)
        ff.usersettings.Classifiers.Name.Hand = opti{2,2};
        saveCompactModel(Mod{1,mod},[destdirectory,filesep,'best',filesep,...
           theseFeats,'-CompactHan.mat']);
        ff.usersettings.Classifiers.MdlHand = loadCompactModel([destdirectory,filesep,'best',filesep,...
           theseFeats,'-CompactHan.mat']);
        ff.usersettings.Classifiers.Feats.Hand = opti{2,1};
         else
        ff.usersettings.Classifiers.Name.Hand = list_models{saveThisClassifier};
        saveCompactModel(Mod{1,mod},[destdirectory,filesep,...
            theseFeats,...
            '-CompactHan.mat']);
        ff.usersettings.Classifiers.MdlHand = loadCompactModel([destdirectory,filesep,...
            theseFeats,...
            '-CompactHan.mat']);
        ff.usersettings.Classifiers.Feats.Hand = theseFeats;
         end

     end
     
      % Save the "best" model to export!!!
      % save(Mdl,'SVMIris');

end

out = ff;

%% Clean up folders with intermediate results!
% rmdir('Epochs', 's');
% rmdir('Features','s');
