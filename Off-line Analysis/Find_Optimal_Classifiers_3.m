% clear all; close all; clc;

% [filename,pathname] = uigetfile('*mat','Select the classification file to be optimized');
function out = Find_Optimal_Classifiers_3(ff,isCont,dt)


% list_models = {'LDA','Coarse Decision Tree','Quadratic SVM','Linear SVM',...
%     'Medium Gaussian SVM','Medium KNN','Ensemble Bagged Trees','Kalman','Linear'};


dirlist = dir('Features');
dirlist = dirlist(3:end);

KalmanStructure = cell(length(dirlist),1);
LinearStructure = cell(length(dirlist),1);

labels = [];

saveThisClassifier = 0;

for d=1:length(dirlist)
    % uigetdir('Features','Select the folder with features to be classified');
    %if(length(dirlist)>1)
    filelist = dir([pwd,filesep,'Features',filesep,dirlist(d).name]);
    %else
    %filelist = dir(['Features',filesep,dirlist.name]);
    %end
    filelist = filelist(3:end);
    
    for fi=1:length(filelist)
        
        file = [filelist(fi).folder,filesep,filelist(fi).name];
        [~,nm,~] = fileparts(file);
        
        tp = load(file);
        tpp = fieldnames(tp);
        feats = tp.(tpp{1});
        
        if(contains(file,'Han'))
            saveThisClassifier = ff.usersettings.DiscreteModelHandIndex;
        else
            if(contains(file,'Elb'))
            saveThisClassifier = ff.usersettings.DiscreteModelElbowIndex;
            end
            
        end
        
        [Mod, Accuracy,list_models] = Classifiers(feats,5,isCont,dt,...
            saveThisClassifier);
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
    
    if isCont == 0
        
        theseFeats = nm;
        
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
            
            [temp,newmod] = max(Accuracy);
            if(temp > maxAcc)
                maxAcc = temp;
                mod = newmod;
                BestFeats = nm(1:end-8);
                BestMod = list_models{mod};
            end
            
        end
        
        % Store the optimal conditions for this trial
        destdirectory = [pwd,filesep,'Results',filesep,dirlist(d).name];
        if(exist(destdirectory,'dir') ~= 7)
            mkdir(destdirectory);   %create the directory
            if(exist([destdirectory,filesep,'best'],'dir') ~= 7)
                mkdir([destdirectory,filesep,'best']);   %create the directory
            end
        end
        opti(1,:) = {'BestFeats','BestClassifier','BestAccuracy'};
        opti(2,:) = {BestFeats,BestMod,maxAcc};
        save([destdirectory,filesep,...
            'Optima','.mat'],'opti')
        
        % Save the model to export!!!
        if(contains(destdirectory,'Elb'))
            
            if(saveThisClassifier < 1 )
                ff.usersettings.Classifiers.Name.Elbow = opti{2,2};
                saveCompactModel(Mod{1,mod},[destdirectory,filesep,'best',filesep,...
                    BestFeats,'-CompactElb.mat'])
                %ff.usersettings.Classifiers.MdlElbow = Mod{mod};
                ff.usersettings.Classifiers.MdlElbow = loadCompactModel([destdirectory,filesep,'best',filesep,...
                    BestFeats,'-CompactElb.mat']);
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
        
        
        if(contains(destdirectory,'Han'))
            
            if(saveThisClassifier < 1)
                ff.usersettings.Classifiers.Name.Hand = opti{2,2};
                saveCompactModel(Mod{1,mod},[destdirectory,filesep,'best',filesep,...
                    BestFeats,'-CompactHan.mat']);
                ff.usersettings.Classifiers.MdlHand = loadCompactModel([destdirectory,filesep,'best',filesep,...
                    BestFeats,'-CompactHan.mat']);
                ff.usersettings.Classifiers.Feats.Hand = opti{2,1};
            else
                ff.usersettings.Classifiers.Name.Hand = list_models{saveThisClassifier};
                saveCompactModel(Mod{1,saveThisClassifier},[destdirectory,filesep,...
                    theseFeats,...
                    '-CompactHan.mat']);
                ff.usersettings.Classifiers.MdlHand = loadCompactModel([destdirectory,filesep,...
                    theseFeats,...
                    '-CompactHan.mat']);
                ff.usersettings.Classifiers.Feats.Hand = theseFeats;
            end
            
        end
        
    else
        if contains(destdirectory, "Cmb")
            if contains(destdirectory, "Rest")
                labels = [labels, "RestCombined"];
            else
                labels = [labels, "Combined"];
            end
        elseif contains(destdirectory, "Spks")
            if contains(destdirectory, "Rest")
                labels = [labels, "RestSpks"];
            else
                labels = [labels, "Spks"];
            end
        elseif contains(destdirectory, "Volts")
            if contains(destdirectory, "Rest")
                labels = [labels, "RestVolts"];
            else
                labels = [labels, "Volts"];
            end
        end
        
        Kalman = Mod{1,1};
        Linear = Mod{1,2};
        
        KalmanStructure{d,1} = Kalman;
        LinearStructure{d,1} = Linear;
        

        
    end
    % Save the "best" model to export!!!
    % save(Mdl,'SVMIris');
    
end

if isCont == 1
    
    ff.usersettings.Classifiers.KalmanTable = cell2table(KalmanStructure','VariableNames',labels);
    ff.usersettings.Classifiers.LinearTable = cell2table(LinearStructure','VariableNames',labels);
    
    a = ff.usersettings.Classifiers.KalmanTable;
    b = ff.usersettings.Classifiers.LinearTable;
    
    save([destdirectory,filesep,...
        'Kalman','.mat'],'a')
    
    save([destdirectory,filesep,...
        'Linear','.mat'],'b')
    
    clear a;
    clear b;
    
end

out = ff;

%% Clean up folders with intermediate results!
rmdir('Epochs', 's');
rmdir('Features','s');
