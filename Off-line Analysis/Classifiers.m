%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* This function takes the feature set as a matrix with the response labels
%* stored in the last data column. The kfold parameter expresses the
%* number of kfold to be used for cross-validation testing.
%* It returns the accuracy of a set of classifiers.

%* Order of the classification models
% 1) LDA;
% 2) Decision Tree Coarse;
% 3) SVM: Quadratic;
% 4) SVM: Linear;
% 5) SVM: Medium Gaussian;
% 6) KNN: Medium;
% 7) Ensembled Bagged Trees;
%

function [Mod,accuracy,list_models] = Classifiers(data,kfold,isCont,...
    dt,UseThisClassifier)

%% Classify
%load('classify_this.mat')

if isCont == 0
    positions = data(:,end);
    x = data(:,1:end-1);
    %kfold = 3;
    
    
    list_models = {'LDA','Coarse Decision Tree','Quadratic SVM','Linear SVM',...
        'Medium Gaussian SVM','Medium KNN','Ensemble Bagged Trees','Kalman','Linear'};
    
    fig = uifigure;
    
    d = uiprogressdlg(fig,'Title','Please wait Training Classifiers',...
                'Message','Training Classifiers',...
                'Indeterminate',1);
            
    switch UseThisClassifier
                   
        case 1    
            d.Message = ['Training ', list_models{UseThisClassifier}];
            
    % Linear Discriminant Analysis
    Mod{1} = fitcdiscr(x,positions,'DiscrimType', 'linear', ...
        'Gamma', 0, 'FillCoeffs', 'off');
    
            close(d);
    
        case 2  
           d.Message = ['Training ', list_models{UseThisClassifier}];

    % Decision Tree coarse
    Mod{2} = fitctree(x,positions,'MaxNumSplits',4);
    % Optmized Decision Tree
    %Mod{2} = fitctree(x,responses,'OptimizeHyperparameters','auto');
    
               close(d);
        case 3
           d.Message = ['Training ', list_models{UseThisClassifier}];
    % SVM Quadratic
    template = templateSVM(...
        'KernelFunction', 'polynomial', ...
        'PolynomialOrder', 2, ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    Mod{3} = fitcecoc(x, positions, 'Learners', template, 'Coding', 'onevsone');
    
               close(d);
    
        case 4
           d.Message = ['Training ', list_models{UseThisClassifier}];
    % SVM Linear
    %Mod{3} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Linear',...
    %    'BoxConstraint',1);
    template = templateSVM(...
        'KernelFunction', 'linear', ...
        'PolynomialOrder', [], ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    Mod{4} = fitcecoc(x, positions,'Learners', template, 'Coding', 'onevsone');
    
                   close(d);

        case 5
           d.Message = ['Training ', list_models{UseThisClassifier}];
    % SVM Medium Gaussian
    %Mod{4} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Gaussian',...
    %    'BoxConstraint',1,'KernelScale','auto');
    template = templateSVM(...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 4.4, ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    Mod{5} = fitcecoc(x, positions, 'Learners', template, 'Coding', 'onevsone');
    
               close(d);
    
        case 6
           d.Message = ['Training ', list_models{UseThisClassifier}];
    % KNN Medium
    Mod{6} = fitcknn(x,positions,'Standardized',true,'Distance','Euclidean',...
        'DistanceWeight','equal','NumNeighbors',10);
               close(d);
    
        case 7
           d.Message = ['Training ', list_models{UseThisClassifier}];

    % Ensemble Bagged Trees
    %Mod{6} = TreeBagger(30,x,responses);
    template = templateTree(...
        'MaxNumSplits', 25);
    Mod{7} = fitcensemble(x,positions,'Method', 'Bag','NumLearningCycles', 30, ...
        'Learners', template);
    
               close(d);
    
    
        otherwise
 
            d.Message = 'Training All Classifiers';
            
      d.Message = ['Training  ',list_models{1}];
      
    % Linear Discriminant Analysis
    Mod{1} = fitcdiscr(x,positions,'DiscrimType', 'linear', ...
        'Gamma', 0, 'FillCoeffs', 'off');
    
    % Decision Tree coarse
    d.Message = ['Training  ',list_models{2}];
    Mod{2} = fitctree(x,positions,'MaxNumSplits',4);
    % Optmized Decision Tree
    %Mod{2} = fitctree(x,responses,'OptimizeHyperparameters','auto');
         

    % SVM Quadratic
    d.Message = ['Training  ',list_models{3}];
    
    template = templateSVM(...
        'KernelFunction', 'polynomial', ...
        'PolynomialOrder', 2, ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    Mod{3} = fitcecoc(x, positions, 'Learners', template, 'Coding', 'onevsone');

    
    % SVM Linear
    %Mod{3} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Linear',...
    %    'BoxConstraint',1);
    template = templateSVM(...
        'KernelFunction', 'linear', ...
        'PolynomialOrder', [], ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
        d.Message = ['Training  ',list_models{4}];

    Mod{4} = fitcecoc(x, positions,'Learners', template, 'Coding', 'onevsone');

    
    % SVM Medium Gaussian
    %Mod{4} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Gaussian',...
    %    'BoxConstraint',1,'KernelScale','auto');
       d.Message = ['Training  ',list_models{5}];

    template = templateSVM(...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 4.4, ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    Mod{5} = fitcecoc(x, positions, 'Learners', template, 'Coding', 'onevsone');

    % KNN Medium
            d.Message = ['Training  ',list_models{6}];

    Mod{6} = fitcknn(x,positions,'Standardized',true,'Distance','Euclidean',...
        'DistanceWeight','equal','NumNeighbors',10);
    
    % Ensemble Bagged Trees
    %Mod{6} = TreeBagger(30,x,responses);
    template = templateTree(...
        'MaxNumSplits', 25);
            d.Message = ['Training  ',list_models{7}];

    Mod{7} = fitcensemble(x,positions,'Method', 'Bag','NumLearningCycles', 30, ...
        'Learners', template);
    
    close(d);
    
    end
    
    accuracy = zeros(length(list_models),1);
    
     d = uiprogressdlg(fig,'Title','Please wait. Evaluating Classifiers',...
                'Message','Cross-validation in Progress',...
                'Indeterminate',1);
            
    for model=1:length(Mod)
        if(~isempty(Mod{model}))
            cvmodel = crossval(Mod{model},'kfold',kfold);
            accuracy(model) = (1-kfoldLoss(cvmodel))*100;
        else
            accuracy(model) = 0;
        end
        
        %LL = kfoldPredict(cvmodel);
        
        %correct = 0;
        % for i=1:length(LL)
        %    if(LL(i)==responses(i))
        %     correct = correct + 1;
        %    end
        % end
        % % Performance
        % sc = (correct/length(LL))*100;
        
    end
    
    close(d);
    
    close(fig);
else
    positions = data(:,end-1:end);
    velocity = [[0,0];(diff(positions)./dt)];
    accel = [[0,0];(diff(velocity)./dt)];
    responses = [positions, velocity, accel];
    responses = responses';
    positions = positions';
    x = data(:,1:end-2)';
    step = 10;
    
    xx = [];
    
    for element=step:size(x,2)
      
        temp = x(:,element-step+1:element);
        
        xx = [xx,temp(:)];
        
    end
    
    list_models = {'Kalman Cont','Linear Cont'};
    
    Mod{1} = fitKalmanContinuous(responses, x, length(positions));
    
    Mod{2} = fitLinearContinuous(positions(:,step:size(positions,2)), xx);
    
    accuracy = zeros(length(Mod),1);
    
    
end


end