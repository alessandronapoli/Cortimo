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

function [Mod,accuracy] = classify_code_2(data,kfold)

%% Classify
%load('classify_this.mat')
responses = data(:,end);
x = data(:,1:end-1);
%kfold = 3;

% Linear Discriminant Analysis
Mod{1} = fitcdiscr(x,responses,'DiscrimType', 'linear', ...
    'Gamma', 0, 'FillCoeffs', 'off');


% Decision Tree coarse
Mod{2} = fitctree(x,responses,'MaxNumSplits',4);
% Optmized Decision Tree
%Mod{2} = fitctree(x,responses,'OptimizeHyperparameters','auto');


% SVM Quadratic
template = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 2, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);
Mod{3} = fitcecoc(x, responses, 'Learners', template, 'Coding', 'onevsone');



% SVM Linear
%Mod{3} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Linear',...
%    'BoxConstraint',1);
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);
Mod{4} = fitcecoc(x, responses,'Learners', template, 'Coding', 'onevsone');

% SVM Medium Gaussian
%Mod{4} = fitcsvm(x,responses,'Standardize',true,'KernelFunction','Gaussian',...
%    'BoxConstraint',1,'KernelScale','auto');
template = templateSVM(...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 4.4, ...
    'BoxConstraint', 1, ...
    'Standardize', true);
Mod{5} = fitcecoc(x, responses, 'Learners', template, 'Coding', 'onevsone');

% KNN Medium
Mod{6} = fitcknn(x,responses,'Standardized',true,'Distance','Euclidean',...
    'DistanceWeight','equal','NumNeighbors',10);

% Ensemble Bagged Trees
%Mod{6} = TreeBagger(30,x,responses);
template = templateTree(...
    'MaxNumSplits', 25);
Mod{7} = fitcensemble(x,responses,'Method', 'Bag','NumLearningCycles', 30, ...
    'Learners', template);

accuracy = zeros(length(Mod),1);

for model=1:length(Mod)
    if(~isempty(Mod{model}))
    cvmodel = crossval(Mod{model},'kfold',kfold);
    accuracy(model) = (1-kfoldLoss(cvmodel))*100;
    else
        accuracy(model) = [];
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


end