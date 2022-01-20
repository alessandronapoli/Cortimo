

MM = loadCompactModel(['\Results\20190822-131119-001-HanEpochVolts',...
    filesep,'CompactHan.mat']);
thistrial = 17;
len = size(feats,2);
X = feats(thistrial,1:len-1);
prediction = predict(MM,X)
real = feats(thistrial,len)

