function [ElbFeats,HanFeats] = FindFeats(RsqMatElbZ,RsqMatHanZ,maxFeats)


% Identify the top maxFeats features in terms of discrimination power between
% conditions.
ElbFeats = zeros(3,maxFeats);
HanFeats = zeros(3,maxFeats);

for cond=1:3
   [~,ElbFeats(cond,:)] = maxk(RsqMatElbZ{cond}(:),maxFeats); 
   [~,HanFeats(cond,:)] = maxk(RsqMatHanZ{cond}(:),maxFeats); 

end






end