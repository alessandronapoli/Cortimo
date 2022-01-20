function [TrainingLabels] = loadMyoProTrainingFile(filenameT)

fid = fopen(filenameT,'rb');

training_data = fread(fid,inf,'double');

fclose(fid);



% %  This is the stream coming in from the C# unity 3D app
% % // Send UDP update
% %         // Build data stream
% %         short[] trainingStream = {(short) trial_cond, (short) target_pos[0],
% %                 (short) target_pos[1], (short) target_labels[0], (short) target_labels[1],
% %                 (short) targetSuccess[0], (short) targetSuccess[1], (short) Myomo.floatsToreply[0],
% %                 (short) Myomo.floatsToreply[1], (short) Myomo.MyomoPosition[0], (short) Myomo.MyomoPosition[1]};
% % 
% % 

% Rememeber to add the additional elements that are added by the Matlab 
% Cortimo suite at run time
% % {'AcqCycle','TimeStampStart','TimeStampStop'}
% % [length(feedbackVR)+4;handles.data.cycles;...
% %        handles.time;...
% %        handles.data.ts(handles.data.cycles);...
% %        feedbackVR],'double');

% Labels for Cortimo training data 
% stream_labels = {'Trial',...
%     'WriTarget','ElbTarget',...
%     'WriLab','ElbLab',...
%     'WriHit','ElbHit',...
%     'VRWriPos','VRElbPos',...
%     'VRMyomoWriPos','VRMyomoElbPos',...
%     'MyoProWriPos','MyoProElbPos'};

stream_labels = {'Trial',...
    'WriTarget','ElbTarget',...
    'WriLab','ElbLab',...
    'WriHit','ElbHit',...
    'VRWriPos','VRElbPos',...
    'MyoProWriPos','MyoProElbPos'};


stream_length = length(stream_labels);

% elements = floor(length(training_data)/stream_length);
% cor_len = stream_length * elements;

step = 1;
thismatrix = [];
temp_4 = training_data;

while(~isempty(temp_4))
   datlen = temp_4(step);
   cycle = temp_4(step+1);
   AbsTime1 = temp_4(step+2);
   AbsTime2 = temp_4(step+3);
   
   matrix = reshape(temp_4(step+4:step+datlen-1),stream_length,[])';
   thismatrix = [thismatrix; [cycle.*ones(size(matrix,1),1),...
       AbsTime1.*ones(size(matrix,1),1),...
       AbsTime2.*ones(size(matrix,1),1), matrix]];
   
%    thismatrix = [thismatrix; cycle.*ones(size(matrix),1),...
%        AbsTime1.*ones(size(matrix),1),...
%        AbsTime2.*ones(size(matrix),1),...
%        temp_4(step+4:step+datlen-1)];
    
    temp_4 = temp_4(step+datlen:length(temp_4));
    
end

stream_labels = [{'MatlabCycle','NSPTimeStart','NSPTimeStop'},stream_labels];

% Try to build a Table variable
 TrainingLabels = array2table(thismatrix,'VariableNames',stream_labels');

