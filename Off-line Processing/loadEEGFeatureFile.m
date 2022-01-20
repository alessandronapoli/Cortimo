function [EEGfeats] = loadEEGFeatureFile(filenameT)


fid = fopen(filenameT,'rb');

training_data = fread(fid,inf,'double');

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* This is the data stream saved by the real-time Cortimo suite.
% 
%    handles.data.EEGfeatures = [handles.data.cycles, handles.data.ts(handles.data.cycles),...
%    handles.settings.EEGtimewin,length(str2num(handles.settings.ch_to_freq)),...
%    str2num(handles.settings.ch_to_freq),handles.data.Fsdown,length(updatedEEGfeatures),...
%    updatedEEGfeatures(:)'];


% Frequencies = 0:step:handles.data.Fsdown/2-1/handles.data.Fsdown/2;
% step = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Fourth position in the data stream is the number of channels
if(training_data(5)<training_data(4))
    error('OLD FILE');
end

temp_4 = training_data;
step = 1;

EEGfeats = cell(1);
count = 1;
% Rebuild the data in a cell variable with matrices containing the channels
% and the spectral power for each cycle
while(~isempty(temp_4))
  datlen = temp_4(step);
    cycle = temp_4(step+1);
    NSPTime = temp_4(step+2);
    FFTWin = temp_4(step+3);
    RefreshWin = temp_4(step+4);
    cn = temp_4(step+5);
    thisid = step+5;
    channels = temp_4(thisid+1:thisid+cn);
    thisid = thisid+cn+1;
    Fs = temp_4(thisid);
    Fstep = temp_4(thisid+1);
    Fstart = temp_4(thisid+2);
    Fstop = temp_4(thisid+3);
    thisid = thisid+4;
    freqnum = temp_4(thisid);
    thisid = thisid+1;
    EEGfeats{count,1} = {'CortimoCycles','NSPTimeStamps','FFTWin','UpdateWin','Fs'};
    EEGfeats{count,2} =  [cycle,NSPTime,FFTWin,RefreshWin,Fs];
    EEGfeats{count,3} = ['Ch Labels     ', num2str(channels')];
    EEGfeats{count,4} = ['Freq   ', num2str(Fstart:Fstep:Fstop-1)];
    EEGfeats{count,5} = reshape(temp_4(thisid:thisid+freqnum*cn-1),[freqnum,cn]);
    
    
    temp_4 = temp_4(step+datlen:length(temp_4));
    count = count +1;
    
end








% elements = 3;
% step = 1;
% thismatrix = [];
% temp_4 = training_data;
% num_channels = (temp_4(1)-4)/elements;
% while(~isempty(temp_4))
%     
%     datlen = temp_4(step);
%     cycle = temp_4(step+1);
%     AbsTime1 = temp_4(step+2);
%     AbsTime2 = temp_4(step+3);
%     num_channels = (datlen-4)/10;
% 
%    matrix = reshape(temp_4(step+4:step+datlen-1),elements,[]);
%    
%    thismatrix = [thismatrix; [cycle.*ones(size(matrix,1),1),...
%        AbsTime1.*ones(size(matrix,1),1),...
%        AbsTime2.*ones(size(matrix,1),1), matrix]];
%    
%        temp_4 = temp_4(step+datlen:length(temp_4));
%    
%         
% end
% 
% % tp = {1:num_channels'};
% % stream_labels = {'MatlabCycle','NSPTimeStart','NSPTimeStop',{1:num_channels};
% 
% % Try to build a Table variable
% % T = array2table(thismatrix,'VariableNames',stream_labels);
% 
