function EMGMatrix = load_MyoEMG(filenameT)

% [fileT,pathT] = uigetfile('*.bin','Select the file containing MyoPro EMG data');
% filenameT = fullfile(pathT,fileT);


fid = fopen(filenameT,'rb');

training_data = fread(fid,inf,'double');

fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The EMG data stream is  nSx5 matrix, where nS is the number of EMG
% samples and the 5 columns are:
% 1: Hand extensors, corresponding to the MyoPro YELLOW forearm sensor;
% 2: Hand flexors, corresponding to the MyoPro GREEN forearm sensor;
% 3: Elbow extensors, corresponding to the MyoPro RED arm sensor;
% 4: Elbow flexors, corresponding to the MyoPro BLUE arm sensor;
% 5: Time stamps derived at runtime from the C++ Bluetooth Controller application;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the data stream saved by the real-time Cortimo suite.
% [length(emg_2)+4;Cortimo.data.cycles; handles.time; Cortimo.data.ts(Cortimo.data.cycles);...
%  double(emg_2)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stream_labels = {'HandExtend','HandFlex','ElbowExtend','ElbowFlex','Time'};

stream_length = length(stream_labels);

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
    temp_4 = temp_4(step+datlen:length(temp_4));

end

stream_labels = [{'MatlabCycle','MatTime','NSPTime'},stream_labels];
% Try to build a Table variable
EMGMatrix = array2table(thismatrix,'VariableNames',stream_labels');



