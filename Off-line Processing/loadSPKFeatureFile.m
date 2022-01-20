function [SPKfeats] = loadSPKFeatureFile(filenameT)

if(exist(filenameT,'file'))
    
fid = fopen(filenameT,'rb');

training_data = fread(fid,inf,'double');

fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%* This is the data stream saved by the real-time Cortimo suite.
% 
%    handles.data.spikefeatures = [handles.data.cycles,...
%    handles.data.ts(handles.data.cycles),handles.settings.spktimewin,...
%    str2num(handles.settings.ch_to_spi),tpe,TheseSpikeWins(:)'...
%    ];

%    fwrite(handles.settings.fid2,[length(handles.data.spikefeatures)+1,...
%    handles.data.spikefeatures],'double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SPKfeats = [];
    
try
    
% Fourth position in the data stream is the number of channels
nc = training_data(5);
% This version will only load data for unsorted spikes
% sorting = {'Unsorted','Unit1','Unit2','Unit3',...
%     'Unit4','Unit5'};
sorting = {'Unsorted'};

sortrep =  repmat(sorting,1,nc);

for cch=1:nc
    for un=1:length(sorting)
    sortrep{(cch-1)*length(sorting)+un} = ['Ch_',num2str(cch),'_',sortrep{(cch-1)*length(sorting)+un}];
    end
end


stream_labels = [{'CortimoCycle','NSPts','AnalysisTimewinDur',...
    'AnalysisTS'}, sortrep];
   

% Rebuild the data matrix
step = 1;
thismatrix = [];
temp_4 = training_data;
while(~isempty(temp_4))
   datlen = temp_4(step);
   cycle = temp_4(step+1);
   NSPTime = temp_4(step+2);
   TimeWin = temp_4(step+3);
   cn = temp_4(step+4);
   thisid = step+5;
   ch = temp_4(thisid:thisid+cn-1);
   thisid = thisid+cn;
   rows = temp_4(thisid);
   thisid = thisid+1;
   tt = temp_4(thisid:thisid+rows-1);
   
   thisid= thisid+rows;
   matrix = reshape(temp_4(thisid:datlen),cn,length(sorting)*rows);
  % matrix = reshape(temp_4(thisid:datlen)',rows,6*cn);
  % nmatrix = reshape(matrix',6,2*32);
  col = length(sorting);
  mat = zeros(rows,length(sorting)*cn);
  for row=1:rows
      ttp = matrix(:,((row-1)*length(sorting))+1:((row-1)*length(sorting))+length(sorting))';
      mat(row,:) = ttp(:);
  
  end
  
  
  thismatrix = [thismatrix;[cycle*ones(rows,1),...
       ones(rows,1)*NSPTime,ones(rows,1)*TimeWin,...
       tt,mat]];
   
    temp_4 = temp_4(step+datlen:length(temp_4));
    
end

SPKfeats = array2table(thismatrix,'VariableNames',stream_labels);


catch ME
    disp('Problems loading the spike feature file');
end

end