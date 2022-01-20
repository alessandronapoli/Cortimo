function TT = SynchTrainingLabels(VRTrainingLabels)

Mtime = VRTrainingLabels{:,2};
NSPtime = VRTrainingLabels{:,3};
Mcycle = VRTrainingLabels{:,1};
Cframe = VRTrainingLabels{:,4};

% D = seconds(NSPtime);
% D = milliseconds(D);

TT = table2timetable(VRTrainingLabels,'SampleRate',50);

%% Calculate the number of C# frames for each Matlab cycle

LastCycleId = find(diff(Mcycle) > 0);
newV = zeros(size(TT{:,2}));

for ti = 2:3


for el=2:length(LastCycleId)
   
    thesevalues = TT{LastCycleId(el-1)+1:LastCycleId(el),ti};
    subthese = sort((0:length(thesevalues)-1) * 0.02,'descend')';
    % TT{LastCycleId(el-1)+1:LastCycleId(el),2} = thesevalues - subthese;
    newV(LastCycleId(el-1)+1:LastCycleId(el)) = thesevalues - subthese;
end

if(length(newV) > LastCycleId(end))
    
    iss = length(newV) - LastCycleId(end);
    thesevalues = TT{LastCycleId(end)+1:LastCycleId(end)+iss,ti};
    subthese = sort((0:length(thesevalues)-1) * 0.02,'descend')';
    % TT{LastCycleId(el-1)+1:LastCycleId(el),2} = thesevalues - subthese;
    newV(LastCycleId(end)+1:LastCycleId(end)+iss) = thesevalues - subthese;
    
end

TT{:,ti} = newV;

end


end