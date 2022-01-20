% Multi-channel ring buffer
function newBuffer = cycleBuffermulti(oldBuffer, newData)
N = size(newData,1);
if N >= length(oldBuffer)
    newBuffer = newData(end-N+1:end,:);
else
    newBuffer = [oldBuffer(N+1:end,:); newData];
end
