function spRasterPlot(spikeMatrix,Fs)


figure
% Transform the data for the plot
el_num = 1:64;

% ToPlot = spikeMatrix .*el_num;

% Make sure to count at most 1 spike per ms per channel and then reduce the
% data size
tstep = 0.001;
num_samps = tstep * Fs;
len = floor(length(spikeMatrix) ./ num_samps);
ToPlot2 = nan(len,length(el_num));
cycle = 1;

for sam=1:num_samps:length(spikeMatrix)-num_samps
    
    thistep = any(spikeMatrix(sam:sam+num_samps-1,:)) == 1;
    ToPlot2(cycle,thistep) = 1;
    cycle = cycle + 1;
    
end

%%%% Original Implementation keeping all the data points
% tvector = (0:length(ToPlot)-1)/Fs;
% h = stem(tvector,ToPlot,'linestyle','none');
% for ch=1:length(el_num)
% h(ch).Color = 'blue';
% h(ch).Marker = 'square';
% h(ch).MarkerFaceColor = 'blue';
% h(ch).MarkerSize = 1;
% end
%xlabel('Time (sec)');
%ylabel('Electrode Number');
%zlabel('Number of Spikes');


% Build a new time vector
time = (0:len-1) .* tstep;
% % Build the matrix to plot
% ToPlot2 = ToPlot2 .*el_num;
% 
% h = stem(time,ToPlot2,'linestyle','none');
% for ch=1:length(el_num)
% h(ch).Color = 'blue';
% h(ch).Marker = 'square';
% h(ch).MarkerFaceColor = 'blue';
% h(ch).MarkerSize = 1;
% end
% 
% 
% % Let's try plotting lines
% 
% xx = [time;time]';
% figure
% hold on
% for i=1:5
% yy = [ToPlot2(:,i),ToPlot2(:,i)+1];
% line(xx',yy');
% % plot(xx',yy')
% end
% 
% end

% Convert ToPlot2 to logical matrix
ToPlot2(isnan(ToPlot2))=0;
spikes = logical(ToPlot2)';

[~,~] = plotSpikeRaster(spikes,'PlotType','vertline',...
    'TimePerBin',num_samps/Fs);
% Adjust X axis
% xticklabels(test);
