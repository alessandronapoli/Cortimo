function scrollThisFig(fi)

figure(fi);
a = gca;
dx = 10;

yy = ylim;
xx = xlim;
%xlim([0 dx])

% Set appropriate axis limits and settings
set(gcf,'doublebuffer','on');
%% This avoids flickering when updating the axis
set(a,'xlim',[0 dx]);
set(a,'ylim',[min(yy) max(yy)]);

% Generate constants for use in uicontrol initialization
pos=get(a,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
%% This will create a slider which is just underneath the axis
%% but still leaves room for the axis labels above the slider
xmax=max(xx);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
%% Setting up callback string to modify XLim of axis (gca)
%% based on the position of the slider (gcbo)

% Creating Uicontrol
h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-dx);

end