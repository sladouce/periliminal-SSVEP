function PlotAxisAtOrigin
%PlotAxisAtOrigin Plot 2D axes through the origin

% GET TICKS
X=get(gca,'Xtick');
Y=get(gca,'Ytick');

% GET LABELS
XL=get(gca,'XtickLabel');
YL=get(gca,'YtickLabel');

% GET OFFSETS
Xoff=diff(get(gca,'XLim'))./60;
Yoff=diff(get(gca,'YLim'))./60;

% DRAW AXIS LINEs
plot(get(gca,'XLim'),[0 0],'k');
plot([556 556],get(gca,'YLim'),':k');

% Plot new ticks  


box off;
% axis square;
axis off;
set(gcf,'color','w');
