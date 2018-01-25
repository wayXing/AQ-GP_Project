function updateElevTestPlot(source,callbackdata,pltData,pltPred,pltVar,px,py,x,y,times,yData,yPred,yVar,txt1,txt2,minYP,maxYP,minVP,maxVP)

newval = source.Value;
newval = round(newval);
set(source, 'Value', newval);
t=source.Value+1;

tscat = ceil(t/2);
scatter(pltData,px,py,30,yData(:,tscat),'Filled');
set(pltData,'FontSize',16,'FontWeight','bold');
axis(pltData,[px(1) px(2) py(1) py(3)])

colorbar(pltData,'FontSize',16,'FontWeight','bold');
% minPM=min(PM2p5(:));
% maxPM=max(PM2p5(:));
caxis(pltData,[minYP maxYP])
grid(pltData,'on');
title(pltData,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(pltData,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(pltData,'y(lat) [km]','FontSize',16,'FontWeight','bold');

set(txt1,'String',['data t=',num2str(times(2*tscat-1))]);

pcolor(pltPred,x,y,yPred(:,:,t));
shading(pltPred,'interp');
hold(pltPred,'on');
scatter(pltPred,px,py,30,'o','filled','MarkerFaceColor','r');
set(pltPred,'FontSize',16,'FontWeight','bold');
colorbar(pltPred,'FontSize',16,'FontWeight','bold');
caxis(pltPred,[minYP maxYP]);
xlabel(pltPred,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(pltPred,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(pltPred,'predicted PM_{2.5} (PurpleAir data)','FontSize',16,'FontWeight','bold')
hold(pltPred,'off');

% minPM=min(PM2p5(:));
% maxPM=max(PM2p5(:));
% caxis(plt,[minPM maxPM])
pcolor(pltVar,x,y,sqrt(yVar(:,:,t)));
shading(pltVar,'interp');
hold(pltVar,'on')
scatter(pltVar,px,py,30,'o','filled','MarkerFaceColor','r');
set(pltVar,'FontSize',16,'FontWeight','bold');
colorbar(pltVar,'FontSize',16,'FontWeight','bold')
caxis(pltVar,[minVP 30]);
xlabel(pltVar,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(pltVar,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(pltVar,'uncertainty (PurpleAir data)','FontSize',16,'FontWeight','bold')
hold(pltVar,'off')

set(txt2,'String',['reg. t=',num2str(times(t))]);

end