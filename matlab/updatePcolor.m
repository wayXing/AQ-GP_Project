function updatePcolor(source,callbackdata,pltPred,pltVar,px,py,x,y,yPred,yVar,times,txt,minYP,maxYP,minVP,maxVP)

newval = source.Value;
newval = round(newval);
set(source, 'Value', newval);
t=source.Value+1;

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
caxis(pltVar,[minVP 35]);
xlabel(pltVar,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(pltVar,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(pltVar,'uncertainty (PurpleAir data)','FontSize',16,'FontWeight','bold')
hold(pltVar,'off')

set(txt,'String',['reg. t=',num2str(times(t))]);

end