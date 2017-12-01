function updateScatter(source,callbackdata,plt,x,y,PM2p5,times,txt)

newval = source.Value;
newval = round(newval);
set(source, 'Value', newval);
t=source.Value+1;
scatter(plt,x,y,30,PM2p5(:,t),'Filled');
set(plt,'FontSize',16,'FontWeight','bold');
axis(plt,[0 20.7 0 15.75])
colorbar(plt,'FontSize',16,'FontWeight','bold');
% minPM=min(PM2p5(:));
% maxPM=max(PM2p5(:));
% caxis(plt,[minPM maxPM])
grid(plt,'on');
title(plt,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(plt,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(plt,'y(lat) [km]','FontSize',16,'FontWeight','bold');
set(txt,'String',['data t=',num2str(times(t))]);

end