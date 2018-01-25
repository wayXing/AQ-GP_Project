%%
clc; clear; close all;
display('Extracting the data...');
   
minLat = 40.598850;
maxLat = 40.810476;
minLong = -112.001349;
maxLong =  -111.713403;

[tr_lat,tr_long,tr_time,tr_pm25] = readQueryFile('data/queriedData.csv');
tr_lat = [tr_lat(1:30);tr_lat(32:end)];
tr_long = [tr_long(1:30);tr_long(32:end)];
tr_pm25 = [tr_pm25(:,1:30),tr_pm25(:,32:end)];

% Downsampling the data
% lat_tr  = lat_tr(:,1:6:end);
% long_tr = long_tr(:,1:6:end);
% tr_pm25  = tr_pm25(16:end,:);
% tr_time = tr_time(16:end);
tr_pm25  = tr_pm25(1:24,:);
tr_time = tr_time(1:24);
tr_pm25  = tr_pm25(1:2:end,:);
tr_time = tr_time(1:2:end);
% nt = size(PMS_tr,2);
nt  = length(tr_time);
nID = length(tr_lat);

%preprocess sensor data
display('Preprocessing the data...');
[tr_xh,tr_xv] = longLat2Meter(tr_long,tr_lat);
tr_xh=tr_xh/1000;
tr_xv=tr_xv/1000;
tr_elevs = longLat2Elevs(tr_long,tr_lat);

tr_xv  = repmat(tr_xv,1,nt);
tr_xv = tr_xv(:);
tr_xh = repmat(tr_xh,1,nt);
tr_xh = tr_xh(:);
tr_elevs  = repmat(tr_elevs,1,nt);
tr_elevs = tr_elevs(:);
tr_time = repmat(tr_time',nID,1);
tr_time = tr_time(:);

%Building the grid of the query points
nPtsMin = 20;
stepSize = min((maxLong-minLong),(maxLat-minLat))/(nPtsMin-1);
minT = min(tr_time);
maxT = max(tr_time);
t=minT:6:maxT;
longs = minLong:stepSize:maxLong;
lats = minLat:stepSize:maxLat;
nLongs = length(longs);
nLats = length(lats);
nts = length(t);

[Q_Long,Q_Lat,Q_T]=meshgrid(longs,lats,t);
[Q_xh,Q_xv] = longLat2Meter(reshape(Q_Long(:,:,1),nLongs*nLats,1),reshape(Q_Lat(:,:,1),nLongs*nLats,1));
Q_xh = Q_xh/1000;
Q_xv = Q_xv/1000;
Q_El = longLat2Elevs(reshape(Q_Long(:,:,1),nLongs*nLats,1),reshape(Q_Lat(:,:,1),nLongs*nLats,1));

Q_T =reshape(Q_T,nLongs*nLats*nts,1);
Q_xh = repmat(Q_xh,nts,1);
Q_xv = repmat(Q_xv,nts,1);
Q_El = repmat(Q_El,nts,1);

Q_X=[Q_xh,Q_xv,Q_El,Q_T];

tr_X=[tr_xh(:),tr_xv(:),tr_elevs,tr_time(:)];

clear lats longs Q_xh Q_xv tr_elevs tr_time tr_xv tr_xh 
clear stepSize nt nPtsMin nID maxT minT minLat maxLat minLong maxLong
%% Applying the regression
% L0 = [1.5, 0.05, 4];sigmaF0 = std; L= [0.1791,0.1791,3.0691,1.4822];sigmaF =   11.2996;
sigmaF0 = std(tr_pm25(:));
L0 = [1.5, 0.05, 4];
sigmaN = 4.0268;

optL       = true;
optSigmaF  = true;
optSigmaN  = false;
center     = true;
effOpt     = true;
isARD      = true;
isSpatIsot = true;
learnRate  = 1e-3;
tol        = 1e-5;
maxIt      = 400;
basisFnDeg = 1;
[yPred,yVar] = gpRegression(tr_X,tr_pm25(:),Q_X,tr_X,tr_pm25(:),...
  sigmaF0,optSigmaF,L0,optL,sigmaN,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,center);

% Plotting the results
YPRED = reshape(yPred,nLats,nLongs,nts);
YVAR = reshape(yVar,nLats,nLongs,nts);

%%
% interactive scatter plot of the data

scrsize = get(0,'Screensize');
figure('Position',[scrsize(3)/6,10,scrsize(3)*2/3,scrsize(4)-100]);
p2=subplot(2,2,2);
Q_ElMat = reshape(Q_El(1:nLats*nLongs,1),nLats,nLongs);
pcolor(p2,Q_Long(:,:,1),Q_Lat(:,:,1),Q_ElMat);
colorbar(p2,'FontSize',16,'FontWeight','bold')
shading(p2,'interp');
set(p2,'FontSize',16,'FontWeight','bold');
colorbar(p2,'FontSize',16,'FontWeight','bold');
xlabel(p2,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p2,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p2,'Elevation map','FontSize',16,'FontWeight','bold')

minYP = min(min(YPRED(:)),min(tr_pm25(:)));
maxYP = max(max(YPRED(:)),max(tr_pm25(:)));

p1=subplot(2,2,1);
scatter(p1,tr_long,tr_lat,30,tr_pm25(1,:)','Filled');
set(p1,'FontSize',16,'FontWeight','bold');
colorbar(p1,'FontSize',16,'FontWeight','bold')
grid(p1,'on');
axis(p1,[Q_Long(1,1) Q_Long(1,end) Q_Lat(1,1) Q_Lat(end,1)])
title(p1,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(p1,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p1,'y(lat) [km]','FontSize',16,'FontWeight','bold');
caxis(p1,[minYP maxYP]);
txt1=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100-60)/2 200 25],'HorizontalAlignment','center',...
        'String','data t=0');

p3=subplot(2,2,3);
pcolor(p3,Q_Long(:,:,1),Q_Lat(:,:,1),YPRED(:,:,1));

shading(p3,'interp');
hold(p3,'on');
scatter(p3,tr_long,tr_lat,30,'o','filled','MarkerFaceColor','r');
set(p3,'FontSize',16,'FontWeight','bold');
colorbar(p3,'FontSize',16,'FontWeight','bold');
xlabel(p3,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p3,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p3,'predicted Y','FontSize',16,'FontWeight','bold')
caxis(p3,[minYP maxYP]);


p4=subplot(2,2,4);
pcolor(p4,Q_Long(:,:,1),Q_Lat(:,:,1),sqrt(YVAR(:,:,1)));
minVP = min(sqrt(YVAR(:)));
maxVP = max(sqrt(YVAR(:)));
shading(p4,'interp');
hold(p4,'on')
scatter(p4,tr_long,tr_lat,30,'o','filled','MarkerFaceColor','r');
set(p4,'FontSize',16,'FontWeight','bold');
colorbar(p4,'FontSize',16,'FontWeight','bold')
xlabel(p4,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p4,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p4,'uncertainty','FontSize',16,'FontWeight','bold')
caxis(p4,[minVP maxVP]);


txt2=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100+70)/2 200 25],'HorizontalAlignment',...
        'center','String','reg. t=0');

sld = uicontrol('Style', 'slider',...
  'Min',0,'Max',nts-1,'Value',0,'SliderStep',[1/(nts-1) 2/(nts-1)],...
  'Position', [(1100)/2 (1100)/2 200 30],...
  'Callback', {@updateAllPlot,p1,p3,p4,tr_long,tr_lat,...
               Q_Long(:,:,1),Q_Lat(:,:,1),t,tr_pm25',YPRED,YVAR,...
               txt1,txt2,minYP,maxYP,minVP,maxVP});



