%%
clc; clear; close all;
display('Extracting the data...');
% trainData = readJson('../../data/dataLast1h.JSON');
% data2D = readJson('../../data/dataLast5m_old.JSON');
% PMS = data2D{3};
% lat = data2D{1};
% long = data2D{2};
% 
% PMS_tr = trainData{3};
% lat_tr = trainData{1};
% long_tr = trainData{2};

% data = readCSV('../../data/inversion06-01To19-01.csv');
load('data/invDataReady.mat');
PMS_tr  = data{3};
lat_tr  = data{1};
long_tr = data{2};
lat_tr  = lat_tr(1:32,:);
long_tr = long_tr(1:32,:);
PMS_tr  = PMS_tr(1:32,:);

for i=1:size(PMS_tr,1)
  if sum(isnan(PMS_tr(i,:)))~=0
    for j=1:size(PMS_tr,2)
      if isnan(PMS_tr(i,j)) && j==1
        PMS_tr(i,j) = PMS_tr(i,find(~isnan(PMS_tr(i,:)),1));
      elseif isnan(PMS_tr(i,j))
        z=j+1;
        while isnan(PMS_tr(i,z))
          z=z+1;
        end
        PMS_tr(i,j) = PMS_tr(i,j-1) + (PMS_tr(i,z)-PMS_tr(i,j-1))/(z-(j-1));
      end
    end
  end
end
PMS = PMS_tr(:,165);
lat = lat_tr(:,165);
long = long_tr(:,165);
lat_tr  = [lat_tr(:,1:164),lat_tr(:,166:end)];
long_tr = [long_tr(:,1:164),long_tr(:,166:end)];
PMS_tr  = [PMS_tr(:,1:164),PMS_tr(:,166:end)];


display('Preprocessing the data...');
%preprocess 2D sensor data
% calibration fit on sensor reading
PM2p5=6.698*exp(0.02758*PMS);
[xh,xv] = longLat2Meter(long(:,1),lat(:,1));
xh=xh/1000;
xv=xv/1000;
%preprocess multi-time sensor data
% calibration fit on sensor reading
PM2p5_tr=6.698*exp(0.02758*PMS_tr);
[xh_tr,xv_tr] = longLat2Meter(long_tr(:,1),lat_tr(:,1));
xh_tr=xh_tr/1000;
xv_tr=xv_tr/1000;
% scatter plot of the 2D data
figure;subplot(1,3,1);
scatter(xh,xv,30,PM2p5,'Filled');
set(gca,'FontSize',16,'FontWeight','bold');
axis([0 20.7 0 15.75])
colorbar('FontSize',16,'FontWeight','bold')
grid on
title('Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel('x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel('y(lat) [km]','FontSize',16,'FontWeight','bold');
p2=subplot(1,3,2);
p3=subplot(1,3,3);
%%
% display('Calculating the negative Leave-One-Out log likelihood objective and ploting its curve...');
% yMean_tr = mean(PM2p5_tr);
% PM2p5_tr = PM2p5_tr - ones(size(PM2p5_tr,1),1)*yMean_tr;
% ss = [10.^(-6:1)];
% L = 1;
% for i=1:16
%   looNLL(i) = looNegLL([xh_tr,xv_tr],PM2p5_tr,L,ss(i),6.3095);
% end
% PM2p5_tr = PM2p5_tr + ones(size(PM2p5_tr,1),1)*yMean_tr;
% figure;
% plot(ss,looNLL)
% ylabel('negative LOO log likelihood');
% % xlabel('L [m]');
% xlabel('\sigma_F');

% yMean = mean(PM2p5);
% PM2p5 = PM2p5 - yMean;
% sigmaF = linspace(20,60,40);
% Ls = linspace(0,20,30);
% [l,s] = meshgrid(Ls,sigmaF);
% for i=1:40
%   for j=1:30
%     looNLL(i,j) = looNegLL([l(i,j),s(i,j)],[xh,xv],PM2p5,6.3095);
%   end
% end
% PM2p5 = PM2p5 + yMean;
% surf(l,s,looNLL)
%%
% Buildng our grid
minX = min(xh);
maxX = max(xh);
% rangeX = maxX-minX;
% minX = minX -0.1*rangeX;
% maxX = maxX +0.1*rangeX;
minY = min(xv);
maxY = max(xv);
% rangeY = maxY-minY;
% minY = minY -0.1*rangeY;
% maxY = maxY +0.1*rangeY;
nPtsMin = 50;
stepSize = min(maxX,maxY)/(nPtsMin-1);
x1 = minX:stepSize:maxX;
x2 = minY:stepSize:maxY;
nPtsX = length(x1);
nPtsY = length(x2);
[X1,X2]=meshgrid(x1,x2);
X1vec=reshape(X1,nPtsX*nPtsY,1);
X2vec=reshape(X2,nPtsX*nPtsY,1);
Xtest=[X1vec,X2vec];

% Applying the regression
sigmaF0 = std(PM2p5_tr(:));
L0 = [1,1];
sigmaN = 5.81;
optL = true;
optSigmaF = true;
optSigmaN = false;
center=true;
effOpt = true;
learnRate = 0.001;
tol = 1e-5;
maxIt=400;
isARD = true;
[yPred,yVar] = gpRegression([xh,xv],PM2p5,Xtest,[xh_tr,xv_tr],PM2p5_tr,...
  sigmaF0,optSigmaF,L0,optL,sigmaN,optSigmaN,isARD,learnRate,tol,maxIt,effOpt,center);

% Plotting the results
YPRED = reshape(yPred,nPtsY,nPtsX);
YVAR = reshape(yVar,nPtsY,nPtsX);
% subplot(1,3,2);
pcolor(p2,X1,X2,YPRED);
shading(p2,'interp');
hold(p2,'on');
scatter(p2,xh,xv,30,'o','filled','MarkerFaceColor','r');
set(p2,'FontSize',16,'FontWeight','bold');
colorbar(p2,'FontSize',16,'FontWeight','bold');
xlabel(p2,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p2,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p2,'predicted PM_{2.5} (PurpleAir data)','FontSize',16,'FontWeight','bold')

% subplot(1,3,3);
pcolor(p3,X1,X2,sqrt(YVAR));
shading(p3,'interp');
hold(p3,'on')
scatter(p3,xh,xv,30,'o','filled','MarkerFaceColor','r');
set(p3,'FontSize',16,'FontWeight','bold');
colorbar(p3,'FontSize',16,'FontWeight','bold')
xlabel(p3,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p3,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p3,'uncertainty (PurpleAir data)','FontSize',16,'FontWeight','bold')

%%
sigma0 = 3;
sigmaF0 = 1;
d = 2;
sigmaM0 = 1*ones(d,1);
gprMdl = fitrgp([xh,xv],PM2p5,'KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0; sigmaF0],'Sigma',sigma0,...
'FitMethod','exact','PredictMethod','exact');

[ypred,ysd,yci] = predict(gprMdl,Xtest);

YPRED = reshape(ypred,nPtsY,nPtsX);
YVAR = reshape(ysd.^2,nPtsY,nPtsX);
figure;surf(X1,X2,YPRED,YVAR);
hold on
scatter3(xh,xv,PM2p5,10,'ro','MarkerFaceColor','r');
colorbar
%%

% Buildng our grid
minLong = min(long);
maxLong = max(long);

minLat = min(lat);
maxLat = max(lat);


nPtsMin = 50;
stepSize = min(maxLong-minLong,maxLat-minLat)/(nPtsMin-1);
xGPS1 = linspace(minLong,maxLong,66);
xGPS2 = linspace(minLat,maxLat,50);
[XGPS1,XGPS2]=meshgrid(xGPS1,xGPS2);

figure;pcolor(X1,X2,YPRED);
figure;pcolor(XGPS1,XGPS2,YPRED);

