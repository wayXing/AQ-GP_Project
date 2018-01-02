%%
clc; clear; close all;
display('Extracting the data...');

% data = readCSV(data/inversion06-01To19-01.csv');
load('data/invDataReady.mat');
PMS_tr  = data{3};
lat_tr  = data{1};
long_tr = data{2};
time_tr = data{4};

lat_tr  = lat_tr(1:32,:);
long_tr = long_tr(1:32,:);
PMS_tr  = PMS_tr(1:32,:);
time_tr = time_tr(1:32,:);

nt = size(PMS_tr,2);

for i=1:size(PMS_tr,1)
  if sum(isnan(PMS_tr(i,:)))~=0
    for j=1:nt
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

% Downsampling the data
lat_tr  = lat_tr(:,1:6:end);
long_tr = long_tr(:,1:6:end);
PMS_tr  = PMS_tr(:,1:6:end);
time_tr = time_tr(:,1:6:end);
nt = size(PMS_tr,2);

% PMS = PMS_tr(:,165);
% lat = lat_tr(:,165);
% long = long_tr(:,165);
% time = time_tr(:,165);
% 
% lat_all  = lat_tr;
% long_all = long_tr;
% PMS_all  = PMS_tr;
% time_all = time_tr;
% 
% lat_tr  = [lat_tr(:,1:164),lat_tr(:,166:end)];
% long_tr = [long_tr(:,1:164),long_tr(:,166:end)];
% PMS_tr  = [PMS_tr(:,1:164),PMS_tr(:,166:end)];
% time_tr  = [time_tr(:,1:164),time_tr(:,166:end)];

display('Preprocessing the data...');
%preprocess sensor data
% calibration fit on sensor reading
PM2p5_tr=6.698*exp(0.02758*PMS_tr);
[xh_tr,xv_tr] = longLat2Meter(long_tr,lat_tr);
xh_tr=xh_tr/1000;
xv_tr=xv_tr/1000;
%%
% interactive scatter plot of the data
scrsize = get(0,'Screensize');
figure('Position',[10,1.2*scrsize(4)/2.5,scrsize(3)-10,scrsize(4)/2.5]);
p1=subplot(1,3,1);
scatter(p1,xh_tr(:,1),xv_tr(:,1),30,PM2p5_tr(:,1),'Filled');
set(p1,'FontSize',16,'FontWeight','bold');
axis(p1,[0 20.7 0 15.75])
colorbar(p1,'FontSize',16,'FontWeight','bold')
% minPM=min(PM2p5_all(:));
% maxPM=max(PM2p5_all(:));
% caxis(p1,[minPM maxPM])
grid(p1,'on');
title(p1,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(p1,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p1,'y(lat) [km]','FontSize',16,'FontWeight','bold');
txt=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[10 35 150 25],'HorizontalAlignment','center',...
        'String','data t=0');
sld = uicontrol('Style', 'slider',...
        'Min',0,'Max',nt-1,'Value',0,'SliderStep',[1/(nt-1) 5/(nt-1)],...
        'Position', [10 10 150 20],...
        'Callback', {@updateScatter,p1,xh_tr(:,1),xv_tr(:,1),PM2p5_tr,time_tr(1,:),txt});
p2=subplot(1,3,2);
p3=subplot(1,3,3);

%%
% Buildng our grid
minX = min(xh_tr(:,1));
maxX = max(xh_tr(:,1));
% rangeX = maxX-minX;
% minX = minX -0.1*rangeX;
% maxX = maxX +0.1*rangeX;
minY = min(xv_tr(:,1));
maxY = max(xv_tr(:,1));
% rangeY = maxY-minY;
% minY = minY -0.1*rangeY;
% maxY = maxY +0.1*rangeY;
nPtsMin = 50;
stepSize = min(maxX,maxY)/(nPtsMin-1);
minT = min(time_tr(1,:));
maxT = max(time_tr(1,:));
t=minT:1/3:48;
x1 = minX:stepSize:maxX;
x2 = minY:stepSize:maxY;
nPtsX = length(x1);
nPtsY = length(x2);
nts = length(t);
[X1,X2,T]=meshgrid(x1,x2,t);
X1vec=reshape(X1,nPtsX*nPtsY*nts,1);
X2vec=reshape(X2,nPtsX*nPtsY*nts,1);
Tvec=reshape(T,nPtsX*nPtsY*nts,1);

sigmaF0 = 11.1303;%std(PM2p5_tr(:));
L0 = [1.2373, 10.2606];
Xtest=[X1vec,X2vec,Tvec];

Xtrain=[xh_tr(:),xv_tr(:),time_tr(:)];

% Applying the regression
sigmaF0 = 8.3779;%std(PM2p5_tr(:));
L0 = [4.7273,  7.5732];
% L = [2.9590, 6.1135, 1.0001];sigmaF = 25.4060; %init [1 1 1] objfinal=5228
% L = [1.4413, 4.9556, 7.4372];sigmaF = 33.3028; %init [1 1 10] objfinal=5145
% L = [1.5047, 5.0918, 6.4988];sigmaF = 31.6470; %init [1 1 7] objfinal=5138
% L = [1.5006, 5.2384, 6.4557];sigmaF = 67.1683; %init [1 1 6.5] objfinal=5138
% L = [1.5006, 5.2384, 6.4557];sigmaF = 59.9835; %init [8,8,7] objfinal=5138
% sigmaF = 7.8230; % L = [1.5006, 5.2384, 6.4557]; objfinal=5113
% SGD
% L=[1.5485, 6.5134, 7.3109];sigmaF = 9.4564 %init [1,5,6] objfinal=5108
% L=[2.8006, 5.4034, 6.7623];sigmaF = 8.6876 %init [1,5,6] objfinal=5137
% L=[1.4489, 5.0547, 6.1137];sigmaF = 9.5859 %init [1,5,6] and LR=0.0001 objfinal=5126
% Spatially isotropic 
% L = [1.1988, 1.1988, 9.8074];sigmaF = 10.3506 %init[1 6] objfinal=5124
% L = [1.2334, 1.2334, 10.2087];sigmaF = 11.0507 %init[1.1988 9.8074] objfinal=5123
% L = [1.2373, 1.2373, 10.2606];sigmaF = 11.1303 % init[1.2334, 10.2087] objfinal=5122.8
 
% with basis function means L = [4.7273, 4.7273,  7.5732]; sigmaF = 8.3779

sigmaN = 5.81;
optL = false;
optSigmaF = false;
optSigmaN = false;
center=true;
effOpt = true;
learnRate = 1e-3;
tol = 1e-5;
maxIt=400;
isARD = true;
isSpatIsot = true;
basisFnDeg = 1;
[yPred,yVar] = gpRegression(Xtrain,PM2p5_tr(:),Xtest,Xtrain,PM2p5_tr(:),...
  sigmaF0,optSigmaF,L0,optL,sigmaN,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,center);

% Plotting the results
YPRED = reshape(yPred,nPtsY,nPtsX,nts);
YVAR = reshape(yVar,nPtsY,nPtsX,nts);
%%
% subplot(1,3,2);
pcolor(p2,X1(:,:,1),X2(:,:,1),YPRED(:,:,1));
minYP = min(YPRED(:));
maxYP = max(YPRED(:));
shading(p2,'interp');
hold(p2,'on');
scatter(p2,xh_tr(:,1),xv_tr(:,1),30,'o','filled','MarkerFaceColor','r');
set(p2,'FontSize',16,'FontWeight','bold');
colorbar(p2,'FontSize',16,'FontWeight','bold');
xlabel(p2,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p2,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p2,'predicted PM_{2.5} (PurpleAir data)','FontSize',16,'FontWeight','bold')
caxis(p2,[minYP maxYP]);
% subplot(1,3,3);
pcolor(p3,X1(:,:,1),X2(:,:,1),sqrt(YVAR(:,:,1)));
minVP = min(YVAR(:));
maxVP = max(YVAR(:));
shading(p3,'interp');
hold(p3,'on')
scatter(p3,xh_tr(:,1),xv_tr(:,1),30,'o','filled','MarkerFaceColor','r');
set(p3,'FontSize',16,'FontWeight','bold');
colorbar(p3,'FontSize',16,'FontWeight','bold')
xlabel(p3,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p3,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p3,'uncertainty (PurpleAir data)','FontSize',16,'FontWeight','bold')
caxis(p3,[minVP 35]);

txt2=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[10 105 150 25],'HorizontalAlignment','center',...
        'String','reg. t=0');
sld2 = uicontrol('Style', 'slider',...
        'Min',0,'Max',nts-1,'Value',0,'SliderStep',[1/(nts-1) 5/(nts-1)],...
        'Position', [10 80 150 20],...
        'Callback', {@updatePcolor,p2,p3,xh_tr(:,1),xv_tr(:,1),X1(:,:,1),X2(:,:,1),YPRED,YVAR,t,txt2,minYP,maxYP,minVP,maxVP});

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
