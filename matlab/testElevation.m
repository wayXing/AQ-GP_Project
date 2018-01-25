%% Elevation test (4D Gaussian)

clc;clear; close all;
rng(0,'twister'); % For reproducibility
n = 20;
nt = 20;
long = linspace(-112.105912,-111.795549,n)';
lat = linspace(40.700310,40.856297,n)';
t = linspace(0,24,nt)';
[LONG,LAT,T]=meshgrid(long,lat,t);
Q_LongVec = reshape(LONG,n*n*nt,1);
Q_LatVec = reshape(LAT,n*n*nt,1);
Q_TVec = reshape(T,n*n*nt,1);

data = load('../python/elevation_map/elevationMap.mat');
elevs = double(data.elevs);
endLong = data.endLong;
endLat = data.endLat;
initLong = data.initLong;
initLat = data.initLat;
gridLongs = data.gridLongs;
gridLats = data.gridLats;
minEl = min(elevs(:));

Q_ElVec = [];
for i=1:length(Q_LongVec)
  lo = Q_LongVec(i);
  la = Q_LatVec(i);
  assert(lo>=initLong && lo<=endLong,'The longitude is out of bound for elevation look-up!')
  assert(la<=initLat && la>=endLat, 'The latitude is out of bound for elevation look-up!')
  Q_ElVec = [Q_ElVec;interp2(gridLongs,gridLats,elevs,lo,la)/1000];
end
Q_El = reshape(Q_ElVec,n,n,nt);
[Q_X1Vec,Q_X2Vec] = longLat2Meter(Q_LongVec,Q_LatVec);
Q_X1Vec = Q_X1Vec/1000;
Q_X2Vec = Q_X2Vec/1000;
Q_X1 = reshape(Q_X1Vec,n,n,nt);
Q_X2 = reshape(Q_X2Vec,n,n,nt);
Q_T = T;

% endmu=exp(-(LONG.^2+LAT.^2+El.^2+T.^2)/(2*3^2));
% endmu=sin(Q_X1/2).*sin(Q_X2/2).*cos(Q_El/2).*cos(Q_T/2);
endmu=60*(1+sin(Q_X1/4).*sin(Q_X2/4).*cos(Q_T/4))-(1/20)*(Q_El-minEl);

sampleX = [1;1;n;n;randi(n,36,1)];
sampleY = [1;n;1;n;zeros(36,1)];
i=5;
while (i<=40)
  r = randi(n);
  if ~ismember(r,sampleY(sampleX==sampleX(i)))
    sampleY(i)=r;
    i=i+1;
  end
end
sample=[sampleX,sampleY];

muSample   = [];
LongVec= [];
LatVec= [];
TVec = [];
ElVec = [];
% for k=1:2:nt
%   for i=1:2:n
%     for j=1:2:n
%       muSample   = [muSample;endmu(i,j,k)];
%       LongVec= [LongVec;LONG(i,j,k)];
%       LatVec= [LatVec;LAT(i,j,k)];
%       TVec = [TVec;T(i,j,k)];
%       ElVec = [ElVec;interp2(gridLongs,gridLats,elevs,LONG(i,j,k),LAT(i,j,k))/1000];
%     end
%   end
% end
for i=1:2:nt
  for j=1:length(sample)
    muSample   = [muSample;endmu(sample(j,1),sample(j,2),i)];
    LongVec= [LongVec;LONG(sample(j,1),sample(j,2),i)];
    LatVec= [LatVec;LAT(sample(j,1),sample(j,2),i)];
    TVec = [TVec;Q_T(sample(j,1),sample(j,2),i)];
    ElVec = [ElVec;Q_El(sample(j,1),sample(j,2),i)];
  end
end
[X1Vec,X2Vec] = longLat2Meter(LongVec,LatVec);
X1Vec = X1Vec/1000;
X2Vec = X2Vec/1000;
% 
% mu = reshape(endmu,n*n*nt,1);
% X1vec=reshape(X1,n*n*nt,1);
% X2vec=reshape(X2,n*n*nt,1);
XQ = [Q_X1Vec,Q_X2Vec,Q_ElVec,Q_TVec];
X=[X1Vec,X2Vec,ElVec,TVec];
sigmaF = 8.3779;
sigmaN = 5.81;
L=[4.7273,4.7273,0.020,7.5732];
nSmplPts = length(sample)*nt/2;
K = zeros(nSmplPts);
for i=1:nSmplPts
  for j=1:nSmplPts
    K(i,j) = kerFunc(X(i,:),X(j,:),sigmaF,L) + sigmaN^2*(i==j);
  end
end

GNoise = mvnrnd(muSample,K);
Y = GNoise';
%L0=[4.1616, 4.1616, 0.1841,6.2142]
% sigmaF0 = 22.9404;
sigmaF0 = std(Y);
L0 = [4,.5, 7];
sigmaN0 = sigmaN;
optSigmaN = false;

effOpt = true;
learnRate = 0.001;
tol = 1e-5;
isARD = true;
isSpatIsot = true;
maxIt=400;
basisFnDeg =1;

[yPred,yVar] = gpRegression(X,Y,XQ,X,Y,sigmaF0,true,L0,true,sigmaN0,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,true);
YPRED = reshape(yPred,n,n,nt);
YVAR = reshape(yVar,n,n,nt);
%%
% figure;
% scatter3(X(:,1),X(:,2),Y,10,'ro','Filled');
% hold on
% YN = reshape(yPred,20,20);
% surf(LONG,LAT,YN,YN);
scrsize = get(0,'Screensize');
figure('Position',[scrsize(3)/6,10,scrsize(3)*2/3,scrsize(4)-100]);
p2=subplot(2,2,2);
pcolor(p2,Q_X1(:,:,1),Q_X2(:,:,2),Q_El(:,:,1));
colorbar(p2,'FontSize',16,'FontWeight','bold')
shading(p2,'interp');
set(p2,'FontSize',16,'FontWeight','bold');
colorbar(p2,'FontSize',16,'FontWeight','bold');
xlabel(p2,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p2,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p2,'Elevation map','FontSize',16,'FontWeight','bold')

minYP = min(min(YPRED(:)),min(Y));
maxYP = max(max(YPRED(:)),max(Y));

p1=subplot(2,2,1);
scatter(p1,X1Vec(1:40,1),X2Vec(1:40,1),30,Y(1:40,1),'Filled');
set(p1,'FontSize',16,'FontWeight','bold');
% axis(p1,[0 20.7 0 15.75])
colorbar(p1,'FontSize',16,'FontWeight','bold')
grid(p1,'on');
axis(p1,[X1Vec(1) X1Vec(2) X2Vec(1) X2Vec(3)])
title(p1,'Data set that we do regression on','FontSize',16,'FontWeight','bold')
xlabel(p1,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p1,'y(lat) [km]','FontSize',16,'FontWeight','bold');
caxis(p1,[minYP maxYP]);
txt1=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100-60)/2 200 25],'HorizontalAlignment','center',...
        'String','data t=0');

p3=subplot(2,2,3);
pcolor(p3,Q_X1(:,:,1),Q_X2(:,:,1),YPRED(:,:,1));

shading(p3,'interp');
hold(p3,'on');
scatter(p3,X1Vec(1:40,1),X2Vec(1:40,1),30,'o','filled','MarkerFaceColor','r');
set(p3,'FontSize',16,'FontWeight','bold');
colorbar(p3,'FontSize',16,'FontWeight','bold');
xlabel(p3,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p3,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p3,'predicted Y','FontSize',16,'FontWeight','bold')
caxis(p3,[minYP maxYP]);

p4=subplot(2,2,4);
pcolor(p4,Q_X1(:,:,1),Q_X2(:,:,1),sqrt(YVAR(:,:,1)));
minVP = min(YVAR(:));
maxVP = max(YVAR(:));
shading(p4,'interp');
hold(p4,'on')
scatter(p4,X1Vec(1:40,1),X2Vec(1:40,1),30,'o','filled','MarkerFaceColor','r');
set(p4,'FontSize',16,'FontWeight','bold');
colorbar(p4,'FontSize',16,'FontWeight','bold')
xlabel(p4,'x(long) [km]','FontSize',16,'FontWeight','bold');
ylabel(p4,'y(lat) [km]','FontSize',16,'FontWeight','bold');
title(p4,'uncertainty','FontSize',16,'FontWeight','bold')
caxis(p4,[minVP 30]);

txt2=uicontrol('Style','text','FontSize',16,'FontWeight','bold',...
        'Position',[1100/2 (1100+70)/2 200 25],'HorizontalAlignment','center',...
        'String','reg. t=0');

sld = uicontrol('Style', 'slider',...
  'Min',0,'Max',nt-1,'Value',0,'SliderStep',[1/(nt-1) 2/(nt-1)],...
  'Position', [(1100)/2 (1100)/2 200 30],...
  'Callback', {@updateElevTestPlot,p1,p3,p4,X1Vec(1:40,1),X2Vec(1:40,1),Q_X1(:,:,1),Q_X2(:,:,1),t,reshape(Y,40,10),YPRED,YVAR,txt1,txt2,minYP,maxYP,minVP,maxVP});
