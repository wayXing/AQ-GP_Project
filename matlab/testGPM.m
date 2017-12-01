%% LOO Negative Likelihood test (SIN)

clc;clear; close all;
rng(0,'twister'); % For reproducibility
n = 98;
x = linspace(-10,10,n)';
sigmaF = 5;
sigmaN = 2;
L=2;
K = zeros(n);
for i=1:n
  for j=1:n
    K(i,j) = kerFunc(x(i,:),x(j,:),sigmaF,L) + sigmaN^2*(i==j);
  end
end
endmu = 5*sin(x');
GNoise = mvnrnd(endmu,K);

y = GNoise';%+ 0.2*randn(n,1);
sigmaF0 = std(y);
L0 = 4;
sigmaN0 = sigmaN;
optSigmaN = false;
% npts=100;
% Ls = linspace(0,50,npts);
% for i=1:npts
%   looNLL(i) = looNegLL(x,y,Ls(i),sigmaF0,0.2,false);
% end
% plot(Ls,looNLL)

% yMean = mean(y);
% y = y - yMean;
% sigmaF = linspace(0,2,10);
% Ls = linspace(0,8,10);
% [l,s] = meshgrid(Ls,sigmaF);
% for i=1:10
%   for j=1:10
%     looNgLL(i,j) = looNegLL(x,y,l(i,j),s(i,j),0.2,false);
%   end
% end
% y = y + yMean;
% figure;surf(l,s,looNgLL)

isARD = false;
isSpatIsot = true;
effOpt = false;
learnRate = 1e-3;
tol = 1e-4;
basisFnDeg =1;
maxIt=500;
[yPred,yVar] = gpRegression(x,y,x,x,y,sigmaF0,true,L0,true,sigmaN0,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,true);
figure;
plot(x,y,'.',x,yPred)
ulim = yPred + 1.96*yVar;
blim = yPred - 1.96*yVar;
hold on
plot(x,ulim,'k--',x,blim,'k--')

%% LOO Negative Likelihood test (2D Gaussian)

clc;clear; close all;
rng(0,'twister'); % For reproducibility
n = 20;
x1 = linspace(-10,10,20)';
x2 = linspace(-10,10,20)';
[X1,X2]=meshgrid(x1,x2);
endmu=exp(-(X1.^2+X2.^2)/(2*3^2));
endmu = reshape(endmu,400,1);
X1vec=reshape(X1,400,1);
X2vec=reshape(X2,400,1);
X=[X1vec,X2vec];                    

sigmaF = 5;
sigmaN = 2;
L=[2,5];
K = zeros(n^2);
for i=1:n^2
  for j=1:n^2
    K(i,j) = kerFunc(X(i,:),X(j,:),sigmaF,L) + sigmaN^2*(i==j);
  end
end

GNoise = mvnrnd(endmu,K);
Y = GNoise';

sigmaF0 = std(Y);
L0 = [4];
sigmaN0 = sigmaN;
optSigmaN = false;

effOpt = true;
learnRate = 1e-4;
tol = 1e-4;
isARD = true;
isSpatIsot = true;
basisFnDeg =1;
maxIt=400;
[yPred,yVar] = gpRegression(X,Y,X,X,Y,sigmaF0,true,L0,true,sigmaN0,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,true);
figure;
scatter3(X(:,1),X(:,2),Y,10,'ro','Filled');
hold on
YN = reshape(yPred,20,20);
surf(X1,X2,YN,YN);

%% LOO Negative Likelihood test (3D Gaussian)

clc;clear; close all;
rng(0,'twister'); % For reproducibility
n = 20;
nt = 20;
x1 = linspace(-10,10,n)';
x2 = linspace(-10,10,n)';
t = linspace(0,20,nt)';
[X1,X2,T]=meshgrid(x1,x2,t);
endmu=exp(-(X1.^2+X2.^2+T.^2)/(2*3^2));
sample = [2,3;3,13;4,7;5,5;5,19;6,11;7,14;8,6;9,4;9,17;10,8;11,11;12,14;13,4;13,18;15,12;16,2;17,9;18,17;19,12];

muSample   = [];
X1vec= [];
X2vec= [];
Tvec = [];
for i=1:nt
  for j=1:length(sample)
  muSample   = [muSample;endmu(sample(j,1),sample(j,2),i)];
  X1vec= [X1vec;X1(sample(j,1),sample(j,2),i)];
  X2vec= [X2vec;X2(sample(j,1),sample(j,2),i)];
  Tvec = [Tvec;T(sample(j,1),sample(j,2),i)];
  end
end
% 
% mu = reshape(endmu,n*n*nt,1);
% X1vec=reshape(X1,n*n*nt,1);
% X2vec=reshape(X2,n*n*nt,1);
XQ = [reshape(X1,n*n*nt,1),reshape(X2,n*n*nt,1),reshape(T,n*n*nt,1)];
X=[X1vec,X2vec,Tvec];
sigmaF = 5;
sigmaN = 2;
L=[2,2,3];
K = zeros(length(sample)*nt);
for i=1:n^2
  for j=1:n^2
    K(i,j) = kerFunc(X(i,:),X(j,:),sigmaF,L) + sigmaN^2*(i==j);
  end
end

GNoise = mvnrnd(muSample,K);
Y = GNoise';

sigmaF0 = std(Y);
L0 = [4,4];
sigmaN0 = sigmaN;
optSigmaN = false;

effOpt = true;
learnRate = 0.001;
tol = 1e-4;
isARD = true;
isSpatIsot = true;
maxIt=400;
basisFnDeg =0;

[yPred,yVar] = gpRegression(X,Y,XQ,X,Y,sigmaF0,true,L0,true,sigmaN0,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,true);
figure;
scatter3(X(:,1),X(:,2),Y,10,'ro','Filled');
hold on
YN = reshape(yPred,20,20);
surf(X1,X2,YN,YN);
%% LOO Negative Likelihood test

rng(0,'twister'); % For reproducibility
n = 300;
x = linspace(-20,20,n)';
sigmaF = 3.7;
l = 5;
y = sigmaF^2*exp(-(x.^2)/(2*l^2)) + 0.2*randn(n,1);
sigmaF0 = std(y);
npts=100;
Ls = linspace(0,14,npts);
for i=1:npts
  looNLL(i) = looNegLL(x,y,Ls(i),sigmaF0,0.2);
end
figure; 
plot(Ls,looNLL,'LineWidth',2)
xlabel('L')
ylabel('LOO negative log likelihood')


L0 = 1;
sigmaN = 0.2;
optimize = true;
[yPred,yVar] = gpRegression(x,y,x,x,y,sigmaF0,true,L0,true,sigmaN,false);
figure;
yN = sigmaF^2*exp(-(x.^2)/(2*l^2));
plot(x,y,'.',x,yPred,x,yN,'--','LineWidth',2)
% hold on
% ulim = yPred + 1.96*yVar;
% blim = yPred - 1.96*yVar;
% plot(x,ulim,'k--',x,blim,'k--')
xlabel('x');
ylabel('f(x)');
legend('given data points','regressed line','true function','Location','NorthWest');
%%
clear; clc; close all;
rng(0,'twister'); % For reproducibility
n = 300;
x1=linspace(-20,20,20);
x2=linspace(-20,20,20);
[X1,X2]=meshgrid(x1,x2);
L=5;
Y=exp(-(X1.^2+X2.^2)/(2*L^2))+ 0.2*randn(20);
Y = reshape(Y,400,1);
X1vec=reshape(X1,400,1);
X2vec=reshape(X2,400,1);
X=[X1vec,X2vec];

sigmaF0 = 0.14;
npts=100;
Ls = linspace(0,14,npts);
for i=1:npts
  looNLL(i) = looNegLL([Ls(i),sigmaF0],X,Y,0.2);
end
figure; 
plot(Ls,looNLL,'LineWidth',2)
xlabel('L')
ylabel('LOO negative log likelihood')

L0 = 1;
sigmaN = 0.2;

optimize=true;
[yPred,yVar] = gpRegression(X,Y,X,sigmaF0,L0,sigmaN,optimize);
YPRED = reshape(yPred,20,20);
YVAR = reshape(yVar,20,20);
figure;surf(X1,X2,YPRED,YPRED);
hold on
scatter3(X(:,1),X(:,2),Y,10,'ro','Filled');
colorbar

YN=exp(-(X1.^2+X2.^2)/(2*L^2));

figure;surf(X1,X2,YN,YN);
axis([-20 20 -20 20 -1 1.5])

%% compare with matlab

load(fullfile(matlabroot,'examples','stats','gprdata2.mat'))

gprMdl1 = fitrgp(x,y,'KernelFunction','squaredexponential');

sigma0 = 0.2;
kparams0 = [3.5, 6.2];
gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
     'KernelParameters',kparams0,'Sigma',sigma0);

ypred1 = resubPredict(gprMdl1);
ypred2 = resubPredict(gprMdl2);

learnRate = 0.001;
tol = 1e-4;
maxIt=200;
[yPred,yVar] = gpRegression(x,y,x,x,y,...
  kparams0(2),true,kparams0(1),true,sigma0,true,learnRate,tol,maxIt,true,true);

figure();
plot(x,y,'r.');
hold on
plot(x,ypred1,'b');
plot(x,ypred2,'g');
plot(x,yPred,'c');
xlabel('x');
ylabel('y');
legend({'data','default kernel parameters',...
'kparams0 = [3.5,6.2], sigma0 = 0.2','My regression'},...
'Location','Best');
title('Impact of initial kernel parameter values');
hold off

%% compare with matlab multidimensional

load(fullfile(matlabroot,'examples','stats','gprdata.mat'))

sigma0 = std(ytrain);
sigmaF0 = sigma0;
d = size(Xtrain,2);
sigmaM0 = 10*ones(d,1);

gprMdl = fitrgp(Xtrain,ytrain,'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigmaF0],'Sigma',sigma0,'Standardize',1);

sigmaM = gprMdl.KernelInformation.KernelParameters(1:end-1,1)
sigmaF = gprMdl.KernelInformation.KernelParameters(end)
sigma  = gprMdl.Sigma

X = [Xtrain(:,1:3) Xtrain(:,6)];
sigma0 = std(ytrain);
sigmaF0 = sigma0;
d = size(X,2);
sigmaM0 = 10*ones(d,1);

gprMdl = fitrgp(X,ytrain,'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigmaF0],'Sigma',sigma0,'Standardize',1);

sigmaM = gprMdl.KernelInformation.KernelParameters(1:end-1,1)
sigmaF = gprMdl.KernelInformation.KernelParameters(end)
sigma  = gprMdl.Sigma

xtest = [Xtest(:,1:3) Xtest(:,6)];
[ypred,yvar] = predict(gprMdl,xtest);
figure;
plot(ytest,'r');
hold on;
plot(ypred,'b');
legend('True response','GPR predicted values','Location','Best');
hold off

[yPred,yVar] = gpRegression(X,ytrain,xtest,X,ytrain,...
  sigmaF0,true,10,true,sigma0,true,center);
hold on;
plot(yvar,'r');
