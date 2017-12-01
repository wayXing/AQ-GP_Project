function [yPred,yVar] = gpRegression(x,y,xQuery,x_tr,y_tr,sigmaF,optSigmaF,L,optL,sigmaN,optSigmaN,basisFnDeg,isARD,isSpatIsot,learnRate,tol,maxIt,effOpt,center)
% assert(size(y,1)>=size(y,2),'The observed values shold be in a column vector');
assert(size(x,1)>=size(x,2),['The independent variables should be in the',...
  'columns, and the observations in the rows']);
% assert(size(y_tr,1)>=size(y_tr,2),'The observed values shold be in a column vector');
assert(size(x_tr,1)>=size(x_tr,2),['The independent variables should be in the',...
  'columns, and the observations in the rows']);
assert(size(xQuery,1)>=size(xQuery,2),['The independent variables should be in the',...
  'columns, and the tests in the rows']);
nL = length(L);
if isARD
  if isSpatIsot
    assert(length(L)==size(x,2)-1,'The number of length scales should be the same as the independent variables minus 1 in spatially isotropic ARD mode');
    L = [L(1),L(1),L(2:end)];
    nL = length(L);
  else
    assert(length(L)==size(x,2),'The number of length scales should be the same as the independent variables in ARD mode');
  end
else
  assert(nL==1,'If the mode is not ARD you just need one length scale for all the independent variables');
end

% data preprocessing
nObs = size(x,1);
nQuery = size(xQuery,1);
if basisFnDeg<0
  yMean = mean(y);
  yMean_tr = mean(y_tr);
  if (center)
    y_tr=y_tr- ones(size(y_tr,1),1)*yMean_tr;
    y=y-yMean;
  end
end

if (optL || optSigmaF || optSigmaN)
  if (optL && ~optSigmaF && ~optSigmaN)
    display('Finding The optimized model parameter L...');
    logFun = @(theta) looNegLL(x_tr,y_tr,theta(1:nL),sigmaF,sigmaN,effOpt,basisFnDeg,isARD,isSpatIsot);
    theta0 = [L];
    %   options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    %   options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','GradObj','off');
    %   options = optimset('Display','iter','PlotFcns','optimplotfval');
    %   [theta,fval,exitflag,output] = fminsearch(logFun,theta0,options);
    % using gradient descent
    gradLogFun = @(theta) gradLOONLL(x_tr,y_tr,theta(1:nL),sigmaF,sigmaN,optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot);
    theta = gradDescent(gradLogFun,logFun,theta0',[true(1,nL)&optL, optSigmaF, optSigmaN],tol,learnRate,maxIt,true);
    L=abs(theta(1:nL))
  elseif (~optL && optSigmaF && ~optSigmaN)
    display('Finding The optimized model parameter sigma_F...');
    logFun = @(theta) looNegLL(x_tr,y_tr,L,theta(1),sigmaN,effOpt,isARD,basisFnDeg,isSpatIsot);
    theta0 = [sigmaF];
    %   options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    %   options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','GradObj','off');
    %   options = optimset('Display','iter','PlotFcns','optimplotfval');
    %   [theta,fval,exitflag,output] = fminsearch(logFun,theta0,options);
    % using gradient descent
    gradLogFun = @(theta) gradLOONLL(x_tr,y_tr,L',theta(1),sigmaN,optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot);
    theta = gradDescent(gradLogFun,logFun,theta0',[true(1,nL)&optL, optSigmaF, optSigmaN],tol,learnRate,maxIt,true);
    sigmaF=abs(theta(1))
  elseif (optL && optSigmaF && ~optSigmaN)
    display('Finding The optimized model parameters L and sigma_F...');
%     logFun = @(theta) looNegLL(x_tr,y_tr,theta(1:nL),theta(nL+1),sigmaN,effOpt,basisFnDeg,isARD,isSpatIsot);
    logFun = @(theta) marginalNegLL(x_tr,y_tr,theta(1:nL),theta(nL+1),sigmaN,basisFnDeg,isARD,isSpatIsot);
    theta0 = [L,sigmaF];
    % using fminunc
    %   options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','GradObj','on','PlotFcns','optimplotfval');
    %   [theta,fval,exitflag,output] = fminunc(logFun,theta0,options);
    % using fminsearch
    %   options = optimset('Display','iter','PlotFcns','optimplotfval');
    %   [theta,fval,exitflag,output] = fminsearch(logFun,theta0,options);
    % using gradient descent
%     gradLogFun = @(theta) gradLOONLL(x_tr,y_tr,theta(1:nL),theta(nL+1),sigmaN,optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot);
    gradLogFun = @(theta) gradMNLL(x,y,theta(1:nL),theta(nL+1),sigmaN,optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot);
    theta = gradDescent(gradLogFun,logFun,theta0',[true(1,nL)&optL, optSigmaF, optSigmaN],tol,learnRate,maxIt,true);
    % using SGD
%     gradLogFun = @(theta,x,y) gradLOONLL(x,y,theta(1:nL),theta(nL+1),sigmaN,optL,optSigmaF,optSigmaN,isARD,isSpatIsot);
%     gradLogFun = @(theta,x,y) gradMNLL(x,y,theta(1:nL),theta(nL+1),sigmaN,optL,optSigmaF,optSigmaN,isARD,isSpatIsot);
%     theta = SGD(gradLogFun,logFun,x_tr,y_tr,theta0',[true(1,nL)&optL, optSigmaF, optSigmaN],5,tol,learnRate,maxIt,true);
    L=abs(theta(1:nL))
    sigmaF=abs(theta(nL+1))
  elseif (optL && optSigmaF && optSigmaN)
    display('Finding The optimized model parameters L, sigma_F and sigma_N...');
    logFun = @(theta) looNegLL(x_tr,y_tr,theta(1:nL),theta(nL+1),theta(nL+2),effOpt,basisFnDeg,isARD,isSpatIsot);
    theta0 = [L,sigmaF,sigmaN];
    %   options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    %   options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','GradObj','off');
    %   options = optimset('Display','iter','PlotFcns','optimplotfval');
    %   [theta,fval,exitflag,output] = fminsearch(logFun,theta0,options);
    % using gradient descent
    gradLogFun = @(theta) gradLOONLL(x_tr,y_tr,theta(1:nL),theta(nL+1),theta(nL+2),optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot);
    theta = gradDescent(gradLogFun,logFun,theta0',[true(1,nL)&optL, optSigmaF, optSigmaN],tol,learnRate,maxIt,true);

    L=abs(theta(1:nL))
    sigmaF=abs(theta(nL+1))
    sigmaN=abs(theta(nL+2))
  else
    error('This combination is not supported for optimization');
  end
end

display('Applying the regression model...');
K = zeros(nObs);
KStar = zeros(nQuery,nObs);
Kss = zeros(nQuery,1);
delta = eye(nObs);
if basisFnDeg>=0
  nIvar = size(x,2);
  nBasis = nchoosek(basisFnDeg+nIvar,basisFnDeg);
  H     = zeros(nBasis,nObs);
  HStar = zeros(nBasis,nQuery);
end

for i=1:nObs
  for j=1:nObs
    K(i,j) = kerFunc(x(i,:),x(j,:),sigmaF,L) + sigmaN^2 * delta(i,j);
  end
  for j=1:nQuery
    KStar(j,i) = kerFunc(xQuery(j,:),x(i,:),sigmaF,L);
  end
  if basisFnDeg>=0
    H(:,i) = basisTerms(x(i,:),basisFnDeg,1,[]);
  end
end

for i=1:nQuery
  Kss(i,1) = kerFunc(xQuery(i,:),xQuery(i,:),sigmaF,L);% + sigmaN^2;
  
  if basisFnDeg>=0
    HStar(:,i) = basisTerms(xQuery(i,:),basisFnDeg,1,[]);
  end
end

invKy = K\y;
yPred = KStar * invKy;
invKKsTr= K\KStar';
if (basisFnDeg<0)
  if (center)
    yPred = yPred + yMean;
  end
else
  invKH = K\H';
  R = HStar - H * invKKsTr;
  Beta = (H*invKH) \ (H*invKy);
  yPred = yPred + R'*Beta;
  tmpTerm = (H*invKH)\R;
end


yVar = zeros(nQuery,1);
for i=1:nQuery
  yVar(i,1) = Kss(i,1) - KStar(i,:) * invKKsTr(:,i);
end
% yVar = Kss - diag(KStar * (K\KStar'));

if (basisFnDeg>=0)
  for i=1:nQuery
    yVar(i,1) = yVar(i,1) + R(:,i)' * tmpTerm(:,i);
  end
  % yVar = yVar + diag(R'*((H*invKH)\R));
end

end

