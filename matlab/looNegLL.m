function [LNLL,grad] = looNegLL(x,y,L,sigmaF,sigmaN,efficient,isARD,isSpatIsot)

if isARD
  assert(length(L)==size(x,2),'The number of length scales should be the same as the independent variables in ARD mode');
else
  assert(length(L)==1,'If the ARD mode is off, you just need one length scale for all the independent variables');
end

nObs = length(x);
nts = size(y,2);

K = zeros(nObs);
delta = eye(nObs);
if ~isARD
  Lnew = ones(1,size(x,2))*L;
else
  Lnew = L;
end

for i=1:nObs
  for j=1:nObs
    K(i,j) = kerFunc(x(i,:),x(j,:),sigmaF,Lnew) + sigmaN^2 * delta(i,j);
  end
end

sum = 0;
if (efficient)
  invK = inv(K);
  invKy = K\y;
  for t=1:nts
    for i=1:nObs
      yiPred = y(i,t) - invKy(i,t)/invK(i,i);
      yiVar = 1/invK(i,i);
      %     (y(i,t)-yiPred)^2/(2*yiVar)
      %     log(yiVar)
      sum = sum + 0.5 * log(yiVar) + (y(i,t)-yiPred)^2/(2*yiVar) + 0.5 * log(2*pi);
    end
  end
else
  KStarns = zeros(nObs,nObs-1);
  Kss = zeros(nObs,1);
%   invK = inv(K);
%   invKy = K\y;
  for i=1:nObs
    Kns{i} = K([1:i-1,i+1:nObs],[1:i-1,i+1:nObs]);
    KStarns(i,:) = K(i,[1:i-1,i+1:nObs]);
    Kss(i)=K(i,i);
  end
  for t=1:nts
    for i=1:nObs
%       yP = y(i,t) - invKy(i,t)/invK(i,i);
%       yV = 1/invK(i,i);      
      yiPred = KStarns(i,:) * (Kns{i}\y([1:i-1,i+1:nObs],t));
      yiVar = Kss(i) - diag(KStarns(i,:) * (Kns{i}\KStarns(i,:)')); 
      sum = sum + 0.5 * log(yiVar) + (y(i,t)-yiPred)^2/(2*yiVar) + 0.5 * log(2*pi);
    end
  end
end
LNLL = sum;

if nargout > 1 % gradient required
  nL = length(L);
  dKdSF = zeros(nObs);
  dKdL = zeros(nObs,nL,nObs);
  Z = zeros(nObs,nObs,nL+1);

  if isARD
    if isSpatIsot
      M = diag(L.^-2);
      L3 = L.^-3;
      d = x(i,:)-x(j,:);
      d2 = d(1:2)*d(1:2)';
      dKdL(i,1,j) = sigmaF^2 * exp(-(d*M*d')/2) * (d2*L3(1));
      dKdL(i,2,j) = dKdL(i,1,j);
      for k=3:nL
        dKdL(i,k,j) = sigmaF^2 * exp(-(d*M*d')/2) * ((d(k)^2)*L3(k));
      end
      dKdSF(i,j) = 2*sigmaF * exp(-(d*M*d')/2);
    else
      M = diag(L.^-2);
      L3 = L.^-3;
      for i=1:nObs
        for j=1:nObs
          d = x(i,:)-x(j,:);
          dKdL(i,:,j) = sigmaF^2 * exp(-(d*M*d')/2) * ((d.^2).*L3);
          dKdSF(i,j) = 2*sigmaF * exp(-(d*M*d')/2);
        end
      end
    end
  else
    for i=1:nObs
      for j=1:nObs
        d = x(i,:)-x(j,:);
        d2 = d*d';
        dKdL(i,1,j) = sigmaF^2 * exp(-d2/(2*L^2)) * (d2/L^3);
        dKdSF(i,j) = 2*sigmaF * exp(-d2/(2*L^2));
      end
    end
  end
  dKdL = permute(dKdL,[1 3 2]);
  for i=1:nL
    Z(:,:,i) = invK * dKdL(:,:,i);
    ZinvK(:,:,i) = Z(:,:,i) * invK;
  end
  Z(:,:,nL+1) = invK * dKdSF;
  ZinvK(:,:,nL+1) = Z(:,:,nL+1) * invK;
  
  grad = zeros(nL+1,1);
  for t=1:nts
    for j=1:nL
      ZAlpha(:,j) = Z(:,:,j) * invKy(:,t);
    end
    ZAlpha(:,nL+1) = Z(:,:,nL+1) * invKy(:,t);
    for i=1:nObs
      for j=1:nL
        grad(j) = grad(j) + (invKy(i,t)*ZAlpha(i,j) - 0.5*(1+(invKy(i,t)^2/invK(i,i)))*ZinvK(i,i,j))/invK(i,i);
      end
      grad(nL+1) = grad(nL+1) + (invKy(i,t)*ZAlpha(i,nL+1) - 0.5*(1+(invKy(i,t)^2/invK(i,i)))*ZinvK(i,i,nL+1))/invK(i,i);
    end
  end
  grad=-grad;
end

end