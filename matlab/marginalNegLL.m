function [NLL,gradLL] = marginalNegLL(x,y,L,sigmaF,sigmaN,basisFnDeg,isARD,isSpatIsot)

if isARD
  assert(length(L)==size(x,2),'The number of length scales should be the same as the independent variables in ARD mode');
else
  assert(length(L)==1,'If the ARD mode is off, you just need one length scale for all the independent variables');
end

nObs = size(x,1);
nts = size(y,2);

if basisFnDeg>=0
  nIvar = size(x,2);
  nBasis = nchoosek(basisFnDeg+nIvar,basisFnDeg);
  H     = zeros(nBasis,nObs);
end

delta = eye(nObs);
K = zeros(nObs);
if ~isARD
  Lnew = ones(1,size(x,2))*L;
else
  Lnew = L;
end

for i=1:nObs
  for j=1:nObs
    K(i,j) = kerFunc(x(i,:),x(j,:),sigmaF,Lnew) + sigmaN^2 * delta(i,j);
  end
  if basisFnDeg>=0
    H(:,i) = basisTerms(x(i,:),basisFnDeg,1,[]);
  end
end
invKy=K\y;

if basisFnDeg>=0
  m =rank(H);
end

if basisFnDeg>=0
  invKH = K\H';
  A = H * invKH;
  C = invKH * (A\(H/K));
end
NLL=0;
for t=1:nts
  NLL = NLL + 0.5*(y(:,t)'*invKy(:,t) + logdet(K,'chol') + nObs*log(2*pi));
  if basisFnDeg>=0
    NLL = NLL - 0.5*( y(:,t)'*C*y(:,t) - log(det(A)) + m*log(2*pi));
  end
end

if nargout > 1 % gradient required
  nL = length(L);
  dKdL = zeros(nObs);
  dKdSF = zeros(nObs);
  
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

  dKdSN = 2*sigmaN*delta;
  gradLL = zeros(nL+2,1);
  for t=1:nts
    term = invKy(:,t)*invKy(:,t)'-inv(K);
    gradLL = gradLL - 0.5*[trace(term*dKdSF);trace(term*dKdL);trace(term*dKdSN)];
  end
end

end