function gLOONLL = gradLOONLL(x,y,L,sigmaF,sigmaN,optL, optSigmaF, optSigmaN,isARD,isSpatIsot)

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

invK = inv(K);
invKy = K\y;

nL = length(L);
dKdSF = zeros(nObs);
dKdSN = zeros(nObs);
dKdL = zeros(nObs,nL,nObs);
Z = zeros(nObs,nObs,nL+1);

if isARD
  if isSpatIsot
    M = diag(L.^-2);
    L3 = L'.^-3;
    for i=1:nObs
      for j=1:nObs
        d = x(i,:)-x(j,:);
        if optL
          d2 = d(1:2)*d(1:2)';
          dKdL(i,1,j) = sigmaF^2 * exp(-(d*M*d')/2) * (d2*L3(1));
          dKdL(i,2,j) = dKdL(i,1,j);
          for k=3:nL
            dKdL(i,k,j) = sigmaF^2 * exp(-(d*M*d')/2) * ((d(k)^2)*L3(k));
          end
        end
        if optSigmaF
          dKdSF(i,j) = 2*sigmaF * exp(-(d*M*d')/2);
        end
        if optSigmaN
          dKdSN(i,j) = 2*sigmaN * delta(i,j);
        end
      end
    end
  else
    M = diag(L.^-2);
    L3 = L'.^-3;
    for i=1:nObs
      for j=1:nObs
        d = x(i,:)-x(j,:);
        if optL
          dKdL(i,:,j) = sigmaF^2 * exp(-(d*M*d')/2) * ((d.^2).*L3);
        end
        if optSigmaF
          dKdSF(i,j) = 2*sigmaF * exp(-(d*M*d')/2);
        end
        if optSigmaN
          dKdSN(i,j) = 2*sigmaN * delta(i,j);
        end
      end
    end
  end
else
  for i=1:nObs
    for j=1:nObs
      d = x(i,:)-x(j,:);
      d2 = d*d';
      if optL
        dKdL(i,1,j) = sigmaF^2 * exp(-d2/(2*L^2)) * (d2/L^3);
      end
      if optSigmaF
        dKdSF(i,j) = 2*sigmaF * exp(-d2/(2*L^2));
      end
      if optSigmaN
        dKdSN(i,j) = 2*sigmaN * delta(i,j);
      end
    end
  end
end

dKdL = permute(dKdL,[1 3 2]);

for i=1:nL
  Z(:,:,i) = invK * dKdL(:,:,i);
  ZinvK(:,:,i) = Z(:,:,i) * invK;
end

Z(:,:,nL+1) = K \ dKdSF;
Z(:,:,nL+2) = K \ dKdSN;
ZinvK(:,:,nL+1) = Z(:,:,nL+1) / K;
ZinvK(:,:,nL+2) = Z(:,:,nL+2) / K;

gLOONLL = zeros(nL+2,1);
for t=1:nts
  for j=1:nL
    ZAlpha(:,j) = Z(:,:,j) * invKy(:,t);
  end
  ZAlpha(:,nL+1) = Z(:,:,nL+1) * invKy(:,t);
  ZAlpha(:,nL+2) = Z(:,:,nL+2) * invKy(:,t);
  for i=1:nObs
    for j=1:nL
      gLOONLL(j) = gLOONLL(j) + (invKy(i,t)*ZAlpha(i,j) - 0.5*(1+(invKy(i,t)^2/invK(i,i)))*ZinvK(i,i,j))/invK(i,i);
    end
    gLOONLL(nL+1) = gLOONLL(nL+1) + (invKy(i,t)*ZAlpha(i,nL+1) - 0.5*(1+(invKy(i,t)^2/invK(i,i)))*ZinvK(i,i,nL+1))/invK(i,i);
    gLOONLL(nL+2) = gLOONLL(nL+2) + (invKy(i,t)*ZAlpha(i,nL+2) - 0.5*(1+(invKy(i,t)^2/invK(i,i)))*ZinvK(i,i,nL+2))/invK(i,i);
  end
end
gLOONLL=-gLOONLL;

end