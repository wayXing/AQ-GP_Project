function gradNLL = gradMNLL(x,y,L,sigmaF,sigmaN,optL,optSigmaF,optSigmaN,basisFnDeg,isARD,isSpatIsot)

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
nL = length(L);
dKdSF = zeros(nObs);
dKdL = zeros(nObs,nL,nObs);
dKdSN = zeros(nObs);

for i=1:nObs
  for j=1:nObs
    K(i,j) = kerFunc(x(i,:),x(j,:),sigmaF,Lnew) + sigmaN^2 * delta(i,j);
    if isARD
      if isSpatIsot
        M = diag(L.^-2);
        L3 = L'.^-3;
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
      else
        M = diag(L.^-2);
        L3 = L'.^-3;
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
    else
      d2 = (x(i,:)-x(j,:))*(x(i,:)-x(j,:))';
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
  if basisFnDeg>=0
    H(:,i) = basisTerms(x(i,:),basisFnDeg,1,[]);
  end
end
dKdL = permute(dKdL,[1 3 2]);
invKy = K\y;

if basisFnDeg>=0
  invKH = K\H';
  A = H * invKH;
  invAH = A\H;
  HinvA = H'/A;
  invA = inv(A);
end

gradNLL = zeros(nL+2,1);
for t=1:nts
  term = invKy(:,t)*invKy(:,t)'-inv(K);
  if optL
    for j=1:nL
      gradNLL(j) = gradNLL(j) - 0.5*trace(term*dKdL(:,:,j));
      if basisFnDeg>=0
        gradNLL(j) = gradNLL(j) ...
          - basisFnGradTerms(invAH,HinvA,dKdL(:,:,j),invKy(:,t),invKH,invA);
      end
    end
  end
  if optSigmaF
    gradNLL(nL+1) = gradNLL(nL+1) - 0.5*trace(term*dKdSF);
    if basisFnDeg>=0
      gradNLL(nL+1) = gradNLL(nL+1) ...
        - basisFnGradTerms(invAH,HinvA,dKdSF,invKy(:,t),invKH,invA);
    end
  end
  if optSigmaN
    gradNLL(nL+2) = gradNLL(nL+2) - 0.5*trace(term*dKdSN);
    if basisFnDeg>=0
      gradNLL(nL+2) = gradNLL(nL+2) ...
        - basisFnGradTerms(invAH,HinvA,dKdSN,invKy(:,t),invKH,invA);
    end
  end
end
end

function bFnGrTerms = basisFnGradTerms(invAH,HinvA,dKdTheta,invKy,invKH,invA)
bFnGrTerms = 0.5*(trace(invA*invKH'*dKdTheta*invKH) ...
                  - invKy'* ...
                    (dKdTheta*invKH*invAH ...
                     - HinvA*invKH'*dKdTheta*invKH*invAH ...
                     + HinvA*invKH'*dKdTheta)*...
                    invKy);
end