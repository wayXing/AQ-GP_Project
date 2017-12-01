function K = kerFunc(x1,x2,sigmaF,L)
M = diag(L.^(-2));
K = sigmaF^2 * exp((-(x1-x2)*M*(x1-x2)')/2);
end