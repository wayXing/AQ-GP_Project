function theta = SGD(gradFun,Fun,x,y,theta0,optInd,batchSize,tol,gamma0,maxIpoc,plotObj)

nObs = length(x);
nts = size(y,2);

theta = theta0;
err=1;
it = 0;
T=0;
obj = Fun(theta0);
display (['Objective = ',num2str(obj)] );

if plotObj
  figure;
  hold on;
  plot(it,obj,'bo','MarkerFaceColor','blue');
  drawnow;
end

while err>tol && T<=maxIpoc
  T=T+1;
  shflObs = randperm(nObs);
  shflTs = randperm(nts);
  for i=1:nts
    t=shflTs(i);
    for j=1:batchSize:nObs
      obs = shflObs(j:min(j+batchSize-1,nObs));
      it = it+1;
      gamma = gamma0;
      grad  = gradFun(theta,x(obs,:),y(obs,t));
      theta   = theta - gamma * grad(optInd);
    end
  end
  obj_n = Fun(theta);
  if plotObj
    plot(T,obj_n,'bo','MarkerFaceColor','blue');
    drawnow;
  end
  err   = abs(obj-obj_n);
  obj = obj_n;
  display (['Objective change = ',num2str(err)]);
end

end
