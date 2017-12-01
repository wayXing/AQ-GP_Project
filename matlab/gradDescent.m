function theta = gradDescent(gradFun,Fun,theta0,optInd,tol,gamma0,maxIt,plotObj)

theta = theta0;
err=1;
it = 0;

obj = Fun(theta0);
display (['Objective = ',num2str(obj)] );

if plotObj
  figure;
  hold on;
  plot(it,obj,'bo','MarkerFaceColor','blue');
  xlabel('#iterations');
  ylabel('objective value');
  drawnow;
end

while err>tol && it<=maxIt
  it = it+1;
  gamma = gamma0;
  grad  = gradFun(theta);
%   theta(1) = theta(1);
%   theta(2) = theta(2) - gamma * grad(2);
  theta = theta - gamma * grad(optInd);
  obj_n = Fun(theta);
  if plotObj
    plot(it,obj_n,'bo','MarkerFaceColor','blue');
    drawnow;
  end
  err = abs(obj-obj_n);
  obj = obj_n;
  display (['it #',num2str(it),': Objective change = ',num2str(err)] );
end

end
