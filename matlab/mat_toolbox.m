n = 5000;
rng(1) % For reproducibility
xTr = linspace(0.5,2.5,n)';
xTest = linspace(0.5,2.5,1000)';
y = sin(10*pi.*xTr) ./ (2.*xTr)+(xTr-1).^4 + 1.5*rand(n,1);

sigma0 = std(y);
sigmaF0 = sigma0;
d = size(xTr,2);
sigmaM0 = 10*ones(d,1);
gprMdl = fitrgp(xTr,y,'KernelFunction','squaredexponential',...
'KernelParameters',[sigmaM0; sigmaF0],'Sigma',sigma0,...
'ActiveSetSize',100,'FitMethod','sr','PredictMethod','fic');

[ypred,ysd,yci] = predict(gprMdl,xTest);

plot(xTr,y,'r.');
hold on
plot(xTest,ypred);
plot(xTest,yci(:,1),'k--');
plot(xTest,yci(:,2),'k--');
xlabel('x');
ylabel('y');