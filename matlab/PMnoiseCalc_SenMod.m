clc;clear;close all;

%////////////////////////////////////////////////////////////////
% Calculating the standard deviation of the noise in 1003 sensors
d=xlsread('data/SenMod1003.xlsx');
%d=d(d(:,1)~=0|(d(:,2)<7),:);
PMS = d(:,2);
DMS = d(:,1);

f1003fit=fit(PMS,DMS,'exp1');
figure;
plot(f1003fit,PMS,DMS);
title('MATLAB exponential fit for PurpleAir 1003');
 
% modelfun = @(b,x)(b(1)*log(b(2)*x+b(3)));
% b0=[-50;-0.001;1];
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% opts.TolFun= 1e-12;
% opts.TolX= 1e-12;
% opts.MaxIter = 1000;
% beta = nlinfit(PMS,DMS,modelfun,b0,opts);
% figure;
% plot(PMS,DMS,'.');
% hold on
% plot(PMS,modelfun(beta,PMS),'r.');

f1003lin = @(x)(0.5431*x+1.0607);
f1003log = @(x)(-54.9149*log(-0.00765*x+0.981971));
[sx,ind]=sort(PMS);
figure;
plot(PMS,DMS,'b.',sx,f1003lin(PMS(ind)),'r-',sx,f1003log(PMS(ind)),'g-','LineWidth',2)
title('Our fits for PurpleAir 1003')

sqer = (DMS-f1003lin(PMS)).^2;
er1003lin = abs(DMS-f1003lin(PMS));
%figure;
% plot(PMS,sqer,'.')
mse = sum(sqer)/length(sqer);
varse = var(sqer);
stdEr1003lin = std(er1003lin);
display(['std of the error in PurpleAire 1003 sensors = ',num2str(stdEr1003lin)]);

%////////////////////////////////////////////////////////////////
% Calculating the standard deviation of the noise in 5003 sensors
d=xlsread('data/SenMod5003.xlsx');
%d=d(d(:,1)~=0|(d(:,2)<7),:);
PMS = d(:,2);
DMS = d(:,1);

f5003fit=fit(PMS,DMS,'exp1');
figure;
plot(f5003fit,PMS,DMS);
title('MATLAB exponential fit for PurpleAir 5003')
 
% modelfun = @(b,x)(b(1)*log(b(2)*x+b(3)));
% b0=[-50;-0.001;1];
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% opts.TolFun= 1e-12;
% opts.TolX= 1e-12;
% opts.MaxIter = 1000;
% beta = nlinfit(PMS,DMS,modelfun,b0,opts);
% figure;
% plot(PMS,DMS,'.');
% hold on
% plot(PMS,modelfun(beta,PMS),'r.');

f5003lin = @(x)(0.7778*x+2.6536);
f5003log = @(x)(-67.0241*log(-0.00985*x+0.973658));
[sx,ind]=sort(PMS);
figure;
plot(PMS,DMS,'b.',sx,f5003lin(PMS(ind)),'r-',sx,f5003log(PMS(ind)),'g-','LineWidth',2)
title('Our fits for PurpleAir 5003')

sqer = (DMS-f5003lin(PMS)).^2;
er5003lin = abs(DMS-f5003lin(PMS));
%figure;
% plot(PMS,sqer,'.')
mse = sum(sqer)/length(sqer);
varse = var(sqer);
stdEr5003lin = std(er5003lin);
display(['std of the error in PurpleAire 5003 sensors = ',num2str(stdEr5003lin)]);

%////////////////////////////////////////////////////////////////
% Calculating the standard deviation of the noise in 3003 sensors
d=xlsread('data/SenMod3003.xlsx');
%d=d(d(:,1)~=0|(d(:,2)<7),:);
PMS = d(:,2);
DMS = d(:,1);

f3003fit=fit(PMS,DMS,'exp1');
figure;
plot(f3003fit,PMS,DMS);
title('MATLAB exponential fit for airU 3003')
 
% modelfun = @(b,x)(b(1)*log(b(2)*x+b(3)));
% b0=[-50;-0.001;1];
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% opts.TolFun= 1e-12;
% opts.TolX= 1e-12;
% opts.MaxIter = 1000;
% beta = nlinfit(PMS,DMS,modelfun,b0,opts);
% figure;
% plot(PMS,DMS,'.');
% hold on
% plot(PMS,modelfun(beta,PMS),'r.');

f3003lin = @(x)(0.4528*x+3.526);
[sx,ind]=sort(PMS);
figure;
plot(PMS,DMS,'b.',sx,f3003lin(PMS(ind)),'r-','LineWidth',2)
title('Our fit for airU 3003')

sqer = (DMS-f3003lin(PMS)).^2;
er3003lin = abs(DMS-f3003lin(PMS));
%figure;
% plot(PMS,sqer,'.')
mse = sum(sqer)/length(sqer);
varse = var(sqer);
stdEr3003lin = std(er3003lin);
display(['std of the error in airU 3003 sensors       = ',num2str(stdEr3003lin)]);

display('-------------------------------------------------------');
display(['Average std of the error over all type of the sensors =',num2str(mean([stdEr1003lin,stdEr3003lin,stdEr5003lin]))]);
display(['std of the error over all type of the sensors         =',num2str(std([er1003lin;er3003lin;er5003lin]))]);

