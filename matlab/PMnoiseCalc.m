d=xlsread('data/AirUFitData.xlsx');
hd=d(:,11:12);
d=d(d(:,1)~=0|(d(:,2)<7),:);
PMS = d(:,1);
DMS = d(:,2);

f1=fit(PMS,DMS,'exp1');
figure;
title('MATLAB fit')
plot(f1,d(:,1),d(:,2));

f2 = @(x)(6.698*exp(0.02758*x));
[sx,ind]=sort(PMS);
figure;
plot(sx,f2(PMS(ind)),'r-',PMS,DMS,'b.')
title('Our fit')
sqer = (DMS-f2(PMS)).^2;
er = abs(DMS-f2(PMS));
%figure;
% plot(PMS,sqer,'.')

mse = sum(sqer)/length(sqer);
var(sqer)
std(er)
