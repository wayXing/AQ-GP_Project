function [xh,xv]=longLat2Meter(long,lat)
% converting the lat degrees to km
minLat = min(lat(:));
meanLat = mean(lat(:,1))*pi/180;
latDiff = lat-minLat;
xv = latDiff*(111132.954 - 559.822 * cos(2*meanLat) + 1.175*cos(4*meanLat));
% converting the lat degrees to km
minLong = min(long(:));
longDiff = long-minLong;
a = 6378137.0;
b = 6356752.3142;
psi = atan((b/a) * tan(lat*pi/180));
xh = longDiff .* ((pi/180)* a * cos(psi));

end