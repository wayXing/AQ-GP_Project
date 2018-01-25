function elevations = longLat2Elevs(longs, lats)
data = load('../python/elevation_map/elevationMap.mat');
elevs = double(data.elevs);
endLong = data.endLong;
endLat = data.endLat;
initLong = data.initLong;
initLat = data.initLat;
gridLongs = data.gridLongs;
gridLats = data.gridLats;

elevations = zeros(length(longs),1);
for i=1:length(longs)
  long = longs(i);
  lat = lats(i);
  assert(long>=initLong && long<=endLong,'The longitude is out of bound for elevation look-up!')
  assert(lat<=initLat && lat>=endLat, 'The latitude is out of bound for elevation look-up!')
  elevations(i) = interp2(gridLongs,gridLats,elevs,long,lat)/1000;
end

end