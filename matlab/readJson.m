function data = readJson(fileName)

fid = fopen(fileName);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
data = JSON.parse(str);

outputNames = data.results{1,1}.series{1,1}.columns;

vals = data.results{1,1}.series{1,1}.values;

long = zeros(1,length(vals));
lat = zeros(1,length(vals));
pm2_5 = zeros(1,length(vals));
for i=1:length(vals)
  if ~isempty(vals{i}{4})
    lat(i) = vals{i}{4};
    long(i) = vals{i}{5};
    pm2_5(i) = vals{i}{length(vals{i})};
    ts(i) = datetime(vals{i}{1},'InputFormat','uuuu-MM-dd''T''HH:mm:ss''Z','TimeZone','UTC');
    id{i} = vals{i}{3};
  else
    lat(i) = -1;
    long(i) = -1;
    pm2_5(i)=-1;
    id{i}=-1;
  end
end

% filter out the sensors without latitude longitude information
lat = lat(lat~=-1);
long = long(long~=-1);
pm2_5 = pm2_5(lat~=-1);
ind = find(~isnat(ts));
ts = ts(ind);
id = id(ind);


% Roughly picking the sensors inside salt lake city
SLCind = find(lat<=40.856297 & lat>=40.700310 & long<=-111.795549 & long>=-112.105912);
lat = lat(SLCind);
long = long(SLCind);
pm2_5 = pm2_5(SLCind);
ts = ts(SLCind);
id = id(SLCind);

% ind = find(lat<=55 & long<=-50);
% lat = lat(ind);
% long = long(ind);
% pm2_5 = pm2_5(ind);

% separating sensors time slices
ts_n=[];
lat_n=[];
long_n=[];
pm2_5_n=[];
while (~isempty(id) && length(id)>=size(ts_n,2))
  ts_tmp = ts();
  ts_n  = [ts_n;ts(ismember(id,id(1)))];
  lat_n   = [lat_n;lat(ismember(id,id(1)))];
  long_n  = [long_n;long(ismember(id,id(1)))];
  pm2_5_n = [pm2_5_n;pm2_5(ismember(id,id(1)))];
  lat   = lat(~ismember(id,id(1)));
  long  = long(~ismember(id,id(1)));
  pm2_5 = pm2_5(~ismember(id,id(1)));
  ts = ts(~ismember(id,id(1)));
  id = id(~ismember(id,id(1)));
end
clear data;
data{1} = lat_n;
data{2} = long_n;
data{3} = pm2_5_n;
  
end

