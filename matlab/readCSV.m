function data = readCSV(fileName)

tmp=readtable(fileName);
len=size(tmp,1);
clear tmp;
fid = fopen(fileName);

tline=fgetl(fid);
outputNames = textscan(tline,'%q%q%q%q%q','Delimiter',',');

load('inversionData.mat');
% long = zeros(1,len);
% lat = zeros(1,len);
% pm2_5 = zeros(1,len);
% id = zeros(1,len);
% l=1;
% tline = fgetl(fid);
% while ischar(tline)
%   t = tline(2:21);
%   d=textscan(tline(24:end),'"%d","%f","%f","%f"');
%   
%   if ~isempty(d{2})
%     lat(l) = d{2};
%     long(l) = d{3};
%     pm2_5(l) = d{4};
%     ts(l) = datetime(t,'InputFormat','uuuu-MM-dd''T''HH:mm:ss''Z','TimeZone','UTC');
%     id(l) = d{1};
%   else
%     lat(l) = -1;
%     long(l) = -1;
%     pm2_5(l)=-1;
%     id(l)=-1;
%   end
%   if mod(l,10000)==0
%     display(l)
%   end
%   l=l+1;
%   tline = fgetl(fid);
% end

fclose(fid);

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

% separating sensors time slices

nID = length(unique(id));
for i=1:nID
  ts_n{i}  = ts(ismember(id,id(1)));
  lat_n{i}   = lat(ismember(id,id(1)));
  long_n{i}  = long(ismember(id,id(1)));
  pm2_5_n{i} = pm2_5(ismember(id,id(1)));
  lat   = lat(~ismember(id,id(1)));
  long  = long(~ismember(id,id(1)));
  pm2_5 = pm2_5(~ismember(id,id(1)));
  ts = ts(~ismember(id,id(1)));
  id = id(~ismember(id,id(1)));
end

clear data;
freq = hours(1);
avgPeriod = minutes(10);
data = timeSampler(ts_n,lat_n,long_n,pm2_5_n,freq,avgPeriod);

% ts_n=[];
% lat_n=[];
% long_n=[];
% pm2_5_n=[];
% while (~isempty(id) && length(id)>=size(ts_n,2))
%   ts_tmp = ts();
%   ts_n  = [ts_n;ts(ismember(id,id(1)))];
%   lat_n   = [lat_n;lat(ismember(id,id(1)))];
%   long_n  = [long_n;long(ismember(id,id(1)))];
%   pm2_5_n = [pm2_5_n;pm2_5(ismember(id,id(1)))];
%   lat   = lat(~ismember(id,id(1)));
%   long  = long(~ismember(id,id(1)));
%   pm2_5 = pm2_5(~ismember(id,id(1)));
%   ts = ts(~ismember(id,id(1)));
%   id = id(~ismember(id,id(1)));
% end
  
end

