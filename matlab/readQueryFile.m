function [lats,longs,times,pm2_5] = readQueryFile(fileName, neglectMissing)

if (nargin<2)
  neglectMissing = false;
end
fid = fopen(fileName,'r');
rows = textscan(fid,'%s');
fclose(fid);
rows=rows{1};
nr = length(rows);

IDs = textscan(rows{1},'%q','delimiter', ',');
IDs = IDs{1}(3:end);
nIDs = length(IDs);

models = textscan(rows{2},'%q','delimiter', ',');
models = models{1}(3:end);

latsTxt = textscan(rows{3},'%q','delimiter', ',');
latsTxt = latsTxt{1}(3:end);

longsTxt = textscan(rows{4},'%q','delimiter', ',');
longsTxt = longsTxt{1}(3:end);

lats  = zeros(length(latsTxt),1);
longs = zeros(length(longsTxt),1);
for i=1:length(latsTxt)
  lats(i)  = str2double(latsTxt{i});
  longs(i) = str2double(longsTxt{i});
end
  
nt = nr-4;
pm2_5 = zeros(nt,nIDs);
times = zeros(nt,1);
for r=5:nr
  row = textscan(rows{r},'%q','delimiter', ',');
  date =  datetime(row{1}{1},'InputFormat','uuuu-MM-dd''T''HH:mm:ss''Z','TimeZone','America/Denver');
  if r==5
    startDate = date;
  end
  times(r-4) = (date - startDate) / hours(1);
  txts = row{1}(3:end);
  for c=1:length(txts)
    if (isempty(txts{c}))
      pm2_5(r-4,c) = -1;
    else
      pm2_5(r-4,c) = calibrate(str2double(txts{c}),models{c});
    end
  end
end

if (~neglectMissing)
  for j=1:nIDs
    if sum(pm2_5(:,j)<=0)~=0
      for i=1:nt
        if pm2_5(i,j)<=0 && i==1
          pm2_5(i,j) = pm2_5(find(pm2_5(:,j)>0,1),j);
        elseif pm2_5(i,j)<=0
          z=i+1;
          while z<=nt && pm2_5(z,j)<=0
            z=z+1;
          end
          if z<=nt
            pm2_5(i,j) = pm2_5(i-1,j) + (pm2_5(z,j)-pm2_5(i-1,j))/(z-(i-1));
          else
            for tmp=i:nt
               pm2_5(tmp,j) = pm2_5(i-1,j);
            end
            break;
          end
        end
      end
    end
  end
end

end
