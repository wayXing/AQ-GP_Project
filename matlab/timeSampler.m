function data = timeSampler(ts,lat,long,pm2_5,freq,avgPeriod)
nID = length(ts);
nts = length(ts{1}(1):freq:ts{1}(end));
startTime = datetime(2017,1,6,00,00,00,'TimeZone','UTC');
endTime  = datetime(2017,1,18,23,40,00,'TimeZone','UTC');
avgLat   = zeros(nID,nts);
avgLong  = zeros(nID,nts);
avgTs    = zeros(nID,nts);
avgPm2_5 = zeros(nID,nts);
for id=1:nID
%   if id==33
%     pause;
%   end
%   nts = length(ts{id}(1):freq:ts{id}(end));
  avgLat(id,:)  = ones(1,nts)*lat{id}(1);
  avgLong(id,:) = ones(1,nts)*long{id}(1);
  avgTs(id,:)   = ((startTime:freq:endTime) - startTime) / freq;
  j=1;
  for t= startTime:freq:endTime
    avg = 0;
    avgCount = 0;
    if j<=length(ts{id})
      while (ts{id}(j)-t)<=avgPeriod/2
        if (t-ts{id}(j))<=avgPeriod/2
          avg = avg + pm2_5{id}(j);
          avgCount = avgCount + 1;
        end
        j=j+1;
        if j > length(ts{id})
          break;
        end
      end
    end
    avgPm2_5(id,((t-startTime)/freq+1))=avg/avgCount;
  end
end

data{1} = avgLat;
data{2} = avgLong;
data{3} = avgPm2_5;
data{4} = avgTs;

end