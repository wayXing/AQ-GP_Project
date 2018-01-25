function sortCandidates(fileName,YVAR,gridLongs,gridLats)

fid = fopen(fileName,'r');
dat = textscan(fid,'%q%q%q%q','Delimiter',',');
fclose(fid);

NAMs      = dat{1}(2:end);
categories = dat{2}(2:end);
longitudes  = dat{3}(2:end);
latitudes = dat{4}(2:end);

lats  = [];
longs = [];
outLats=[];
outLongs=[];

nN =1;
noutN = 1;
minLong = gridLongs(1,1);
maxLong = gridLongs(1,end);
minLat  = gridLats(1,1);
maxLat  = gridLats(end,1);
maxNamL = 0;
for i=1:length(categories)
  if strcmp(categories{i},'not yet selected')
    maxNamL = max(maxNamL,length(NAMs{i}));
    lat  = str2double(latitudes{i});
    long = str2double(longitudes{i});
    if (long>maxLong || long<minLong || lat>maxLat || lat<minLat)
      outNames{noutN} = NAMs{i};
      outLats  = [outLats;lat];
      outLongs = [outLongs;long];
      noutN = noutN +1;
    else
      lats=[lats;str2double(latitudes{i})];
      longs=[longs;str2double(longitudes{i})];
      names{nN} = NAMs{i};
      nN = nN+1;
    end
  end
end

YVARmean = mean(YVAR,3);
uncerts = interp2(gridLongs,gridLats,YVARmean,longs,lats);

[sortedUncerts,ind] = sort(uncerts,'descend');

display (['minLong = ',num2str(minLong),' maxLong = ',num2str(maxLong),...
          'minLat = ',num2str(minLat),' maxLat = ',num2str(maxLat)]);
display([sprintf(['%-',num2str(maxNamL),'s'],'Name'),'   ',...
         sprintf('%-13s','Longitude'),'   ',...
         sprintf('%-13s','Latitude'),'   ',...
         'Uncertainty']);
display([char(ones(1,maxNamL)*'-'),'   ',...
         char(ones(1,13)*'-'),'   ',...
         char(ones(1,13)*'-'),'   ',...
         '-----------']);     
for i=1:length(lats)
  display([sprintf(['%-',num2str(maxNamL),'s'],names{ind(i)}),'   ',...
         sprintf('%-13s',num2str(longs(ind(i)))),'   ',...
         sprintf('%-13s',num2str(lats(ind(i)))),'   ',...
         num2str(sqrt(sortedUncerts(i)))]);
end
display(' ');
display('Out of range sensors:');

for i=1:length(outLats)
  display([sprintf(['%-',num2str(maxNamL),'s'],outNames{i}),'   ',...
         sprintf('%-13s',num2str(outLongs(i))),'   ',...
         sprintf('%-13s',num2str(outLats(i)))]);
end

end