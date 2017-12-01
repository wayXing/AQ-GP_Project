import csv
import json
import os
import sys

from datetime import datetime
from datetime import timedelta

from influxdb import InfluxDBClient


TIMESTAMP = datetime.now().isoformat()


def getConfig():
    with open(sys.path[0] + '.\..\config\config.json', 'r') as configfile:
        return json.loads(configfile.read())
    sys.stderr.write('%s\tProblem reading config file.\n' % TIMESTAMP)
    sys.exit(1)


def writeLoggingDataToFile(data):
    fileName = 'oneWeekPurpleAirData.csv'
    with open(fileName, 'ab') as csvFile:
        writer = csv.writer(csvFile, delimiter=',', quoting=csv.QUOTE_ALL)
        writer.writerow(data)


def generateHourlyDates(start, end, delta):

    result = []
    start += delta
    while start < end:
        result.append(start.strftime('%Y-%m-%dT%H:%M:%SZ'))
        start += delta
    result.append(end.strftime('%Y-%m-%dT%H:%M:%SZ'))

    return result


if __name__ == "__main__":

    # using CURL to get the data:
    # curl -G 'http://air.eng.utah.edu:8086/query'
    # --data-urlencode "db=defaultdb" --data-urlencode "chunked=true"
    # --data-urlencode "chunk_size=20000"
    # --data-urlencode "q=SELECT * FROM airQuality
    # WHERE Source = 'Purple Air' AND time >= '2017-04-11T00:00:00.000000000Z'"

    # simplified bbox
    # from: https://gist.github.com/mishari/5ecfccd219925c04ac32
    print "Starting date: " + sys.argv[1]
    print "Ending date: " + sys.argv[2]
    print "Binning frequency: " + sys.argv[3]

    startDate = datetime.strptime(sys.argv[1],'%Y-%m-%d')
    endDate   = datetime.strptime(sys.argv[2],'%Y-%m-%d')

    binFreqT = datetime.strptime(sys.argv[3],'%H:%M:%S')
    binFreq  = binFreqT.hour*3600 + binFreqT.minute*60 + binFreqT.second
    
    if len(sys.argv)>=7:
        print "Geographic area [top left bottom right]: ["+sys.argv[4]+", "+sys.argv[5]+", "+sys.argv[6]+", "+sys.argv[7]+"]"
        utahBbox = {
            'left': float(sys.argv[5]),
            'right': float(sys.argv[7]),
            'bottom': float(sys.argv[6]),
            'top': float(sys.argv[4])
        }
    else:
        print "Geographic area [top left bottom right]: [-109.041069030762 36.9979667663574 -114.053932189941 42.0013885498047]"
        utahBbox = {
            'left': 36.9979667663574,
            'right': 42.0013885498047,
            'bottom': -114.053932189941,
            'top': -109.041069030762
        }

    try:
        os.remove('oneWeekPurpleAirData.csv')
    except OSError:
        pass

    #config = getConfig()
    client = InfluxDBClient(
        'air.eng.utah.edu',
        8086,
        '',
        '',
        'defaultdb',
        ssl=True,
        verify_ssl=True
    )

    hourlyDates = generateHourlyDates(startDate, endDate, timedelta(seconds=500*binFreq))
    initialDate = startDate.strftime('%Y-%m-%dT%H:%M:%SZ')
    finalDate   = endDate.strftime('%Y-%m-%dT%H:%M:%SZ')

    result = client.query('SELECT "pm2.5 (ug/m^3)","ID","Longitude","Latitude" FROM airQuality WHERE "Sensor Source" = \'Purple Air\' AND time >= \'' + initialDate + '\' AND time <= \'' + finalDate + '\';')
    result = list(result.get_points())
    
    uniqueIDs = []
    latitudes = []
    longitudes = []
    for row in result:
        if row['Latitude'] is None or row['Longitude'] is None:
            continue

        if not((float(row['Longitude']) < float(utahBbox['top'])) and (float(row['Longitude']) > float(utahBbox['bottom']))) or not((float(row['Latitude']) > float(utahBbox['left'])) and (float(row['Latitude']) < float(utahBbox['right']))):
            continue

        if row['ID'] not in uniqueIDs:
            uniqueIDs += [row['ID']]
            latitudes += [row['Latitude']]
            longitudes += [row['Longitude']]

    writeLoggingDataToFile(sum([[''],['ID'],uniqueIDs],[]))
    writeLoggingDataToFile(sum([[''],['Latitude'],latitudes],[]))
    writeLoggingDataToFile(sum([['time'],['Longitude'],longitudes],[]))

    #writeLoggingDataToFile([
    #    'time',
    #    'Latitude',
    #    'Longitude',
    #    'pm2.5 (ug/m^3)',
    #    'Temp (*C)',
    #    'Humidity (%)'
    #])

    for anEndDate in hourlyDates:
        print initialDate
        print anEndDate
        data=[]
        for anID in uniqueIDs:
            #print 'SELECT * FROM airQuality WHERE "Sensor Source" = \'Purple Air\' AND time >= ' + initialDate + ' AND time <= ' + anEndDate + ';'
            result = client.query('SELECT MEAN("pm2.5 (ug/m^3)") FROM airQuality WHERE "Sensor Source" = \'Purple Air\' AND time >= \'' + initialDate + '\' AND time < \'' + anEndDate  + '\' AND ID = \'' + anID + '\' group by time('+ str(binFreq)+ 's);')
            result = list(result.get_points())
            if anID==uniqueIDs[0]:
                for row in result:
                    data.append([row['time'],[''],row['mean']])
            else:
                for i in range(len(result)):
                    data[i] += [result[i]['mean']]

        for t in data:
            writeLoggingDataToFile(t)

        initialDate = anEndDate

    print 'DONE'
