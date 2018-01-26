#import csv
import json
import sys

from datetime import datetime

from influxdb import InfluxDBClient


TIMESTAMP = datetime.now().isoformat()


def getConfig():
    with open(sys.path[0] + '/config/config.json', 'r') as configfile:
        return json.loads(configfile.read())
    sys.stderr.write('%s\tProblem reading config file.\n' % TIMESTAMP)
    sys.exit(1)


#def writeLoggingDataToFile(data):
#    fileName = 'queriedData.csv'
#    with open(fileName, 'ab') as csvFile:
#        writer = csv.writer(csvFile, delimiter=',', quoting=csv.QUOTE_ALL)
#        writer.writerow(data)
#
#


if __name__ == "__main__":

    #Default arguments
    timeFram = 3600               # 1 hour
    isSchool = False              # Query the status of all the sensors
    borderBox = {
    'bottom': 36.9979667663574,
    'top': 42.0013885498047,
    'left': -114.053932189941,
    'right': -109.041069030762 }    # Utah border coordinates
    
    # Reading input arguments
    nargv = len(sys.argv)
    if nargv>=1:
        i=1
        while (i<nargv):
            assert(sys.argv[i]=="-timeframe" or sys.argv[i]=="-box" or sys.argv[i]=="-school"),'argument number ' + str(i) + ' is not a valid argument.'
            if (sys.argv[i]=="-timeframe"):
                assert(nargv>=i+1),'-timeframe argument should be followed by the time frame in H:M:S format.'
                timeFrame = datetime.strptime(sys.argv[i+1],'%H:%M:%S')
                timeFrame = timeFrame.hour * 3600 + timeFrame.minute * 60 + timeFrame.second
                i += 2
            elif (sys.argv[i]=="-box"):
                assert(nargv>=i+4),'-box argument should be followed by 4 coordinates of the border box [top left bottom right].'
                borderBox = {
                'left': float(sys.argv[i+2]),
                'right': float(sys.argv[i+4]),
                'bottom': float(sys.argv[i+3]),
                'top': float(sys.argv[i+1])
                }
                i += 5
            elif(sys.argv[i]=="-school"):
                assert(nargv>=i+1 and (sys.argv[i+1]==0 or sys.argv[i+1]==1)),'The -school argument is a boolean argument and should be followed by either 0 or 1.'
                isSchool = sys.argv[i+1]
                i += 2
                
    print ("Time frame: last " + str(timeFrame) +" seconds")
    print ("Geographic area [top left bottom right]: ["+str(borderBox['top'])+", "+str(borderBox['left'])+", "+str(borderBox['bottom'])+", "+str(borderBox['right'])+"]")

    # Reading the config file
    config = getConfig()
    # Purple Air client
    pAirClient = InfluxDBClient(
        config['influx_host'],
        config['influx_port'],
        config['influx_username'],
        config['influx_pwd'],
        config['purpleAir_db'],
        ssl=True,
        verify_ssl=True
    )
    # airU client
    airUClient = InfluxDBClient(
        config['influx_host'],
        config['influx_port'],
        config['influx_username'],
        config['influx_pwd'],
        config['airu_db'],
        ssl=True,
        verify_ssl=True
    )    

    # Querying the Purple Air sensor IDs with their coordinates and sensor model
    result = pAirClient.query('SELECT "pm2.5 (ug/m^3)","ID","Longitude","Latitude","Sensor Model" FROM airQuality WHERE "Sensor Source" = \'Purple Air\' AND time >= now()-' + str(timeFrame) + 's;')
    result = list(result.get_points())
    
    pAirUniqueIDs = []
    pAirLatitudes = []
    pAirLongitudes = []
    pAirSensorModels = []
    for row in result:
#        if row['Latitude'] is None or row['Longitude'] is None:
#            print ("Skipped sensor with ID:" + row['ID'] + " -> Latitude/Longitude information not available!")
#            continue

        if not((float(row['Longitude']) < borderBox['right']) and (float(row['Longitude']) > borderBox['left'])) or not((float(row['Latitude']) > borderBox['bottom']) and (float(row['Latitude']) < borderBox['top'])):
            continue

        if row['ID'] not in pAirUniqueIDs:
            pAirUniqueIDs += [row['ID']]
            if row['Latitude'] is None:
                pAirLatitudes += ['missing']
            else:
                pAirLatitudes += [row['Latitude']]
            if row['Longitude'] is None:
                pAirLongitudes += ['missing']
            else:
                pAirLongitudes += [row['Longitude']]
            if row['Sensor Model'] is None:
                pAirSensorModels += ['missing']
            else:
                pAirSensorModels += [row['Sensor Model'].split('+')[0]]

    # Querying the airU sensor IDs with their coordinates and sensor model
    result = airUClient.query('SELECT "PM2.5","ID","SensorModel" FROM ' + config['airu_pm25_measurement'] + ' WHERE time >= now()-' + str(timeFrame) + 's;')
    result = list(result.get_points())
    
    # Querying the sensor IDs
    tmpIDs = []
    for row in result:
        if row['ID'] not in tmpIDs:
            tmpIDs += [row['ID']]

    # Querying the coordinates and model of each sensor in the queried geographic area
    airUUniqueIDs = []
    airULatitudes = []
    airULongitudes = []
    airUSensorModels = []
    for anID in tmpIDs:
        last = airUClient.query('SELECT LAST(Latitude),"SensorModel" FROM ' \
                 + config['airu_lat_measurement'] + ' WHERE ID=\'' + anID +'\' AND time >= now()-' + str(timeFrame) + 's;')
        last = list(last.get_points())[0]
        senModel = last['SensorModel']
        lat = last['last']

        last = airUClient.query('SELECT LAST(Longitude),"SensorModel" FROM ' \
                 + config['airu_long_measurement'] + ' WHERE ID=\'' + anID +'\' AND time >= now()-' + str(timeFrame) + 's;')
        last = list(last.get_points())[0]
        long = last['last']
        
#        if lat is None or long is None:
#            print ("Skipped sensor with ID:" + anID + " -> Latitude/Longitude information not available!")
#            continue
#        if lat==0 or long==0:
#            print ("Skipped sensor with ID:" + anID + " -> Latitude/Longitude has not been aquired!")
#            continue

        if not((float(long) < borderBox['right']) and (float(long) > borderBox['left'])) or not((float(lat) > borderBox['bottom']) and (float(lat) < borderBox['top'])):
            continue
            
        airUUniqueIDs += [anID]
        if lat is None:
            airULatitudes += ['missing']
        elif lat==0:
            airULatitudes += ['not aquired']
        else:
            airULatitudes += [lat]
        if long is None:
            airULongitudes += ['missing']
        elif long==0:
            airULongitudes += ['not aquired']
        else:
            airULongitudes += [long]
        if senModel is None:
            airUSensorModels += ['missing']
        else:
            airUSensorModels += [senModel.split('+')[0]]


    # Printing the status of the sensors in the required box
    print ('            \t            \t           \t             \t        Query Status         \t             ')
    print ('ID          \tSensor Model\tLatitude   \tLongitude    \toffline/failure/online (total)\tLatest Status')
    print ('------------\t------------\t-----------\t-------------\t------------------------------\t-------------')
    for i, anID in enumerate(airUUniqueIDs):
        result = airUClient.query('SELECT "PM2.5" FROM ' \
                 + config['airu_pm25_measurement'] + ' WHERE time >= now()-' \
                 + str(timeFrame)+'s AND ID = \'' + anID + '\';')
        result = list(result.get_points())
        nFail = 0;
        nOff  = 0;
        for t, res in enumerate(result):
            if res['PM2.5'] is None:
                nOff += 1
            elif res['PM2.5']<=0:
                nFail += 1
            elif t>0 and res['PM2.5']==result[t-1]['PM2.5']:
                isFail = True
                revt = t-2
                while (revt>=0 and revt>t-10):
                    if res['PM2.5']!=result[revt]['PM2.5']:
                        isFail = False
                        break
                    revt -= 1
                if isFail:
                    nFail += 1
  
        nTotal = len(result)
        nFine  = nTotal - nFail - nOff
        status = ('Offline' if (res['PM2.5'] is None) else ('Failed' if res['PM2.5']<=0 else 'Online'))
        print ('%-12s'%anID + '\t' + '%-12s'%airUSensorModels[i] + '\t' + '%-11s'%airULatitudes[i] + '\t' + '%-13s'%airULongitudes[i] \
              + '\t' + format(str(nOff)+'/'+str(nFail)+'/'+str(nFine)+' ('+str(nTotal)+')', '^30') + '\t' + status)

    for i, anID in enumerate(pAirUniqueIDs):
        result = pAirClient.query('SELECT "pm2.5 (ug/m^3)" FROM airQuality WHERE "Sensor Source" = \'Purple Air\' AND time >= now()-' \
                 + str(timeFrame)+'s AND ID = \'' + anID + '\';')
        result = list(result.get_points())
        nFail = 0;
        nOff  = 0;
        for res in result:
            if res['pm2.5 (ug/m^3)'] is None:
                nOff += 1
            elif res['pm2.5 (ug/m^3)']<=0:
                nFail += 1
        nTotal = len(result)
        nFine  = nTotal - nFail - nOff
        status = ('Offline' if (res['pm2.5 (ug/m^3)'] is None) else ('Failed' if res['pm2.5 (ug/m^3)']<=0 else 'Online'))
        print ('%-12s'%anID + '\t' + '%-12s'%pAirSensorModels[i] + '\t' + '%-11s'%pAirLatitudes[i] + '\t' + '%-13s'%pAirLongitudes[i] \
              + '\t' + format(str(nOff)+'/'+str(nFail)+'/'+str(nFine)+' ('+str(nTotal)+')', '^30') + '\t' + status)


