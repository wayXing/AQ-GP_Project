#import scipy.io
import numpy as np
#from negativeLogLikelihood import marginalNegLL
#from NegLLGradient import gradMNLL
import csv
from AQ_API import AQGPR

def readCSVFile(fileName):
    csvFile = open(fileName, "rb")
    reader = csv.reader(csvFile)
    output = []
    for row in reader:
        output.append([float(row[0])])
    csvFile.close()

    return output
    
def main():
    
    long_tr = readCSVFile('data/example_data/LONG_tr.csv')
    lat_tr = readCSVFile('data/example_data/LAT_tr.csv')
    time_tr = readCSVFile('data/example_data/TIME_tr.csv')
    pm2p5_tr = readCSVFile('data/example_data/PM2p5_tr.csv')
    long_Q = readCSVFile('data/example_data/LONG_Q.csv')
    lat_Q = readCSVFile('data/example_data/LAT_Q.csv')
    time_Q = readCSVFile('data/example_data/TIME_Q.csv')

    long_tr = np.matrix(long_tr)
    lat_tr = np.matrix(lat_tr)
    time_tr = np.matrix(time_tr)
    long_Q = np.matrix(long_Q)
    lat_Q = np.matrix(lat_Q)
    time_Q = np.matrix(time_Q)
    
    # This would be y_tr of the AQGPR function
    pm2p5_tr = np.matrix(pm2p5_tr)
    
    # This would be the x_tr of the AQGPR function
    x_tr = np.concatenate((lat_tr, long_tr, time_tr), axis=1)
    # This would be the xQuery of the AQGPR function
    x_Q = np.concatenate((lat_Q, long_Q, time_Q), axis=1)
    
    # we usually initialize the sigmaF0 for training as the standard deviation of the sensor measurements
    #sigmaF0=np.std(pm2p5_tr, ddof=1)
    # If we know the sigmaF from previous training we use the found parameter
    sigmaF0 = 8.3779
    L0 = [4.7273, 7.5732]
    # This is the noise variance and is being calculated from the sensor calibration data. This is hard coded in the AQGPR as well
    sigmaN = 5.81
    # This is the degree of the mean function used in the regression, we would like to have it equal to 1 for now
    basisFnDeg=1
    
    # Indicating wether we want to do training to find model parameters or not
    isTrain=False
    # Indicating wether we want to do the regression and find some estimates or not
    isRegression=True
    
    [yPred, yVar] = AQGPR(x_Q, x_tr, pm2p5_tr, sigmaF0, L0, sigmaN, basisFnDeg, isTrain, isRegression)

if __name__ == '__main__':
    main()
