'''
@author: Elmarie
'''

import numpy as np

# Different kernel functions for generating time-series drawn form Gaussian Processes are defined below.


def whiteNoiseKernel(variance,numberOfSamples):

# This kernel allows us to add uncertainty to our observed data and is typically
# added to other kernels.

# The hyper-parameter "variance" specify the variance of the white noise.
# The other 2 arguments specify the parameters necessary for the simulation

    z=np.ones((numberOfSamples+1,numberOfSamples+1))*variance

    return z


def squaredExponentialKernel(meanValue, height, lamda, numberOfSamples, timeDurationOfSimulation, seedValue):

# This kernel give rather smooth variations with a time-scale of lambda.

# The hyper-parameter "height" govern the output scale/magnitude/gain of the function
# The hyper-parameter "lambda" govern the input, or time, scale.
# The other 2 arguments specify the parameters necessary for the simulation

    mean= np.ones(numberOfSamples+1)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp=np.power(((scalingOfTimeInstances*x-scalingOfTimeInstances*y)/lamda),2)
            cov[x][y]=np.power(height,2)*np.exp(-1*temp)

    np.random.seed(seedValue)
    z = np.random.multivariate_normal(mean,cov)

    return z


def rationalQuadraticKernel(meanValue,height,lamda,alpha,numberOfSamples,timeDurationOfSimulation,seedValue ):

    # This kernel function produces a scale mixture of squared exponential kernels with different length scales.
    # The length scales are distributed according to a Beta distribution with parameters alpha and lambda^-2.
    # This gives variations with a range of time-scales, the distribution peaking around "lambda" but extending
    # to significantly longer period.

    # The hyper-parameter "height" govern the output scale/magnitude/gain of the function
    # The hyper-parameter "lambda" govern the input, or time, scale.
    # The hyper-parameter "alpha" is know as the index parameter.
    # The other 2 arguments specify the parameters necessary for the simulation

    mean= np.ones(numberOfSamples+1)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp1=np.power((scalingOfTimeInstances*x-scalingOfTimeInstances*y),2)
            temp2=alpha*np.power(lamda,2)
            temp3=1+(temp1/temp2)
            temp4=np.power(temp3,(-1*alpha))
            cov[x][y]=np.power(height,2)*temp4

    np.random.seed(seedValue)
    z = np.random.multivariate_normal(mean,cov)


    return z


def periodicSquaredExponential(meanValue,height,period,roughNess,numberOfSamples,timeDurationOfSimulation, seedValue):

    # The hyper-parameter "height" govern the output amplitude of the function
    # The hyper-parameter "period" govern the period of the kernel function
    # The hyper-parameter "rougNess" is similar to the hyper-parameter lambda in the stationary covariance case.
    # The other 2 arguments specify the parameters necessary for the simulation

    mean= np.ones(numberOfSamples+1)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp1=np.abs(((scalingOfTimeInstances*x-scalingOfTimeInstances*y)/period))
            temp2=np.pi*temp1
            temp3=np.power(np.sin(temp2),2)
            temp4=np.exp((-1/(2*np.power(roughNess,2)))*temp3)
            cov[x][y]=np.power(height,2)*temp4

    np.random.seed(seedValue)
    z = np.random.multivariate_normal(mean,cov)

    return z


def quasiPeriodicSquaredExponential(meanValue,height,period,roughNess,lamda,numberOfSamples,timeDurationOfSimulation, seedValue):

    # The hyper-parameter "height" govern the output amplitude of the function
    # The hyper-parameter "period" govern the period of the kernel function
    # The hyper-parameter "rougNess" is similar to the hyper-parameter lambda in the stationary covariance case
    # The hyper-parameter "lamda" specifies the rate of the evolution of the periodic signal
    # The other 2 arguments specify the parameters necessary for the simulation

    mean= np.ones(numberOfSamples+1)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp=(scalingOfTimeInstances*x-scalingOfTimeInstances*y)
            temp1=np.abs((temp/period))
            temp2=np.pi*temp1
            temp3=np.power(np.sin(temp2),2)
            temp4=np.exp((-1/(2*np.power(roughNess,2)))*temp3-(np.power(temp,2)/np.power(lamda,2)))
            cov[x][y]=np.power(height,2)*temp4

    np.random.seed(seedValue)
    z = np.random.multivariate_normal(mean,cov)

    return z
