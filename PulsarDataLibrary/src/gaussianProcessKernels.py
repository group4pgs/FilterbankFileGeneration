#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

import numpy as np
from scipy.linalg import toeplitz

# Different kernel functions for generating kernel matrices to be used for drawing functions from
# Gaussian processes (GP).


########################################################################################################
# 1. White noise kernel function
########################################################################################################


def whiteNoiseKernel(variance,numberOfSamples):

# This kernel allows us to add uncertainty to our observed data and is typically
# added to other kernels.

# The hyper-parameter "variance" specify the variance of the white noise.
# The other 2 arguments specify the parameters necessary for the simulation

    cov=np.eye(numberOfSamples,numberOfSamples)*variance

    return cov



########################################################################################################
# 2. Squared exponential (SE) kernel function
########################################################################################################

def squaredExponentialKernel(height, lamda, numberOfSamples, timeDurationOfSimulation):

# This kernel give rather smooth variations with a time-scale of lambda.

# The hyper-parameter "height" govern the output scale/magnitude/gain of the function
# The hyper-parameter "lambda" govern the input, or time, scale.
# The other 2 arguments specify the parameters necessary for the simulation

    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples
    row=[]


    start=0;
    for x in range(start,numberOfSamples):
        temp=start-scalingOfTimeInstances*(x)
        temp1=np.power((temp/lamda),2)
        temp2=np.power(height,2)*np.exp(-1*temp1)
        row.append(temp2)


    cov=toeplitz(row)


    return cov


########################################################################################################
# 3. Rational quadratic kernel (RQ) kernel function
########################################################################################################

def rationalQuadraticKernel(height,lamda,alpha,numberOfSamples,timeDurationOfSimulation ):

    # This kernel function produces a scale mixture of squared exponential kernels with different length scales.
    # The length scales are distributed according to a Beta distribution with parameters alpha and lambda^-2.
    # This gives variations with a range of time-scales, the distribution peaking around "lambda" but extending
    # to significantly longer period.

    # The hyper-parameter "height" govern the output scale/magnitude/gain of the function
    # The hyper-parameter "lambda" govern the input, or time, scale.
    # The hyper-parameter "alpha" is know as the index parameter.
    # The other 2 arguments specify the parameters necessary for the simulation

    #mean= np.ones(numberOfSamples)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    row=[]
    start=0;

    for x in range(start, numberOfSamples):
        temp1=np.power((start*x-scalingOfTimeInstances*x),2)
        temp2=alpha*np.power(lamda,2)
        temp3=1+(temp1/temp2)
        temp4=np.power(temp3,(-1*alpha))
        temp5=np.power(height,2)*temp4
        row.append(temp5)


    cov=toeplitz(row)


    return cov

########################################################################################################
# 4. Quadratic exponential + Constant + Linear kernel function
########################################################################################################

def quadraticExponential_plusConstant_plusLinear(height,var,constant,linear,numberOfSamples,timeDurationOfSimulation):

    # The hyper-parameter "height" govern the output amplitude of the function
    # The hyper-parameter "var" governs the rate of change of the exponential function
    # The hyper-parameter "constant" serves as an offset for the covariance kernel function
    # The hyper-parameter "linear" specifies the contribution of the linear component
    # The other 2 arguments specify the parameters necessary for the simulation

    #mean= np.ones(numberOfSamples)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples,numberOfSamples))

    for x in range(0, numberOfSamples):
        for y in range(x,numberOfSamples):
            temp=(scalingOfTimeInstances*x-scalingOfTimeInstances*y)
            temp1=(-1*var/2)*np.power(temp,2)
            temp2=height*np.exp(temp1)+constant+(linear*scalingOfTimeInstances*x*scalingOfTimeInstances*y)
            cov[x][y]=temp2
            cov[y][x]=cov[x][y]

    cov=np.matrix(cov)


    return cov


########################################################################################################
# 5. Periodic squared exponential (per-SE) kernel function
########################################################################################################

def periodicSquaredExponential(height,period,roughNess,numberOfSamples,timeDurationOfSimulation):

    # The hyper-parameter "height" govern the output amplitude of the function
    # The hyper-parameter "period" govern the period of the kernel function
    # The hyper-parameter "rougNess" is similar to the hyper-parameter lambda in the stationary covariance case.
    # The other 2 arguments specify the parameters necessary for the simulation

    #mean= np.ones(numberOfSamples)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    row=[]
    start=0;

    for x in range(start, numberOfSamples):
        temp1=np.abs(((start*x-scalingOfTimeInstances*x)/period))
        temp2=np.pi*temp1
        temp3=np.power(np.sin(temp2),2)
        temp4=np.exp((-1/(2*np.power(roughNess,2)))*temp3)
        temp5=np.power(height,2)*temp4
        row.append(temp5)



    cov=toeplitz(row)

    return cov

########################################################################################################
# 6. Quasi-periodic squared exponential (QP,SE) kernel function
########################################################################################################
def quasiPeriodicSquaredExponential(height,period,roughNess,lamda,numberOfSamples,timeDurationOfSimulation):

    # The hyper-parameter "height" govern the output amplitude of the function
    # The hyper-parameter "period" govern the period of the kernel function
    # The hyper-parameter "rougNess" is similar to the hyper-parameter lambda in the stationary covariance case
    # The hyper-parameter "lamda" specifies the rate of the evolution of the periodic signal
    # The other 2 arguments specify the parameters necessary for the simulation

    #mean= np.ones(numberOfSamples)*meanValue
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    row=[]

    start=0

    for x in range(start, numberOfSamples):
        temp=(start*x-scalingOfTimeInstances*x)
        temp1=np.abs((temp/period))
        temp2=np.pi*temp1
        temp3=np.power(np.sin(temp2),2)
        temp4=np.exp((-1/(2*np.power(roughNess,2)))*temp3-(np.power(temp,2)/np.power(lamda,2)))
        temp5=np.power(height,2)*temp4
        row.append(temp5)



    cov=toeplitz(row)

    return cov




########################################################################################################
    #####################    Produce functions drawn from Gaussian Processes   #####################
########################################################################################################


def generateFunctionsFromGaussianPorces(meanValue,covarianceMatrix,numberOfSamples):

    # The hyper-parameter "meanValue", is the mean value of the Gaussian process from which functions are drawn
    # The hyper-parameter "covarianceMatrix", is the covariance matrix of the Gaussian process from which functions are drawn
    # The hyper-parameter "seedValue" serves as the seed value for the random number generator
    # The hyper-parameter "numberOfSamples"  specifies the number of samples to be drawn from the Gaussian Process

    mean= np.ones(numberOfSamples)*meanValue

    z = np.random.multivariate_normal(mean,covarianceMatrix)


    return z