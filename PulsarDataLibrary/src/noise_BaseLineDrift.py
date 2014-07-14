#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter


def noise_BaseLineDrift(height, lamda, numberOfSamples, timeDurationOfSimulation, seedValue):
###########################################################################
# 1. Generate baseline drift per polarization channel by convolving
#    a low-pass filter function with samples drawn from a unit variance Guassian distribution
###########################################################################


    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples
    row=[]


    start=0;
    for x in range(start,numberOfSamples):
        temp=start-scalingOfTimeInstances*(x)
        temp1=np.power((temp/lamda),2)
        temp2=np.power(height,2)*np.exp(-1*temp1)
        row.append(temp2)

    cov1 = np.float64(row)
    cov1[0] = cov1[0]+0.000001

    mask = ( cov1[:]  >1e-9 )

    cov=cov1[mask[:]]
    cov2=np.array(cov[::-1])

    window=[]
    for k in range(0,len(cov[:])):
        window.append(cov2[k])
    for k in range(1,len(cov[:])):
        window.append(cov[k])

    np.random.seed(seedValue)
    unitVarGaussSamples=np.random.normal(0,1,(numberOfSamples+(len(cov)-1)*2))
    z1=lfilter(window,1,unitVarGaussSamples)[len(window)-1::]
    z1=z1.T

###########################################################################
# 2. Add noise, with standard deviation that is proportional to the square
#    root of the mean, to the baseline drift samples.
###########################################################################

    mask=np.abs(z1)
    average=np.mean(np.abs(z1))
    mask[mask<average]=1*average
    sigma=np.sqrt(mask)+ np.abs(np.mean(z1))

    np.random.seed()
    wn1=np.multiply(np.random.normal(0,1,numberOfSamples),sigma)
    z1_noise=z1+wn1
    del wn1

    np.random.seed()
    wn2=np.multiply(np.random.normal(0,1,numberOfSamples),sigma)
    z2_noise=z1+wn2
    del wn2

    np.random.seed()
    wn3=np.multiply(np.random.normal(0,1,numberOfSamples),sigma)
    z3_noise=z1+wn3
    del wn3

    np.random.seed()
    wn4=np.multiply(np.random.normal(0,1,numberOfSamples),sigma)
    z4_noise=z1+wn4
    del wn4

    del sigma


###########################################################################
# 3. Calculate the Total noise power per frequency channel
#    Total power = X_i^2+X_r^2+Y_i^2+Y_r^2
###########################################################################

    z1_pow=np.power(z1_noise,2)
    z2_pow=np.power(z2_noise,2)
    z3_pow=np.power(z3_noise,2)
    z4_pow=np.power(z4_noise,2)

    z_pow=z1_pow+z2_pow+z3_pow+z4_pow
    z_pow=z_pow.T

    return z_pow




# out=noise_BaseLineDrift(1, 30, 3000, 300, 400000)
# plt.plot(out)
# plt.show()






