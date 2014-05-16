'''
@author: Elmarie
'''

import numpy as np
import matplotlib.pyplot as plt


def gaussianProcess_suqaredExponentialKernel(lamda,height,numberOfSamples,timeDurationOfSimulation):

# The hyperparameter "height" govern the output scale/magnitude/gain of the function
# The hyperparameter "lambda" govern the input, or time, scale.
# The other 2 arguments specify the parameters necessary for the simulation

    mean= np.zeros(numberOfSamples+1)
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp=np.power(((scalingOfTimeInstances*x-scalingOfTimeInstances*y)/lamda),2)
            cov[x][y]=np.power(height,2)*np.exp(-1*temp)


    z = np.random.multivariate_normal(mean,cov)

    return z


def gaussianProcess_rationalQuadraticKernel(lamda,height,numberOfSamples,timeDurationOfSimulation):

    # The hyperparameter "height" govern the output scale/magnitude/gain of the function
    # The hyperparameter "lambda" govern the input, or time, scale.
    # The other 2 arguments specify the parameters necessary for the simulation

    mean= np.zeros(numberOfSamples+1)
    scalingOfTimeInstances=timeDurationOfSimulation/numberOfSamples

    cov=np.zeros((numberOfSamples+1,numberOfSamples+1))

    for x in range(0, numberOfSamples+1):
        for y in range(0,numberOfSamples+1):
            temp=np.power(((scalingOfTimeInstances*x-scalingOfTimeInstances*y)/lamda),2)
            cov[x][y]=np.power(height,2)*np.exp(-1*temp)


    z = np.random.multivariate_normal(mean,cov)


    return



l=1;
z1=gaussianProcess_suqaredExponentialKernel(l,1,1000,10)
z2=gaussianProcess_suqaredExponentialKernel(l,1,1000,10)
z3=gaussianProcess_suqaredExponentialKernel(l,1,1000,10)
z4=gaussianProcess_suqaredExponentialKernel(l,1,1000,10)
z5=gaussianProcess_suqaredExponentialKernel(l,1,1000,10)

dt = 0.01
t = np.arange(0, (10+dt), dt)
#plt.xlabel('Steps')
#Plt.ylabel('Magnitude')
#plt.grid(True)
plt.title(r'$\lambda = 1$')
plt.plot(t,z1,'b',t,z2,'g',t,z3,'k',t,z4,'y',t,z5,'r')
plt.show()


