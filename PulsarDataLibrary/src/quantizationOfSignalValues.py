#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported
import numpy as np

import matplotlib.pyplot as plt

from noise_BaseLineDrift import noise_BaseLineDrift
from noise_Impulse import noise_Impulse

def quantizationOfBaseLineSignal(inputSignal):

    outputSignal= inputSignal- np.median(inputSignal)
    maksimum = np.max(outputSignal)
    outputSignal=outputSignal/maksimum*(3.8*24)

    average  = np.mean(outputSignal)

    outputSignal = np.uint8(outputSignal -average + 96)

    return outputSignal

# out=noise_BaseLineDrift(1, 30, 3000, 300, 9000040)
# out2=quantizationOfBaseLineSignal(out)
#
#
# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, 3000)*24
# s=s+96
#
# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out2,'k')
# plt.ylim((0, 255))
# plt.show()


def quantizationOfImpulseNoise(height,inputSignal):

    outputSignal= inputSignal- np.median(inputSignal)
    maksimum = np.max(outputSignal)
    outputSignal=outputSignal/maksimum*(4.0*2.4*height)

    average  = np.mean(outputSignal)

    outputSignal = np.uint8(outputSignal -average + 115)

    return outputSignal


# tsamp=300/3000
# timeDuration=3
# nrOfSamples=np.uint32(timeDuration/tsamp)
#
# out=noise_Impulse(1, nrOfSamples,timeDuration, np.random.seed())
# plt.plot(out)
# plt.show()
# out2=quantizationOfImpulseNoise(10,out)
#
#
# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, nrOfSamples)*24
# s=s+96
#
# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out2,'k')
# plt.ylim((0, 255))
# plt.show()
