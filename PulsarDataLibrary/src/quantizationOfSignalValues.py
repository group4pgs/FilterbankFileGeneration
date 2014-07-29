#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported
import numpy as np

def quantizationOfBaseLineSignal(inputSignal, normalizationFactor):

    #maksimum = np.max(inputSignal)
    outputSignal=inputSignal/normalizationFactor*(3.8*24)
    average  = np.mean(outputSignal)
    outputSignal = np.uint8(outputSignal -average + 96)
    return outputSignal


def quantizationOfImpulseNoise(height,inputSignal,normalizationFactor):

    #maksimum = np.max(inputSignal)
    outputSignal=inputSignal/normalizationFactor*(24*height)
    average  = np.mean(outputSignal)
    outputSignal = np.uint8(outputSignal -average + 125)

    return outputSignal

def quantizationOfNarrowbandNoise(height,inputSignal, normalizationFactor):

    #maksimum = np.max(inputSignal)
    outputSignal=inputSignal/normalizationFactor*(24*height)
    average  = np.mean(outputSignal)
    outputSignal = np.uint8(outputSignal -average + 125)

    return outputSignal

