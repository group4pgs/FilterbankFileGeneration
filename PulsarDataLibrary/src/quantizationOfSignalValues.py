#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported
import numpy as np

# def quantizationOfBaseLineSignal(inputSignal, normalizationFactor, nbits):
#
#     #maksimum = np.max(inputSignal)
#     outputSignal=inputSignal/normalizationFactor*(3.8*(np.power(2,nbits)*0.09375))
#     minValue  = np.min(outputSignal)
#     outputSignal = np.uint16(outputSignal -minValue + (np.power(2,nbits)*0.375))
#     return outputSignal

def quantizationOfBaseLineSignal(inputSignal, nbits):

    outputSignal=inputSignal/(4*np.power(156,2))*np.power(2,nbits)
#    outputSignal = np.array(outputSignal, dtype=np.int8)
    return outputSignal

def quantizationOfBaseLineSignalPlot(inputSignal, nbits):

    maksimum = np.max(inputSignal)
    outputSignal=inputSignal/maksimum*(3.8*(np.power(2,nbits)*0.09375))
    meanValue  = np.mean(outputSignal)
    outputSignal = np.uint16(outputSignal -meanValue + (np.power(2,nbits)*0.375))
    return outputSignal


def quantizationOfImpulseNoise(height,inputSignal,normalizationFactor, nbits):

    #maksimum = np.max(inputSignal)
    outputSignal=inputSignal/normalizationFactor*(height*(np.power(2,nbits)*0.0859375))
    average  = np.min(outputSignal)
    outputSignal = np.uint16(outputSignal -average + (np.power(2,nbits)*0.46875))

    return outputSignal

def quantizationOfImpulseNoisePlot(height,inputSignal,nbits):

    maksimum = np.max(inputSignal)
    outputSignal=inputSignal/maksimum*(height*(np.power(2,nbits)*0.0859375))
    average  = np.mean(outputSignal)
    outputSignal = np.uint16(outputSignal -average + (np.power(2,nbits)*0.46875))


    return outputSignal



def quantizationOfNarrowbandNoise(height,inputSignal, normalizationFactor,nbits):

    #maksimum = np.max(inputSignal)
    outputSignal=inputSignal/normalizationFactor*(height*(np.power(2,nbits)*0.0859375))
    average  = np.min(outputSignal)
    outputSignal = np.uint16(outputSignal -average + (np.power(2,nbits)*0.46875))


    return outputSignal

def quantizationOfNarrowbandNoisePlot(height,inputSignal,nbits):

    maksimum = np.max(inputSignal)
    outputSignal=inputSignal/maksimum*(height*(np.power(2,nbits)*0.0859375))
    average  = np.mean(outputSignal)
    outputSignal = np.uint16(outputSignal -average + (np.power(2,nbits)*0.46875))


    return outputSignal

