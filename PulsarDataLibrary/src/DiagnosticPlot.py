#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported
from matplotlib.pylab import *
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# List of customized libraries that need to be imported
from noise_BaseLineDrift import noise_BaseLineDrift
from noise_Impulse import noise_Impulse
import quantizationOfSignalValues as QZ

def DiagnosticPlot():

#######################################################################################################################
###    Extract all of the values from the input file
#######################################################################################################################


    time=300
    samples=3000
    occurrences=3
    temp=np.uint32(samples/occurrences)
    tsamp=time/3000
    timeDuration=10
    nrOfSamples=np.uint32(timeDuration/tsamp)

#######################################################################################################################
###    Image plot
#######################################################################################################################
    z=[]
    v = np.linspace(0, 2.0, 15, endpoint=True)
    for row in range(1024):
        out=noise_BaseLineDrift(1, 30, samples, time, 400000)
        out2=QZ.quantizationOfBaseLineSignal(out)
#         plt.plot(out2)
#         plt.show()

        for k in range(0,occurrences):
            out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
            out4=QZ.quantizationOfImpulseNoise(3,out3)
            np.random.seed(k)
            temp2=np.random.randint(k*temp,(k+1)*temp)
            out2[temp2:(temp2+nrOfSamples)]=out4[:]
        y=((out2[:].T))#-np.median(out2[:]))/24#*0.02-1.8
        z=np.hstack((z,y))
    z=np.reshape(z, (1024,3000))
    z=z/255*6.67
    print(z.shape)
    plt.imshow(z, vmin=0, vmax=6.67, origin='upper',extent=[0,300,1470,1550])
    plt.xlabel('Time(sec)')
    plt.ylabel('Frequency channels (MHz)')
#     plt.set_xticks(np.arange(0, 10, 0.5))
#     plt.set_xticklabels(np.arange(0, 10, 0.5))
    cbar=plt.colorbar()
    cbar.set_label('$\sigma$ from expected value', rotation=270, labelpad=20, y=0.5)
    plt.show()
    return 0

DiagnosticPlot()


