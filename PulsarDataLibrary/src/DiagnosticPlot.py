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

def DiagnosticPlot(inputFile):
    
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
###    Funky print
#######################################################################################################################
    x = np.arange(samples)


    cmap = mpl.cm.get_cmap('gist_rainbow')  ##I set colomap to 'jet'
    norm = mpl.colors.Normalize(vmin=0, vmax=6.67)

    # Make a figure that is exactly the size of a CD cover (12 cm x 12 cm)
    fig = plt.figure(1, figsize=(12 / 2.54, 12 / 2.54))
    fig.patch.set_facecolor('white')
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)


    for row in range(16):
        out=noise_BaseLineDrift(1, 30, samples, time, 400000)
        out2=QZ.quantizationOfBaseLineSignal(out)

        for k in range(0,occurrences):
            out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
            out4=QZ.quantizationOfImpulseNoise(8,out3)
            np.random.seed(k)
            temp2=np.random.randint(k*temp,(k+1)*temp)
            out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
        y=((out2[:,0].T)-np.median(out2[:,0]))/24#*0.02-1.8
        for j in range(1,2999):
            exampleColor = cmap( norm(y[j]))
            ax.add_patch(Polygon(np.c_[x[j:j+2], 6.67*row + 2 * y[j:j+2]], fc=exampleColor, ec=exampleColor, lw=0.8,closed=False, zorder=-row, alpha=1.0))


        ax.autoscale_view()
        ax.axis('tight')
        # Remove variability of top limit of plot due to peak amplitude of top traces
        y_top = ax.get_ylim()[1]
        ax.set_position([0.285, 0.22, 0.43, 0.56 * y_top / 42.])
        ax.set_ylim(-0.01 * y_top, 1.01 * y_top)
        ax.axis('off')

        cmmapable = mpl.cm.ScalarMappable(norm, cmap)
        cmmapable.set_array(range(0, 7))
        colorbar(cmmapable,label="Number of standard deviation from expected noise level")

        plt.show()

        return 0



