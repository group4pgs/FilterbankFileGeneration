#!/usr/bin/python
 
'''
# Elmarie van Heerden
# 17 June 2014
'''
 
# List of standard libraries that need to be imported
from matplotlib.pylab import *
import matplotlib as mpl
import numpy as np
from scipy.signal import  filtfilt, get_window
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
 
# List of customized libraries that need to be imported
from noise_BaseLineDrift import noise_BaseLineDrift
from noise_Impulse import noise_Impulse
import quantizationOfSignalValues as QZ
 
 
 
time=300
samples=3000
# out=noise_BaseLineDrift(1, 30, samples, time, 400000)
# out2=QZ.quantizationOfBaseLineSignal(out)
#  
#  
# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, 3000)*24
# s=s+96
 
# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out2,'k')
# plt.ylim((0, 255))
# plt.show()
 
 
tsamp=time/3000
timeDuration=3
nrOfSamples=np.uint32(timeDuration/tsamp)
 
# out3=noise_Impulse(1, nrOfSamples,timeDuration, np.random.seed())
# plt.plot(out3)
# plt.show()
# out4=QZ.quantizationOfImpulseNoise(10,out3)
 
#  
# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, samples)*24
# s=s+96
#  
# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out4,'k')
# plt.ylim((0, 255))
# plt.show()
 
occurrences=3
temp=np.uint32(samples/occurrences)
tsamp=time/3000
timeDuration=10
nrOfSamples=np.uint32(timeDuration/tsamp)
output=[]
# for k in range(0,occurrences):
#     out3=noise_Impulse(1, nrOfSamples,timeDuration, np.random.seed())
#     out4=QZ.quantizationOfImpulseNoise(10,out3)
#     np.random.seed(k)
#     temp2=np.random.randint(k*temp,(k+1)*temp)
#     print(temp2)
#     out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
# 
# 
# dt1 = time/3000
# t1 = np.arange(0, (time), dt1)
 
# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out2,'k')
# plt.ylim((0, 255))
# plt.show()
 
 
 
 
# plt.xlabel('Time (sec)')
# plt.ylabel('Total Power')
# plt.grid(True)
# #plt.title(r'$\lambda = 1$')
# plt.title(r'Total Power of the noise per frequency channel with baseline drift and RFI noise')
# plt.plot(t1,s,'b', label=u'Total Noise Power In Existing Software (AWGN)')
# plt.hold(True)
# plt.plot(t1,out2,'k',label=u'Total Noise Power Per Frequency Channel')
# plt.hold(True)
# #plt.plot(t2,(z_pow_noNoise-mean_z_pow_noNoise+96),'r',label=u'Total Mean Noise Power Per Frequency Channel')#,t3,z3_pow,'k',t4,z4_pow,'y')
# #plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$'])
# plt.legend(loc='upper left')
# plt.ylim(0,255)
# plt.show()
# plt.hold(False)
 
 
 
 
#######################################################################################################################
###    Funky print
#######################################################################################################################
points=3000
N=points
# Define smoothing and non-linear operators to tweak each pulse trace
# The edge window gets rid of transient effects at the edges due to filtfilt
edge_window = get_window('hamming', 17)[:8]
edge_window = np.r_[edge_window, np.ones(N - 16), edge_window[::-1]]
# Smooth traces symmetrically with simple moving average filter
smooth = lambda x, M: filtfilt(np.ones(M) / M, [1.], x) * edge_window
# The nonlinearity increases the peakiness of traces
nonlin_knee = 1.
nonlin = lambda y: nonlin_knee * (np.exp(y / nonlin_knee) - 1)
x = np.arange(N)
 
 
# The CD cover has 128 pulse traces
 
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

    for k in range(0,occurrences-1):
        out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
        out4=QZ.quantizationOfImpulseNoise(8,out3)
        np.random.seed(k)
        temp2=np.random.randint(k*temp,(k+1)*temp)
        out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
    if (5 < row < 10):
        k=2
        out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
        out4=QZ.quantizationOfImpulseNoise(8,out3)
        np.random.seed(k)
        temp2=np.random.randint(k*temp,(k+1)*temp)
        out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
    y=((out2[:,0].T)-np.median(out2[:,0]))/24#*0.02-1.8
#     y = np.abs(smooth(y, 4))# + smooth(noise1, 4)
#     y = nonlin(y)
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
c=colorbar(cmmapable,label="Number of standard deviation from expected noise level")

# fig.savefig('unknown_pleasures.pdf', facecolor='k', edgecolor='k')

plt.show()



