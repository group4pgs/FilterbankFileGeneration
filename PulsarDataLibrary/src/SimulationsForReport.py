#!/usr/bin/python

'''
# Elmarie van Heerden
# 17 June 2014
'''

# List of standard libraries that need to be imported

import numpy as np
import time
from scipy.signal import filtfilt, get_window
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# List of customized libraries that need to be imported
import gaussianProcessKernels as GP
from noise_BaseLineDrift import noise_BaseLineDrift
from noise_Impulse import noise_Impulse
import quantizationOfSignalValues as QZ



time=300
samples=3000
out=noise_BaseLineDrift(1, 30, samples, time, 400000)
out2=QZ.quantizationOfBaseLineSignal(out)


mu, sigma = 0, 1 # mean and standard deviation
s = np.random.normal(mu, sigma, 3000)*24
s=s+96

# plt.plot(s,'b')
# plt.hold(True)
# plt.plot(out2,'k')
# plt.ylim((0, 255))
# plt.show()


tsamp=time/3000
timeDuration=3
nrOfSamples=np.uint32(timeDuration/tsamp)

out3=noise_Impulse(1, nrOfSamples,timeDuration, np.random.seed())
# plt.plot(out3)
# plt.show()
out4=QZ.quantizationOfImpulseNoise(10,out3)


mu, sigma = 0, 1 # mean and standard deviation
s = np.random.normal(mu, sigma, samples)*24
s=s+96

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


# Make a figure that is exactly the size of a CD cover (12 cm x 12 cm)
fig = plt.figure(1, figsize=(12 / 2.54, 12 / 2.54))
fig.patch.set_facecolor('black')
fig.clf()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Time (sec)')
ax.set_ylabel('Filterbank data')
# avg=np.mean(z1)
# print(avg)
#signal=out
# The CD cover has 128 pulse traces
for row in range(16):
    out=noise_BaseLineDrift(1, 30, samples, time, 400000)
    out2=QZ.quantizationOfBaseLineSignal(out)
    for k in range(0,occurrences):
        out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
        out4=QZ.quantizationOfImpulseNoise(8,out3)
        np.random.seed(k)
        temp2=np.random.randint(k*temp,(k+1)*temp)
    #print(temp2)
        out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
    # The signal is positive with a mean given by the integrated pulse profile
    #noise1 = np.multiply(np.random.normal(0,1,3001),sigma)
    #noise2 = np.multiply(np.random.normal(0,1,3001),sigma)
    #noise3 = np.multiply(np.random.normal(0,1,3001),sigma)
    #noise4 = np.multiply(np.random.normal(0,1,3001),sigma)
    # The signal and noise have different smoothing factors / bandwidths
    #pow1=np.power((signal+noise1),2)
    #pow2=np.power((signal+noise2),2)
    #pow3=np.power((signal+noise3),2)
    #pow4=np.power((signal+noise4),2)
    y=(out2[:,0].T)*0.3-20
    y = smooth(y, 8)# + smooth(noise1, 4)
    #y = nonlin(y)
    # Create overlapping traces via a manual painter's algorithm
    ax.add_patch(Polygon(np.c_[x, 5*row + 2 * y], fc='k', ec='0.85', lw=0.8,
                         closed=False, zorder=-row, alpha=1.0))
ax.autoscale_view()
ax.axis('tight')
# Remove variability of top limit of plot due to peak amplitude of top traces
y_top = ax.get_ylim()[1]
print(y_top)
ax.set_position([0.185, 0.22, 0.63, 0.56 * y_top / 130.])
ax.set_ylim(-0.01 * y_top, 1.01 * y_top)
# ax.set_xlabel('Frequency Channels')
# ax.set_ylabel('Time (sec)')
ax.axis('off')
plt.show()



























# ###########################################################################
# ###    SE kernel
# ###########################################################################
#
#
# # np.random.seed(587482)
# # seed1 = np.random.seed()
# # seed2 = np.random.seed()
# # seed3 = np.random.seed()
# # seed4 = np.random.seed()
# # seed5 = np.random.seed()
# #
# t1 = time.time()
# seed1=400000
#
#
# #######################################################################################################################
# ###    SE kernel
# #######################################################################################################################
#
# l=30
# tobs=300
# points=3000
# cov1=GP.squaredExponentialKernel(1,l,points,tobs)
# cov1 = cov1.astype(np.float64, copy=False)
# z1=GP.generateFunctionsFromGaussianPorces(0, cov1, seed1, points)
# points=3000
# cov2=GP.squaredExponentialKernel(1,l,points,tobs)
# cov2 = cov2.astype(np.float64, copy=False)
# z2=GP.generateFunctionsFromGaussianPorces(0, cov2, seed1, points)
# points=3000
# cov3=GP.squaredExponentialKernel(1,l,points,tobs)
# cov3 = cov3.astype(np.float64, copy=False)
# z3=GP.generateFunctionsFromGaussianPorces(0, cov3, seed1, points)
# points=3000
# cov4=GP.squaredExponentialKernel(1,l,points,tobs)
# cov4 = cov4.astype(np.float64, copy=False)
# z4=GP.generateFunctionsFromGaussianPorces(0, cov4, seed1, points)
#
# dt1 = tobs/3000
# dt2 = tobs/3000
# # dt3 = tobs/3000
# # dt4 = tobs/3000
# #
# t1 = np.arange(0, (tobs+dt1), dt1)
# t2 = np.arange(0, (tobs+dt2), dt2)
# # t3 = np.arange(0, (tobs+dt3), dt3)
# # t4 = np.arange(0, (tobs+dt4), dt4)
# #
# # print(len(t1), len(z1), len(t2), len(z2), len(t3), len(z3), len(t4), len(z4))
# #
# # plt.xlabel('Input')
# # plt.ylabel('Output')
# # plt.grid(True)
# # #plt.title(r'$\lambda = 1$')
# # plt.title(r'Squared Exponential Kernel ($\lambda = 30$, $t_{obs}=300sec, \# Samples=3000$)')
# # plt.plot(t1,z1,'b',t2,z2,'g')#,t3,z3,'k',t4,z4,'y')
# # plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$'])
# # plt.show()
#
#
#
#
#
# #######################################################################################################################
# ###White noise with standard deviation proportional to sqrt(mean) added to samples drawn from SE kernel Gaussian process
# #######################################################################################################################
# np.random.seed()
# sigma=np.sqrt(np.abs(z1))+24
#
# points=3000
# np.random.seed()
#
# wn1=np.multiply(np.random.normal(0,1,3001),sigma)
# z1_3=z1+wn1*0.1
# np.random.seed()
#
# wn2=np.multiply(np.random.normal(0,1,3001),sigma)
# z2_3=z2+wn2*0.1
# np.random.seed()
#
# wn3=np.multiply(np.random.normal(0,1,3001),sigma)
# z3_3=z3+wn3*0.1
# np.random.seed()
#
# wn4=np.multiply(np.random.normal(0,1,3001),sigma)
# z4_3=z4+wn4*0.1
#
# del wn1, wn2, wn3, wn4
#
#
# # dt1 = tobs/3000
# # dt2 = tobs/3000
# # dt3 = tobs/3000
# # dt4 = tobs/3000
# #
# # t1 = np.arange(0, (tobs+dt1), dt1)
# # t2 = np.arange(0, (tobs+dt2), dt2)
# # t3 = np.arange(0, (tobs+dt3), dt3)
# # t4 = np.arange(0, (tobs+dt4), dt4)
#
# # print(len(t1), len(z1), len(t2), len(z2), len(t3), len(z3), len(t4), len(z4))
# #
# # plt.xlabel('Input')
# # plt.ylabel('Output')
# # plt.grid(True)
# # #plt.title(r'$\lambda = 1$')
# # plt.title(r'White noise added to samples drawn from Squared Exponential Kernel ($\lambda = 30$, $t_{obs}=300sec, \# Samples=3000$)')
# # plt.plot(t1,z1_3,'k', label=u'Polarization channel')
# # plt.hold(True)
# # plt.plot(t2,z1,'r',label=u'Mean')#,t3,z3_3,'k',t4,z4_3,'y')
# # #plt.fill(np.concatenate([t1, t1[::-1]]),
# # #        np.concatenate([z2 - 1.9600 * sigma,
# # #                       (z2 + 1.9600 * sigma)[::-1]]),
# # #        alpha=.5, fc='b', ec='None')#, label='95% confidence interval')
# # # plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$','$95%$ confidence interval'])
# # plt.legend(loc='upper left')
# # plt.show()
# # plt.hold(False)
#
#
# #######################################################################################################################
# ###    Power of noise: White+SE kernel
# #######################################################################################################################
# z1_pow=np.power(z1_3,2)
# z2_pow=np.power(z2_3,2)
# z3_pow=np.power(z3_3,2)
# z4_pow=np.power(z4_3,2)
#
# # plt.xlabel('Time (sec)')
# # plt.ylabel('Power')
# # plt.grid(True)
# # #plt.title(r'$\lambda = 1$')
# # plt.title(r'Power of the noise per polarization channel ($\lambda = 30$, $t_{obs}=300sec, \# Samples=3000$)')
# # plt.plot(t1,z1_pow,'k',label=u'Noise Power Per Polarization')
# # plt.hold(True)
# # plt.plot(t2,z2_pow,'r',label=u'Mean Noise Power Per Polarization')#,t3,z3_pow,'k',t4,z4_pow,'y')
# # #plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$'])
# # plt.legend(loc='upper left')
# # plt.show()
# # plt.hold(False)
#
#
#
# #######################################################################################################################
# ###    Total Power of noise: White+SE kernel samples
# #######################################################################################################################
# z_pow=z1_pow+z2_pow+z3_pow+z4_pow
# z_pow_noNoise=np.power(z1,2)+np.power(z2,2)+np.power(z3,2)+np.power(z4,2)
# #z_pow_mean=np.mean(z_pow)
# #print(z_pow_mean)
# #z_pow=np.log(z_pow)
# #z_pow=z_pow-z_pow_mean
#
#
# mean_z_pow_noNoise=np.mean(z_pow_noNoise)
# mean_z_pow=np.mean(z_pow)
#
# mu, sigma = 0, 1 # mean and standard deviation
# s = np.random.normal(mu, sigma, 3001)*24
# s=s+96
# z=np.ones((1,3001))*96
#
# plt.xlabel('Time (sec)')
# plt.ylabel('Total Power')
# plt.grid(True)
# #plt.title(r'$\lambda = 1$')
# plt.title(r'Total Power of the noise per frequency channel')
# plt.plot(t1,s,'b', label=u'Total Noise Power In Existing Software (AWGN)')
# plt.hold(True)
# plt.plot(t1,(z_pow-mean_z_pow+96),'k',label=u'Total Noise Power Per Frequency Channel ($\chi^{2}_{4}$)')
# plt.hold(True)
# plt.plot(t2,(z_pow_noNoise-mean_z_pow_noNoise+96),'r',label=u'Total Mean Noise Power Per Frequency Channel')#,t3,z3_pow,'k',t4,z4_pow,'y')
# #plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$'])
# plt.legend(loc='upper left')
# plt.ylim(0,255)
# plt.show()
# plt.hold(False)
#
#
# #######################################################################################################################
# ###    Funky print
# #######################################################################################################################
# N=points+1
# # Define smoothing and non-linear operators to tweak each pulse trace
# # The edge window gets rid of transient effects at the edges due to filtfilt
# edge_window = get_window('hamming', 17)[:8]
# edge_window = np.r_[edge_window, np.ones(N - 16), edge_window[::-1]]
# # Smooth traces symmetrically with simple moving average filter
# smooth = lambda x, M: filtfilt(np.ones(M) / M, [1.], x) * edge_window
# # The nonlinearity increases the peakiness of traces
# nonlin_knee = 1.
# nonlin = lambda y: nonlin_knee * (np.exp(y / nonlin_knee) - 1)
# x = np.arange(N)
#
#
# # Make a figure that is exactly the size of a CD cover (12 cm x 12 cm)
# fig = plt.figure(1, figsize=(12 / 2.54, 12 / 2.54))
# fig.patch.set_facecolor('black')
# fig.clf()
# ax = fig.add_subplot(1, 1, 1)
# ax.set_xlabel('Time (sec)')
# ax.set_ylabel('Total Power per Frequency Channel')
# # avg=np.mean(z1)
# # print(avg)
# signal=z1
# # The CD cover has 128 pulse traces
# for row in range(128):
#     # The signal is positive with a mean given by the integrated pulse profile
#     noise1 = np.multiply(np.random.normal(0,1,3001),sigma)
#     #noise2 = np.multiply(np.random.normal(0,1,3001),sigma)
#     #noise3 = np.multiply(np.random.normal(0,1,3001),sigma)
#     #noise4 = np.multiply(np.random.normal(0,1,3001),sigma)
#     # The signal and noise have different smoothing factors / bandwidths
#     pow1=np.power((signal+noise1),2)
#     #pow2=np.power((signal+noise2),2)
#     #pow3=np.power((signal+noise3),2)
#     #pow4=np.power((signal+noise4),2)
#     y=pow1
#     y = smooth(y, 28)# + smooth(noise1, 4)
#     #y = nonlin(y)
#     # Create overlapping traces via a manual painter's algorithm
#     ax.add_patch(Polygon(np.c_[x, row + 2 * y], fc='k', ec='0.85', lw=0.8,
#                          closed=False, zorder=-row, alpha=1.0))
# ax.autoscale_view()
# ax.axis('tight')
# # Remove variability of top limit of plot due to peak amplitude of top traces
# y_top = ax.get_ylim()[1]
# print(y_top)
# ax.set_position([0.185, 0.22, 0.63, 0.56 * y_top / 130.])
# ax.set_ylim(-0.01 * y_top, 1.01 * y_top)
# # ax.set_xlabel('Frequency Channels')
# # ax.set_ylabel('Time (sec)')
# ax.axis('off')
#
#
# #fig.savefig('unknown_pleasures.pdf', facecolor='k', edgecolor='k')
#
# plt.show()

#######################################################################################################################
###    SE kernel+White noise kernel
#######################################################################################################################
# l=30
# tobs=300
# points=3000
# cov12=GP.whiteNoiseKernel(1,points)
# cov12 = cov12.astype(np.float64, copy=False)
# cov1_2=cov1+cov12
# z12=GP.generateFunctionsFromGaussianPorces(0, cov1_2, seed1, points)
# points=3000
# cov22=GP.whiteNoiseKernel(1,points)
# cov22 = cov22.astype(np.float64, copy=False)
# cov2_2=cov2+cov22
# z22=GP.generateFunctionsFromGaussianPorces(0, cov2_2, seed1, points)
# points=3000
# cov32=GP.whiteNoiseKernel(1,points)
# cov32 = cov32.astype(np.float64, copy=False)
# cov3_2=cov3+cov32
# z32=GP.generateFunctionsFromGaussianPorces(0, cov3_2, seed1, points)
# points=3000
# cov42=GP.whiteNoiseKernel(1,points)
# cov42 = cov42.astype(np.float64, copy=False)
# cov4_2=cov4+cov42
# z42=GP.generateFunctionsFromGaussianPorces(0, cov4_2, seed1, points)
#
# del cov12,cov22, cov32, cov42
# del cov1_2, cov2_2, cov3_2, cov4_2

#
# dt1 = tobs/3000
# dt2 = tobs/3000
# dt3 = tobs/3000
# dt4 = tobs/3000
#
# t1 = np.arange(0, (tobs+dt1), dt1)
# t2 = np.arange(0, (tobs+dt2), dt2)
# t3 = np.arange(0, (tobs+dt3), dt3)
# t4 = np.arange(0, (tobs+dt4), dt4)
#
# print(len(t1), len(z12), len(t2), len(z22), len(t3), len(z32), len(t4), len(z42))
#
# plt.xlabel('Input')
# plt.ylabel('Output')
# plt.grid(True)
# #plt.title(r'$\lambda = 1$')
# plt.title(r'Squared Exponential Kernel ($\lambda = 30$, $t_{obs}=300sec, \# Samples=3000$) + unit variance White Noise kernel')
# plt.plot(t1,z12,'b',t2,z22,'g',t3,z32,'k',t4,z42,'y')
# plt.legend(['$X_{real}$','$X_{imag}$','$Y_{real}$','$Y_{imag}$'])
# plt.show()


#######################################################################################################################
###     SE kernel with varying number of samples for the same observation time
#######################################################################################################################
#
#
# l=30
# tobs=300
# points=100
# cov1=GP.squaredExponentialKernel(1,l,points,tobs)
# cov1 = cov1.astype(np.float64, copy=False)
# z1=GP.generateFunctionsFromGaussianPorces(0, cov1, seed1, points)
# points=1001
# cov2=GP.squaredExponentialKernel(1,l,points,tobs)
# cov2 = cov2.astype(np.float64, copy=False)
# z2=GP.generateFunctionsFromGaussianPorces(0, cov2, seed1, points)
# points=2000
# cov3=GP.squaredExponentialKernel(1,l,points,tobs)
# cov3 = cov3.astype(np.float64, copy=False)
# z3=GP.generateFunctionsFromGaussianPorces(0, cov3, seed1, points)
# points=3000
# cov4=GP.squaredExponentialKernel(1,l,points,tobs)
# cov4 = cov4.astype(np.float64, copy=False)
# z4=GP.generateFunctionsFromGaussianPorces(0, cov4, seed1, points)
# points=4000
# cov5=GP.squaredExponentialKernel(1,l,points,tobs)
# cov5 = cov5.astype(np.float64, copy=False)
# z5=GP.generateFunctionsFromGaussianPorces(0, cov5, seed1, points)
#
#
# total=time.time()-t1
# print(total)
#
#
# dt1 = tobs/100
# dt2 = tobs/1000
# dt3 = tobs/2000
# dt4 = tobs/3000
# dt5 = tobs/4000
#
# t1 = np.arange(0, (tobs+dt1), dt1)
# t2 = np.arange(0, (tobs+dt2), dt2)
# t3 = np.arange(0, (tobs+dt3), dt3)
# t4 = np.arange(0, (tobs+dt4), dt4)
# t5 = np.arange(0, (tobs+dt5), dt5)
#
# print(len(t1), len(z1), len(t2), len(z2), len(t3), len(z3), len(t4), len(z4), len(t5), len(z5))
#
# plt.xlabel('Input')
# plt.ylabel('Output')
# plt.grid(True)
# #plt.title(r'$\lambda = 1$')
# plt.title(r'Squared Exponential Kernel ($\lambda = 30$, $t_{obs}=300sec$)')# , \# Samples=4000$)')
# plt.plot(t1,z1,'b',t2,z2,'g',t3,z3,'k',t4,z4,'y',t5,z5,'r')
# plt.legend(['$\# Samples=100$','$\# Samples=1000$','$\# Samples=2000$','$\# Samples=3000$','$\# Samples=4000$'])
# plt.show()