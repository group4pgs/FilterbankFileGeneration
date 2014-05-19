'''
@author: Elmarie
'''

import numpy as np
import matplotlib.pyplot as plt
import gaussianProcesses as GP


#np.random.seed(587482)
seed1 = np.random.seed()
seed2 = np.random.seed()
seed3 = np.random.seed()
seed4 = np.random.seed()
seed5 = np.random.seed()

l=5
z1=GP.squaredExponentialKernel(0,1,l,1000,10,seed1)
z2=GP.squaredExponentialKernel(3,1,l,1000,10,seed2)
z3=GP.squaredExponentialKernel(6,1,l,1000,10,seed3)
z4=GP.squaredExponentialKernel(9,1,l,1000,10,seed4)
z5=GP.squaredExponentialKernel(12,1,l,1000,10,seed5)

dt = 0.01
t = np.arange(0, (10+dt), dt)
plt.xlabel('Input')
plt.ylabel('Output')
plt.grid(True)
#plt.title(r'$\lambda = 1$')
plt.title('SE')
plt.plot(t,z1,'b',t,z2,'g',t,z3,'k',t,z4,'y',t,z5,'r')
plt.show()


l=5;
a=0.5
z1=GP.rationalQuadraticKernel(0,1,l,a,1000,10,seed1)
z2=GP.rationalQuadraticKernel(3,1,l,a,1000,10,seed2)
z3=GP.rationalQuadraticKernel(6,1,l,a,1000,10,seed3)
z4=GP.rationalQuadraticKernel(9,1,l,a,1000,10,seed4)
z5=GP.rationalQuadraticKernel(12,1,l,a,1000,10,seed5)

dt = 0.01
t = np.arange(0, (10+dt), dt)
plt.xlabel('Input')
plt.ylabel('Output')
plt.grid(True)
plt.title('RQ')
plt.plot(t,z1,'b',t,z2,'g',t,z3,'k',t,z4,'y',t,z5,'r')
plt.show()


l=2;
a=5
z1=GP.periodicSquaredExponential(0,1,l,a,1000,10,seed1)
z2=GP.periodicSquaredExponential(3,1,l,a,1000,10,seed2)
z3=GP.periodicSquaredExponential(6,1,l,a,1000,10,seed3)
z4=GP.periodicSquaredExponential(9,1,l,a,1000,10,seed4)
z5=GP.periodicSquaredExponential(12,1,l,a,1000,10,seed5)

dt = 0.01
t = np.arange(0, (10+dt), dt)
plt.xlabel('Input')
plt.ylabel('Output')
plt.grid(True)
plt.title('SE(sin)')
plt.plot(t,z1,'b',t,z2,'g',t,z3,'k',t,z4,'y',t,z5,'r')
plt.show()

