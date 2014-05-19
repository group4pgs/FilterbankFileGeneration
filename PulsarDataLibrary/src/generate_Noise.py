'''
@author: Elmarie
'''

# include predefine libraries
import numpy as np
import matplotlib.pyplot as plt
import csv

# include other PulsarDataLibrary functions
import noise_Impulse
import noise_BaseLineDrift
import noise_Narrowband
import quantizationOfSignalValues


# Set the default values for simulating noise
nchans  = 1024              #Number of frequency channels across observed band
obstime = 20                #Observation time in seconds
tsamp   = 64                #Sampling period microseconds
fch1    = 1550              #Frequency of highest recorded channel in MHz
foff    = -0.078125         #Bandwidth of each channel as negative value in MHz
nbits   = 8                 #Number of bits {1,2,8,16}
seed    = np.random.seed()  #Sets the seed value for the number generator by using the current time

age=[]
filename = 'test.txt'
with open(filename, 'rt') as f:
    reader = csv.reader(f)
    for row in reader:
        age.append(np.int(row[1]))
        print (age)
