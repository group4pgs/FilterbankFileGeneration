#!/usr/bin/python
'''
@author: Elmarie
'''

#include predefine libraries
import numpy as np
import matplotlib.pyplot as plt
import csv
import struct



def fastfake_help():
    print(
    "\nfast_fake [options] > output.fil\n"\
    "\n"\
    "Create a gausian noise fake filterbank file. Writes to stdout.\n"\
    "Michael Keith (2014) - mkeith@pulsarastronomy.net\n"\
    "\n"\
    "OPTIONS:\n"\
    "   --help, -h          This help text\n"\
    "   --tobs,-T           Total observation time, s (def=270)\n"\
    "   --tsamp,-t          Sampling interval, us (def=64)\n"\
    "   --mjd,-m            MJD of first sample (def=56000.0)\n"\
    "   --fch1,-F           Frequency of channel 1, MHz (def=1581.804688)\n"\
    "   --foff,-f           Channel bandwidth, MHz (def=-0.390625)\n"\
    "   --nbits,-b          Output number of bits, 1,2,4,8,32 (def=2)\n"\
    "   --nchans,-c         Output number of channels (def=1024)\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "\n"\
    "Default parameters make a HTRU-style data file.\n"\
    "\n");
    return exit(1);
      
# # include other PulsarDataLibrary functions
# import noise_Impulse
# import noise_BaseLineDrift
# import noise_Narrowband
# import quantizationOfSignalValues
#
#
# Set the default values for simulating noise
telescope_id    = 4
nchans          = 128              #Number of frequency channels across observed band
obstime         = 300.0             #Observation time in seconds
tsamp           = 50*1e-6          #Sampling period microseconds
fch1            = 1550.0            #Frequency of highest recorded channel in MHz
foff            = -0.078125         #Bandwidth of each channel as negative value in MHz
nbits           = 8                 #Number of bits {1,2,8,16}
nifs            = 1
nbeams          = 1
ibeam           = 1
tstart          = 56000.25
seed            = np.random.seed()  #Sets the seed value for the number generator by using the current time
source_name     = "SKA sim v1"
#
# # Read from a text file the sources of noise
#
# age=[]
# filename = 'ska1.fil'
# with open(filename, 'rt') as f:
#     reader = csv.reader(f)
#     for row in reader:
#         age.append(np.int(row[1]))
#         print (age)
#
# f = open('text.txt', 'wb')
# f.write(bytes([13]))


fastfake_help()

#Writing to a binary #

# with open('temp', 'wb') as f:
#     f.write(struct.pack('<I', x[y]))  # Big-endian, unsigned int


# infile = open('ska1.fil', 'rb+')
# infile.seek(258)
# x = infile.read()#.decode("utf-8")
# print(x)
# dt = np.dtype(np.uint8)
#
# with open('ska1.fil', 'rb') as f:
#     f.seek(258)
#     z=np.fromfile(f,dtype=dt)
# print(z[0])

# z=np.arange(128)
# for k in range (0,23):
#     z=np.concatenate((z,z))
# z1=z[0:6000000]
# print(z1[0:130])
# #z1=np.ones(6000000)
# # z2=np.uint8(np.concatenate((z,z1)))
# with open('temp.fil', 'wb') as f:
# #     if isinstance(value, int):
# #         f.write(struct.pack('i', value) # write an int
# #     elif isinstance(value, str):
#     f.write(struct.pack('<I', 12))
#     f.write(bytes("HEADER_START", 'UTF-8')) # write a string
#     f.write(struct.pack('<I', 11))
#     f.write(bytes("source_name", 'UTF-8'))
#     f.write(struct.pack('<I', 10))
#     f.write(bytes(source_name, 'UTF-8'))
#     f.write(struct.pack('<I', 10))
#     f.write(bytes("machine_id",'UTF-8'))
#     f.write(struct.pack('<I', 10))
#     f.write(struct.pack('<I', 12))
#     f.write(bytes("telescope_id",'UTF-8'))
#     f.write(struct.pack('<I', telescope_id))
#     f.write(struct.pack('<I', 9))
#     f.write(bytes("data_type", 'UTF-8'))
#     f.write(struct.pack('<II',1,4))
#     f.write(bytes("fch1", 'UTF-8'))
#     f.write(struct.pack('d',fch1))
#     f.write(struct.pack('<I', 4))
#     f.write(bytes("foff", 'UTF-8'))
#     f.write(struct.pack('d', foff))
#     f.write(struct.pack('<I', 6))
#     f.write(bytes("nchans", 'UTF-8'))
#     f.write(struct.pack('i', nchans))
#     f.write(struct.pack('<I', 5))
#     f.write(bytes("nbits", 'UTF-8'))
#     f.write(struct.pack('i', nbits))
#     f.write(struct.pack('<I', 6))
#     f.write(bytes("nbeams", 'UTF-8'))
#     f.write(struct.pack('i', nbeams))
#     f.write(struct.pack('<I', 5))
#     f.write(bytes("ibeam", 'UTF-8'))
#     f.write(struct.pack('i', ibeam))
#     f.write(struct.pack('<I', 6))
#     f.write(bytes("tstart", 'UTF-8'))
#     f.write(struct.pack('d', tstart))
#     f.write(struct.pack('<I', 5))
#     f.write(bytes("tsamp", 'UTF-8'))
#     f.write(struct.pack('d', tsamp))
#     f.write(struct.pack('<I', 4))
#     f.write(bytes("nifs", 'UTF-8'))
#     f.write(struct.pack('i', nifs))
#     f.write(struct.pack('<I', 6))
#     f.write(bytes("signed", 'UTF-8'))
#     f.write(struct.pack('B', 1))
#     f.write(struct.pack('<I', 10))
#     f.write(bytes("HEADER_END", 'UTF-8'))
#     for k in range (0,128):
#         z2=np.uint8(z1)
#         f.write(bytes(z2))
# #     with open('ska1.fil', 'rb') as f2:
# #         f2.seek(258)
# # #         z=np.fromfile(f,dtype=dt)
# #         f.write(bytes(f2.read()))
# 
# #     else:
# #         raise TypeError('Can only write str or int')
# 
# #Reading binary values from a binary file
# # dt = np.dtype(np.uint8)
# # # Reading binary values from a file
# # with open('ska1.fil', 'rb') as f:
# # 
# #     z=np.fromfile(f,dtype=dt)
# # print((z[0:400]))
# 
# 
# infile = open('temp.fil', 'rb')
# infile.seek(12000258)
# x = infile.read(300)#.decode("utf-8")
# print(x)
# 
# 
# infile1 = open('ska1.fil', 'rb')
# x1 = infile1.read(300)#.decode("utf-8")
# 
# for i in x:
#     print(i, end=', ')
# print()
# for i in x1:
#     print(i,end=', ')
# print()
# print(x1)



mu, sigma = 0, 1 # mean and standard deviation
s = np.random.normal(mu, sigma, 3000)*24
s=s+96
z=np.ones((1,3000))*96
print(z)

plt.plot(s)
plt.plot(z.T,'r')
plt.xlabel('Number of samples')
plt.ylabel('8-bit Quantized noise values')
plt.show()

# count, bins, ignored = plt.hist(s, 30, normed=True)
# plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
#                np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
#          linewidth=2, color='r')
# plt.show()





