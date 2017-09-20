#!/usr/bin/python
'''
@author: Elmarie
'''
###########################################################################
# 1. Include predefined libraries.
###########################################################################
import numpy as np
import struct
from array import array
import argparse
import sys
import os.path
import os
import re
from scipy import signal
import time
from scipy.signal import lfilter
from scipy.interpolate import interp1d
import ctypes
###########################################################################
# 2. Include sub-functions.
###########################################################################
## This class determines what the available memory is for the function
## to run optimally
class memoryCheck():
    """Checks memory of a given system"""

    def __init__(self):

        if os.name == "posix":
            self.value = self.linuxRam()
        elif os.name == "nt":
            self.value = self.windowsRam()
        else:
            print "I only work with Win or Linux :P"

    def windowsRam(self):
        """Uses Windows API to check RAM in this OS"""
        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        class MEMORYSTATUS(ctypes.Structure):
            _fields_ = [
                ("dwLength", c_ulong),
                ("dwMemoryLoad", c_ulong),
                ("dwTotalPhys", c_ulong),
                ("dwAvailPhys", c_ulong),
                ("dwTotalPageFile", c_ulong),
                ("dwAvailPageFile", c_ulong),
                ("dwTotalVirtual", c_ulong),
                ("dwAvailVirtual", c_ulong)
            ]
        memoryStatus = MEMORYSTATUS()
        memoryStatus.dwLength = ctypes.sizeof(MEMORYSTATUS)
        kernel32.GlobalMemoryStatus(ctypes.byref(memoryStatus))

        return int(memoryStatus.dwAvailPhys/1024**2)

    def linuxRam(self):
        """Returns the RAM of a linux system"""
        totalMemory = os.popen("free -m").readlines()[1].split()[6]
        return int(totalMemory)


def noise_BaseLineDriftSmooth(height, lamda, numberOfSamples, timeDurationOfSimulation, seedValue):
###########################################################################
# 3. Generate baseline drift per polarization channel by convolving
#    a low-pass filter function with samples drawn from a unit variance Guassian distribution
###########################################################################

    scalingOfTimeInstances=np.float32(np.float32(timeDurationOfSimulation)/numberOfSamples)
    row=[]
    lamda=np.float32(lamda)
    stopCriteria = np.int64(np.sqrt(np.log(1e-03)*-1)*lamda/scalingOfTimeInstances)
    start=0;
    for x in range(start,stopCriteria):
        temp=np.float32(start-scalingOfTimeInstances*(x))
        temp1=np.float32(np.power((temp/lamda),2))
        temp2=np.power(height,2)*np.exp(-1*temp1)
        row.append(temp2)

    cov1 = np.float32(row)
    cov1[0] = cov1[0]+0.000001

    mask = ( cov1[:]  > 1e-03 )

    cov=cov1[mask[:]]
    cov2=np.array(cov[::-1])

    window=[]
    for k in range(0,len(cov[:])):
        window.append(cov2[k])
    for k in range(1,len(cov[:])):
        window.append(cov[k])
    print('Finished creating the smoothing kernel:'+str(len(window))+'')
    np.random.seed(seedValue)
    print('Number of Gaussian samples:'+str((numberOfSamples+(len(cov)-1)*2))+'')
    unitVarGaussSamples=np.random.normal(0,1,np.uint64(numberOfSamples+(len(cov)-1)*2))
    print('Finished drawing samples from Normal Distribution')
    z1=lfilter(window,1,unitVarGaussSamples)[len(window)-1::]
    z1=z1.T
    print('Finished the convolution of the kernel with the random samples')
    return z1

###########################################################################
# 4. Help function -- how to set the variables in the arguments list.
###########################################################################
def fakenoise_info():
    print(
    "\nfake_noise [options] --outFile NameOfOutputFile.fil\n"\
    "\n"\
    "Create a fake noise filterbank file. Writes to stdout.\n"\
    "Elmarie van Heerden (2014) - elmarie007@hotmail.com\n"\
    "\n"\
    "OPTIONS:\n"\
    "   --info, -i          This info text\n"\
    "   --tobs,-T           Total observation time, s (def=10)\n"\
    "   --tsamp,-t          Sampling interval, us (def=1000)\n"\
    "   --mjd,-m            MJD of first sample (def=56000.0)\n"\
    "   --fch1,-F           Frequency of channel 1, MHz (def=1550.0)\n"\
    "   --foff,-f           Channel bandwidth, MHz (def=-0.078125)\n"\
    "   --nchans,-c         Output number of channels (def=16)\n"\
    "   --nbits,-b          Number of bits (8, 16, 32) (def=8)\n"\
    "   --noiseInput,-n     File name: file containing noise specifications\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "   --header,-H         Write header to output file (Yes/No) (def=Yes)\n"\
    "   --stationary,-y     Produce stationary noise (Yes/No) (def=Yes)\n"\
    "   --bandPass,-p       Bandpass should have the shape specified in the noiseInput file\n"\
    "\n"\
    "Default parameters make a HTRU-style data file.\n"\
    "\n");
    return exit(1)


###########################################################################
# 5. Check if input file exists.
###########################################################################

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!"%arg)
    else:
        return DecipherInputTextFile(arg)



###########################################################################
# 6. Reads from input file and sets the necessary variables.
###########################################################################

def DecipherInputTextFile(inputFile):
    global lamda, amplitude, I_Occurrences, I_tStart, I_tEnd, I_Magnitude, N_Occurrences, N_FStart, N_FEnd, N_tStart, N_tEnd, N_Magnitude, filename
    global P_B_Occurrences, P_B_Period, P_B_DutyCycle, P_B_tStart, P_B_tEnd, P_B_Magnitude, P_N_Occurrences, P_N_Period, P_N_DutyCycle, P_N_FStart, P_N_FEnd
    global P_N_tStart, P_N_tEnd, P_N_Magnitude, bandPass_rampUp, bandPass_rampDown, bandPass_amplitude

# Read from a text file the sources of noise
    filename = inputFile
    with open(filename, 'rt') as f:
        for row in f:
            if (('Baseline' in row) and (('Lambda' in row) or ('lambda' in row))):
                lamda = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Baseline' in row) and ('Amplitude' in row)):
                amplitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

            elif (('Broadband' in row) and ('Occurrences' in row) and not('Periodic' in row)):
                I_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Broadband' in row) and ('t_start' in row) and not('Periodic' in row)):
                I_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Broadband' in row) and ('t_end' in row) and not('Periodic' in row)):
                I_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Broadband' in row) and ('Magnitude' in row) and not('Periodic' in row)):
                I_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

            elif (('Narrowband' in row) and ('Occurrences' in row) and not('Periodic' in row)):
                N_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Narrowband' in row) and ('F_start' in row) and not('Periodic' in row)):
                N_FStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Narrowband' in row) and ('F_end' in row) and not('Periodic' in row)):
                N_FEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Narrowband' in row) and ('t_start' in row) and not('Periodic' in row)):
                N_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Narrowband' in row) and ('t_end' in row) and not('Periodic' in row)):
                N_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Narrowband' in row) and ('Magnitude' in row) and not('Periodic' in row)):
                N_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

            elif (('Periodic' in row) and ('Broadband' in row) and ('Occurrences' in row)):
                P_B_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic Broadband Period' in row)):
            # elif (('Periodic' in row) and ('Broadband' in row) and ('Period' in row)):
                P_B_Period = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Broadband' in row) and ('Duty cycle' in row)):
                P_B_DutyCycle = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Broadband' in row) and ('t_start' in row)):
                P_B_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Broadband' in row) and ('t_end' in row)):
                P_B_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Broadband' in row) and ('Magnitude' in row)):
                P_B_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

            elif (('Periodic' in row) and ('Narrowband' in row) and ('Occurrences' in row)):
                P_N_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic Narrowband Period' in row)):
                P_N_Period = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('Duty cycle' in row)):
                P_N_DutyCycle = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('F_start' in row)):
                P_N_FStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('F_end' in row)):
                P_N_FEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('t_start' in row)):
                P_N_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('t_end' in row)):
                P_N_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('Periodic' in row) and ('Narrowband' in row) and ('Magnitude' in row)):
                P_N_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

            elif (('BandPass' in row) and ('ramp-up' in row)):
                bandPass_rampUp = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('BandPass' in row) and ('ramp-down' in row)):
                bandPass_rampDown = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('BandPass' in row) and ('Amplitude' in row)):
                bandPass_amplitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))

    places1_b_start=[]
    places1_b_end = []
    amp_b = []
    for idx in range(0,P_B_Occurrences):
        # Generate broadband periodic RFI
        period_b = P_B_Period[idx]
        deltaTinTime_b = (period_b)/100*P_B_DutyCycle[idx]
        Instances_b=int(np.floor((P_B_tEnd[idx]-P_B_tStart[idx])/(period_b)))
        I_Occurrences = I_Occurrences + Instances_b
        places1_b_start=np.concatenate((places1_b_start,np.arange(0,Instances_b)*(period_b)+P_B_tStart[idx]), axis=0)
        places1_b_end = np.concatenate((places1_b_end,(np.arange(0,Instances_b)*(period_b)+P_B_tStart[idx]+deltaTinTime_b)), axis=0)
        amp_b= np.concatenate((amp_b,np.ones(Instances_b)*P_B_Magnitude[idx]), axis =0)

    I_tStart = np.concatenate((I_tStart,places1_b_start ), axis=0)
    I_tEnd = np.concatenate((I_tEnd,places1_b_end ), axis=0)
    I_Magnitude = np.concatenate((I_Magnitude,amp_b ), axis=0)


    places1_n_start=[]
    places1_n_end = []
    placesbStart_n =[]
    placesbEnd_n = []
    amp_n = []

    for idx in range(0,P_N_Occurrences):
        # Generate narrowband periodic RFI
        period_n = P_N_Period[idx]
        deltaTinTime_n = (period_n)/100*P_N_DutyCycle[idx]
        Instances_n= int(np.floor((P_N_tEnd[idx]-P_N_tStart[idx])/(period_n)))
        N_Occurrences = N_Occurrences + Instances_n
        places1_n_start = np.concatenate((places1_n_start,np.arange(0,Instances_n)*(period_n)+P_N_tStart[idx]), axis=0)
        places1_n_end = np.concatenate((places1_n_end,(np.arange(0,Instances_n)*(period_n)+P_N_tStart[idx]+deltaTinTime_n)), axis=0)
        placesbStart_n= np.concatenate((placesbStart_n,np.ones(Instances_n)*P_N_FStart[idx]), axis =0)
        placesbEnd_n= np.concatenate((placesbEnd_n,np.ones(Instances_n)*P_N_FEnd[idx]), axis =0)
        amp_n= np.concatenate((amp_n,np.ones(Instances_n)*P_N_Magnitude[idx]), axis =0)


    N_tStart = np.concatenate((N_tStart,places1_n_start ), axis=0)
    N_tEnd = np.concatenate((N_tEnd,places1_n_end ), axis=0)
    N_FStart = np.concatenate((N_FStart,placesbStart_n ), axis=0)
    N_FEnd =np.concatenate((N_FEnd,placesbEnd_n ), axis=0)
    N_Magnitude = np.concatenate((N_Magnitude,amp_n ), axis=0)

    return 0


###########################################################################
# 7. Start of the main function: fake_noise.py
###########################################################################

if __name__ == '__main__':


###########################################################################
# 8. Set the default values for generating fake_noise.
###########################################################################
    telescope_id    = 4
    machine_id      = 10
    nchans          = 16            #Number of frequency channels across observed band
    obstime         = 10              #Observation time in seconds
    tsamp           = 1000           #Sampling period microseconds
    fch1            = 1550            #Frequency of highest recorded channel in MHz
    foff            = -0.078125       #Bandwidth of each channel as negative value in MHz
    nbits           = 8               #Number of bits {8,16}
    nifs            = 1
    nbeams          = 1
    ibeam           = 1
    tstart          = 56000.0
    seedValue       = np.uint32(np.random.randint(100000,size=1))
    source_name     = "Fake"
    outputFile      = "output.fil"
    header          = "Yes"
    stationary      = "Yes"
    bandPass_shape  = 0


    global lamda, amplitude, I_Occurrences, I_tStart, I_tEnd, I_Magnitude, N_Occurrences, N_FStart, N_FEnd, N_tStart, N_tEnd, N_Magnitude

    lamda           =0.0
    amplitude       =0.0
    I_Occurrences   =0.0
    I_tStart        =0.0
    I_tEnd          =0.0
    I_Magnitude     =0.0
    N_Occurrences   =0.0
    N_FStart        =0.0
    N_FEnd          =0.0
    N_tStart        =0.0
    N_tEnd          =0.0
    N_Magnitude     =0.0
    cut             = 0
    positionHeader  = 0

###########################################################################
# 9.Parse the arguments passed to fake_noise from the terminal.
###########################################################################


    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--info", help="Information regarding the function",action="store_true")
    parser.add_argument("-T","--tobs", help="Total observation time, s (def=10)", action="store", type=float )
    parser.add_argument("-t","--tsamp", help="Sampling interval, us (def=1000)", action="store", type=float )
    parser.add_argument("-m","--mjd", help="MJD of first sample (def=56000.0)", action="store", type=float )
    parser.add_argument("-F","--fch1", help="Frequency of channel 1, MHz (def=1550.0)", action="store", type=float )
    parser.add_argument("-f","--foff", help=" Channel bandwidth, MHz (def=-0.078125)", action="store", type=float )
    parser.add_argument("-c","--nchans", help="Output number of channels (def=16)", type=int)
    parser.add_argument("-b","--nbits", help="Number of bits (8,16,32) (def=8)", type=int)
    parser.add_argument("-n","--noiseInput", help="File name: file containing nDroise specifications",
                        metavar="FILE", type=lambda x: is_valid_file(parser,x))
    parser.add_argument("-S","--seed", help="Random seed (def=time())", action="store", type=float )
    parser.add_argument("-s","--name", help="Source name for header (def=FAKE)", action="store")
    parser.add_argument("-o","--outFile", help="Output file name (def=output.fil)", action="store")
    parser.add_argument("-H","--header", help="Write header to output file (Yes/No) (def=Yes)", action="store")
    parser.add_argument("-y","--stationary", help="Produce stationary noise (Yes/No) (def=Yes)", action="store")
    parser.add_argument("-p", "--bandPass", help="Bandpass should have a shape",action="store_true")

    args = parser.parse_args()

    if args.info:
        fakenoise_info()
        sys.exit(1)
    if args.tobs:
        obstime=args.tobs
    if args.tsamp:
        tsamp=args.tsamp
    if args.mjd:
        tstart=args.mjd
    if args.fch1:
        fch1=args.fch1
    if args.foff:
        foff=args.foff
    if args.nchans:
        nchans=args.nchans
    if args.nbits:
        nbits=args.nbits
    if args.seed:
        seedValue=args.seed
    if args.name:
        source_name=args.name
    if args.outFile:
        outputFile=args.outFile
    if args.header:
        header=args.header
    if args.stationary:
        stationary=args.stationary
    if args.bandPass:
        bandPass_shape = 1


    NumberOfLoops =1
    tsamp=tsamp*1e-06
###########################################################################
# 10. Generate a summary of the parameters of the observation to be printed to the screen
###########################################################################
    print
    print('##################################################################')
    print('          Parameter values for the simulation:')
    print('##################################################################')
    print"\n"\
    "   tobs        ", obstime          ,"s\n"\
    "   tsamp       ", tsamp/(1e-06)    ,"us\n"\
    "   mjd,        ", tstart           ,"\n"\
    "   fch1        ", fch1             ,"MHz\n"\
    "   foff        ", foff             ,"MHz\n"\
    "   nchans      ", nchans           ,"\n"\
    "   nbits       ", nbits           ,"\n"\
    "   noiseInput  ", filename         ,"\n"\
    "   seed        ", seedValue        ,"\n"\
    "   name        ", source_name      ,"\n"\
    "   header      ", header           ,"\n"\
    "   outputFile  ", outputFile           ,"\n"\
    "\n"
    print('##################################################################')
    print('          Parameter values from the input file:')
    print('##################################################################')
    print"\n"\
    "   Lambda                    ", lamda          ,"\n\n"\
    "   Amplitude                 ", amplitude      ,"\n"\
    "   Impulse Occurrences       ", I_Occurrences  ,"\n"\
    "   Impulse t_start           ", I_tStart       ,"\n"\
    "   Impulse t_end             ", I_tEnd         ,"\n\n"\
    "   Narrowband Occurrences    ", N_Occurrences  ,"\n"\
    "   Narrowband F_start        ", N_FStart       ,"MHz\n"\
    "   Narrowband F_end          ", N_FEnd         ,"MHz\n"\
    "   Narrowband t_start        ", N_tStart       ,"\n"\
    "   Narrowband t_end          ", N_tEnd         ,\
    "\n"
    print('##################################################################')

####################################################################if __name__ == '__main__':######
# 11. Check that all the parameters are kosher
###########################################################################
    if (np.any(I_tStart>=obstime) or np.any(I_tEnd>obstime)):
        clear = lambda: os.system('clear')
        clear()
        print "\n\n\n"
        print('##################################################################')
        print
        print('ERROR:One of the impulse RFI falls outside the observation time')
        print
        print('HINT:Change the start or end time of one of the impulse instances')
        print
        print('##################################################################')
        exit(1)
    if (np.any(N_tStart>=obstime) or np.any(N_tEnd>obstime)):
        clear = lambda: os.system('clear')
        clear()
        print "\n\n\n"
        print('##################################################################')
        print
        print('ERROR:One of the narrowband RFI falls outside the observation time')
        print
        print('HINT:Change the start or end time of one of the narrowband instances')
        print
        print('##################################################################')
        exit(1)

    if (np.any(N_FStart>=(fch1+(foff*nchans))) and np.any(N_FEnd<=foff)):
        clear = lambda: os.system('clear')
        clear()
        print "\n\n\n"
        print('##################################################################')
        print
        print('ERROR:One of the narrowband RFI falls outside the observation time')
        print
        print('HINT:Change the start or end Frequency of one of the narrowband instances')
        print
        print('##################################################################')
        exit(1)
# ###########################################################################
# # 12. Determine the available memory that can be used by this function
# #     such that it doesn't deplete all the resources.
# ###########################################################################

#    M = memoryCheck()
#    MemoryAvailable = M.value*1e6*0.7    # 1e6 converts the availble memory from Mbytes to bytes
#                                         # 0.7 limits the memory used to 70% of the available memory
#    MemoryNeeded    = obstime/tsamp*nchans*4*5   # 5 is the four polarizations + baseline
#                                                 # 4 is the number of bytes per value i.e. np.float32
#    NumberOfLoops   = np.int64(np.ceil(MemoryNeeded/MemoryAvailable))


###########################################################################
# 13. Write the Header part of the output file
###########################################################################

    if ((header=="Yes") or (header=="yes")):
        with open(outputFile, 'wb') as f:
            f.write(struct.pack('<I', len('HEADER_START')))
            unicode_array=array('b',b'HEADER_START')
            unicode_array.tofile(f)
            f.write(struct.pack('<I', len('source_name')))
            unicode_array=array('b',b'source_name')
            unicode_array.tofile(f)
            f.write(struct.pack('<I', len(source_name)))
            #formatted_source_name = '%80s' % source_name
            w = bytearray(source_name)
            unicode_array=array('b',w)
            unicode_array.tofile(f)
            f.write(struct.pack('<I', len('machine_id')))
            unicode_array=array('b',b'machine_id')
            unicode_array.tofile(f)
            f.write(struct.pack('<I', machine_id))
            f.write(struct.pack('<I', len('telescope_id')))
            unicode_array=array('b',b'telescope_id')
            unicode_array.tofile(f)
            f.write(struct.pack('<I', telescope_id))
            f.write(struct.pack('<I', len('data_type')))
            unicode_array=array('b',b'data_type')
            unicode_array.tofile(f)
            f.write(struct.pack('<II',1,len('fch1')))
            unicode_array=array('b',b'fch1')
            unicode_array.tofile(f)
            f.write(struct.pack('d',fch1))
            f.write(struct.pack('<I', len('foff')))
            unicode_array=array('b',b'foff')
            unicode_array.tofile(f)
            f.write(struct.pack('d', foff))
            f.write(struct.pack('<I', len('nchans')))
            unicode_array=array('b',b'nchans')
            unicode_array.tofile(f)
            f.write(struct.pack('i', nchans))
            f.write(struct.pack('<I', len('nbits')))
            unicode_array=array('b',b'nbits')
            unicode_array.tofile(f)
            f.write(struct.pack('i', nbits))
            f.write(struct.pack('<I', 6))
            unicode_array=array('b',b'nbeams')
            unicode_array.tofile(f)
            f.write(struct.pack('i', nbeams))
            f.write(struct.pack('<I', 5))
            unicode_array=array('b',b'ibeam')
            unicode_array.tofile(f)
            f.write(struct.pack('i', ibeam))
            f.write(struct.pack('<I', 6))
            unicode_array=array('b',b'tstart')
            unicode_array.tofile(f)
            f.write(struct.pack('d', tstart))
            f.write(struct.pack('<I', 5))
            unicode_array=array('b',b'tsamp')
            unicode_array.tofile(f)
            f.write(struct.pack('d', tsamp))
            f.write(struct.pack('<I', 4))
            unicode_array=array('b',b'nifs')
            unicode_array.tofile(f)
            f.write(struct.pack('i', nifs))
#            f.write(struct.pack('<I', 6))
#            unicode_array=array('b',b'signed')
#            unicode_array.tofile(f)
#            f.write(struct.pack('B', 1))
            f.write(struct.pack('<I', 10))
            unicode_array=array('b',b'HEADER_END')
            unicode_array.tofile(f)
            f.close()
            positionHeader = os.path.getsize(outputFile)

            print("Finished writing Header to binary file")

##########################################################################
# 14. Generate Baseline drift noise profile
###########################################################################
    if (nbits==32):
        nbits=8
        mbits=32
    else:
        mbits=nbits

    t=time.time()

    numberOfSamples=np.uint32(obstime/(tsamp))

    listOfListsI=[]
    listOfListsN=[]
    z=np.zeros(nchans).T

    # Generate the smooth baseline drift
    print('Start: Smooth baseline drift profile')
    if (stationary=='No'):
        seedValue=10 #np.uint32(np.random.randint(100000))
        points = obstime/(lamda*2/10)
        out1=noise_BaseLineDriftSmooth(1, lamda, points , obstime , seedValue)
        x = np.linspace(0, obstime,points)
        xi = np.linspace(0, obstime, numberOfSamples)
        linear = interp1d(x, out1)
        y = linear(xi)
        y = y -np.min(y)

        # std = np.round(((np.power(2,nbits))/25.6),0)
        # PosiveOffset= 3.5*std
        # dynamicRangeOfY=std/np.sqrt(nchans)*amplitude
        # y = (y-np.min(y))/(np.max(y)-np.min(y))*dynamicRangeOfY+PosiveOffset-dynamicRangeOfY/2
        # Normalize = np.max(y)+(np.sqrt(np.max(y))*1.0*3.0)
        std = np.round(((np.power(2,nbits))/10.0),0)
        PosiveOffset= 3*std
        dynamicRangeOfY=2*std/np.sqrt(nchans)*amplitude
        y = (y-np.min(y))/(np.max(y)-np.min(y))*dynamicRangeOfY+PosiveOffset-dynamicRangeOfY/2
        Normalize = 4*std+(np.sqrt(4*std)*1.0*3.0)
    else:
        seedValue=np.uint32(np.random.randint(100000))
        std = np.round(((np.power(2,nbits))/10.0),0)
        PosiveOffset= 3*std
        dynamicRangeOfY=2*std/np.sqrt(nchans)*amplitude
        y=np.ones(numberOfSamples)*(PosiveOffset) # + 0.5*dynamicRangeOfY)
        Normalize = 4*std+(np.sqrt(4*std)*1.0*3.0)
        print('End: Baseline drift profile')
        # seedValue=np.uint32(np.random.randint(100000))
        # std = np.round(((np.power(2,nbits))/25.6),0)
        # dynamicRangeOfY=std/np.sqrt(nchans)*amplitude
        # y=np.ones(numberOfSamples)*(PosiveOffset) # + 0.5*dynamicRangeOfY)
        # Normalize = np.max(y)+(np.sqrt(np.max(y))*1.0*3.0)
        # print('End: Baseline drift profile')
###########################################################################
# 15. Generate Impulse noise profile(s)
###########################################################################
    if (np.uint32(I_Occurrences)!=0):
        print('Start: Impulse noise profile(s)')
#        x1=np.arange(-3,3.000072,0.001)
#        y1=np.exp(-np.power(x1,2)/2)/np.sqrt(2*np.pi)/0.4
#        del x1
    for m in range(0,np.uint32(I_Occurrences)):

        TimeDuration=np.float32(I_tEnd[m]-I_tStart[m])
        nrOfSamplesI=np.uint32(np.floor(TimeDuration/(tsamp)))
        out3=np.sqrt(std)*I_Magnitude[m]*np.ones(nrOfSamplesI) #np.sqrt(PosiveOffset)*signal.resample(y1,nrOfSamplesI)
        location=np.uint32(np.round(I_tStart[m]/tsamp))
        y[location:(location+nrOfSamplesI)]=y[location:(location+nrOfSamplesI)]+out3[:]
        del out3
    if (np.uint32(I_Occurrences)!=0):
        print('End: Impulse noise profile(s)')



###########################################################################
# 16. Generate Narrowband noise profile(s)
###########################################################################
    if (np.uint32(N_Occurrences)!=0):
        print('Start: Narrowband noise profile(s)')
#        x1=np.arange(-3,3.000072,0.001)
#        y1=np.exp(-np.power(x1,2)/2)/np.sqrt(2*np.pi)/0.4
#        del x1
        y2=np.multiply(y,np.ones(nchans).reshape((nchans,1)))
        del y
        y=y2
    for n in range(0,np.uint32(N_Occurrences)):
        firstchannel=int((fch1-N_FEnd[n])/np.abs(foff))
        channelsaffected=np.abs(int((N_FEnd[n]-N_FStart[n])/np.abs(foff)))
        y1 = np.ones(channelsaffected)
        # x1=np.arange(1,(np.ceil(channelsaffected/2.0+1)))
        # x2=np.power(x1,(-1/np.sqrt(channelsaffected)))
        # if np.mod(channelsaffected,2)==0:
        #     y1=np.concatenate((x2[:],x2[::-1]))
        # else:
        #     y1=np.concatenate((x2[1::],x2[::-1]))
        # del x1,x2
        TimeDuration=np.float32(N_tEnd[n]-N_tStart[n])
        nrOfSamplesI=np.uint32(np.floor(TimeDuration/(tsamp)))
        out3=np.sqrt(std)*N_Magnitude[n]*np.ones(nrOfSamplesI) #*np.sqrt(PosiveOffset)*signal.resample(y1,nrOfSamplesI)
        location=np.uint32(np.round(N_tStart[n]/tsamp))

        for indx in range(0,channelsaffected):
            #y[firstchannel:(firstchannel+channelsaffected),location:(location+nrOfSamplesI)]=y[firstchannel:(firstchannel+channelsaffected),location:(location+nrOfSamplesI)]+out3[:]
            y[(firstchannel+indx),location:(location+nrOfSamplesI)]=y[(firstchannel+indx),location:(location+nrOfSamplesI)]+out3[:]*y1[indx]
        del out3, y1
    if (np.uint32(N_Occurrences)!=0):
        print("End: Narrowband noise profile(s)")

    c=[]


###########################################################################
# 17. Write the Baseline Drift noise to the output file.
###########################################################################
    # Define the function specifying the standard deviation of the baseline drift  noise

    sigma=y
    sigma = np.sqrt(sigma)
    bandPass = np.ones((nchans))

    if (bandPass_shape):
        print('The bandPass has the shape specified in the noiseInput file')
        numberOfChannelsInRampUp = int(np.floor(nchans*bandPass_rampUp))
        for teller in range(0,numberOfChannelsInRampUp):
            bandPass[teller] = (1-bandPass_amplitude)+teller*bandPass_amplitude/numberOfChannelsInRampUp;

        numberOfChannelsInRampDown = int(np.floor(nchans*bandPass_rampDown))
        for teller in range(-numberOfChannelsInRampDown,0):
            bandPass[nchans+teller] = (1-bandPass_amplitude)+np.abs(teller)*bandPass_amplitude/numberOfChannelsInRampDown;

    if ((header=="Yes") or (header=="yes")):
        f = open(outputFile, 'ab')
    else:
        f = open(outputFile, 'wb')


    NumberOfSamplesPerLoop = np.int64(np.floor(numberOfSamples/NumberOfLoops))
    NumberOfSamplesRemaining = np.int64(numberOfSamples - NumberOfSamplesPerLoop*NumberOfLoops)
    if (NumberOfSamplesRemaining !=0):
        NumberOfLoops = NumberOfLoops + 1

    for loop in range(0,NumberOfLoops):
        if (NumberOfSamplesRemaining==0.0):
            seedValue=np.uint32(np.random.randint(100000))
            np.random.seed(seedValue)
            noise= np.random.normal(0,1,(nchans*NumberOfSamplesPerLoop*4))
            if (np.uint32(N_Occurrences)!=0):
                for k in range(0, NumberOfSamplesPerLoop):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans):(k*4*nchans+nchans)]) + (y[:,idx])
                    z2_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)]) + (y[:,idx])
                    z3_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)]) + (y[:,idx])
                    z4_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)]) + (y[:,idx])
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)

                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2
            else:
                for k in range(0, NumberOfSamplesPerLoop):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= (sigma[idx]*noise[(k*4*nchans):(k*4*nchans+nchans)] + y[idx])
                    z2_noise= (sigma[idx]*noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)] + y[idx])
                    z3_noise= (sigma[idx]*noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)] + y[idx])
                    z4_noise= (sigma[idx]*noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)] + y[idx])
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)

                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2
        elif ((loop<NumberOfLoops-1) and (NumberOfSamplesRemaining>0)):
            seedValue=np.uint32(np.random.randint(100000,size=1))
            np.random.seed(seedValue)
            noise= np.random.normal(0,1,(nchans*NumberOfSamplesPerLoop*4))
            if (np.uint32(N_Occurrences)!=0):
                for k in range(0, NumberOfSamplesPerLoop):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans):(k*4*nchans+nchans)]) + (y[:,idx])
                    z2_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)]) + (y[:,idx])
                    z3_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)]) + (y[:,idx])
                    z4_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)]) + (y[:,idx])
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)
                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2
            else:
                for k in range(0, NumberOfSamplesPerLoop):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= (sigma[idx]*noise[(k*4*nchans):(k*4*nchans+nchans)] + y[idx])
                    z2_noise= (sigma[idx]*noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)] + y[idx])
                    z3_noise= (sigma[idx]*noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)] + y[idx])
                    z4_noise= (sigma[idx]*noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)] + y[idx])
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)

                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;                        
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2
        else:
            seedValue=np.uint32(np.random.randint(100000,size=1))
            np.random.seed(seedValue)
            noise= np.random.normal(0,1,(nchans*NumberOfSamplesRemaining*4))
            if (np.uint32(N_Occurrences)!=0):
                for k in range(0, NumberOfSamplesRemaining):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans):(k*4*nchans+nchans)]) + (y[:,idx]).reshape((1,nchans))
                    z2_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)]) + (y[:,idx]).reshape((1,nchans))
                    z3_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)]) + (y[:,idx]).reshape((1,nchans))
                    z4_noise= np.multiply(sigma[:,idx],noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)]) + (y[:,idx]).reshape((1,nchans))
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)

                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;                        
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2
            else:
                for k in range(0, NumberOfSamplesRemaining):
                    idx = NumberOfSamplesPerLoop*loop+k
                    z1_noise= (sigma[idx]*noise[(k*4*nchans):(k*4*nchans+nchans)] + y[idx])
                    z2_noise= (sigma[idx]*noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)] + y[idx])
                    z3_noise= (sigma[idx]*noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)] + y[idx])
                    z4_noise= (sigma[idx]*noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)] + y[idx])
                    z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
                    z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))
                    z_pow=np.multiply(z_pow,bandPass)
                    # 10.2 Write the values out to the binary file
                    if (mbits==8):
                        z_pow= (z_pow + abs(z_pow)) / 2;
                        z_pow = (z_pow + 255 - abs(z_pow - 255)) / 2;                        
                        z2=np.uint8(z_pow)
                        z2.tofile(f)
                        del z2
                    elif (mbits==16):
                        z2=np.uint16(z_pow)
                        z2.tofile(f)
                        del z2
                    else:
                        z2=array('f',list(z_pow))
                        z2.tofile(f)
                        del z2




    f.close()

    print("Finished generating and writing BaselineDrift to output file.")
    print(time.time()-t)



