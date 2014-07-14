#!/usr/bin/python
'''
@author: Elmarie
'''
###########################################################################
# 1. Include predefined libraries.
###########################################################################
import numpy as np
import matplotlib.pyplot as plt
import struct
import argparse
import sys
import os.path
import re
from matplotlib.pylab import *
import matplotlib as mpl
from scipy.signal import  filtfilt, get_window
from matplotlib.patches import Polygon

###########################################################################
# 2. Include sub-functions.
###########################################################################


from noise_BaseLineDrift import noise_BaseLineDrift
from noise_Impulse import noise_Impulse
from noise_Narrowband import noise_Narrowband
import quantizationOfSignalValues as QZ
from DiagnosticPlot import *


###########################################################################
# 3. Help function -- how to set the variables in the arguments list.
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
    "   --tobs,-T           Total observation time, s (def=270)\n"\
    "   --tsamp,-t          Sampling interval, us (def=64)\n"\
    "   --mjd,-m            MJD of first sample (def=56000.0)\n"\
    "   --fch1,-F           Frequency of channel 1, MHz (def=1581.804688)\n"\
    "   --foff,-f           Channel bandwidth, MHz (def=-0.390625)\n"\
    "   --nchans,-c         Output number of channels (def=1024)\n"\
    "   --noiseInput,-n     File name: file containing noise specifications\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "   --plot,-p           Diagnostic plot (Yes or No) (def=No)\n"\
    "\n"\
    "Default parameters make a HTRU-style data file.\n"\
    "\n");
    return exit(1)


###########################################################################
# 4. Check if input file exists.
###########################################################################

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!"%arg)
    else:
        return DecipherInputTextFile(arg)


###########################################################################
# 5. Reads from input file and sets the necessary variables.
###########################################################################

def DecipherInputTextFile(inputFile):
    global lamda, I_Occurrences, I_tStart, I_tEnd, I_Magnitude, N_Occurrences, N_FStart, N_FEnd, N_tStart, N_tEnd, N_Magnitude

# Read from a text file the sources of noise
    filename = inputFile
    with open(filename, 'rt') as f:
        for row in f:
            if ('>B' and 'Lambda') in row:
                lamda = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('Occurrence' in row)):
                I_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('t_start' in row)):
                I_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('t_end' in row)):
                I_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('Magnitude' in row)):
                I_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('Occurrence' in row)):
                N_Occurrences = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('F_start' in row)):
                N_FStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('F_end' in row)):
                N_FEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('t_start' in row)):
                N_tStart = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('t_end' in row)):
                N_tEnd = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('Magnitude' in row)):
                N_Magnitude = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", row))
                print(N_Magnitude)
    return void

###########################################################################
# 7. Start of the main function: fake_noise.py
###########################################################################

if __name__ == '__main__':

###########################################################################
# 8. Set the default values for generating fake_noise.
###########################################################################


    telescope_id    = 4
    nchans          = 128              #Number of frequency channels across observed band
    obstime         = 270.0            #Observation time in seconds
    tsamp           = 50               #Sampling period microseconds
    fch1            = 1550.0           #Frequency of highest recorded channel in MHz
    foff            = -0.078125        #Bandwidth of each channel as negative value in MHz
    nbits           = 8                #Number of bits {1,2,8,16}
    nifs            = 1
    nbeams          = 1
    ibeam           = 1
    tstart          = 56000.25
    seed            = np.random.seed()  #Sets the seed value for the number generator by using the current time
    source_name     = "SKA sim v1"
    diagnosticplot  = "No"
    outputFile      = "output.fil"

    global lamda, I_Occurrences, I_tStart, I_tEnd, I_Magnitude, N_Occurrences, N_FStart, N_FEnd, N_tStart, N_tEnd, N_Magnitude

    lamda           =0.0
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
###########################################################################
# 9.Parse the arguments passed to fake_noise from the terminal.
###########################################################################


    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--info", help="Information regarding the function",action="store_true")
    parser.add_argument("-T","--tobs", help="Total observation time, s (def=270)", action="store", type=float )
    parser.add_argument("-t","--tsamp", help="Sampling interval, us (def=64)", action="store", type=float )
    parser.add_argument("-m","--mjd", help="MJD of first sample (def=56000.0)", action="store", type=float )
    parser.add_argument("-F","--fch1", help="Frequency of channel 1, MHz (def=1581.804688)", action="store", type=float )
    parser.add_argument("-f","--foff", help=" Channel bandwidth, MHz (def=-0.390625)", action="store", type=float )
    parser.add_argument("-c","--nchans", help="Output number of channels (def=1024)", type=int)
    parser.add_argument("-n","--noiseInput", help="File name: file containing noise specifications",
                        metavar="FILE", type=lambda x: is_valid_file(parser,x))
    parser.add_argument("-S","--seed", help="Random seed (def=time())", action="store", type=float )
    parser.add_argument("-s","--name", help="Source name for header (def=FAKE)", action="store")
    parser.add_argument("-o","--outFile", help="Output file name (def=output.fil)", action="store")
    parser.add_argument("-p","--plot", help="Diagnostic plot (Yes or No) (def=No)", action="store")

    args = parser.parse_args()

    if args.info:
        parser.print_help()
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
    if args.seed:
        seed=args.seed
    if args.name:
        source_name=args.name
    if args.outFile:
        outputFile=args.outFile
    if args.plot:
        diagnosticplot=args.plot
        
###########################################################################
# 10. Write the Header part of the output file
###########################################################################
        
        
        
###########################################################################
# 11. Generate Baseline drift noise, Impulse noise and Narrowband noise
###########################################################################
    numerOfSamples=np.floor(obstime/(tsamp*1e-06))
    seedValueForImpulseNoise=np.random.seed()

    for k in range(0, nchans):
        # 11.1 Generate Baseline drift noise
        noise_BaseLineDrift(1, lamda, numerOfSamples, obstime, seed)
        # 11.2 Generate impulse noise
        for m in range(0,I_Occurrences):
            out3=noise_Impulse(1, nrOfSamples,timeDuration, k)
            out4=QZ.quantizationOfImpulseNoise(8,out3)
            np.random.seed(500*k)
            temp2=np.random.randint(k*temp,(k+1)*temp)
            out2[temp2:(temp2+nrOfSamples),0]=out4[:,0]
        # 11.3 Generate narrowband noise
        for n in range(0,N_Occurrences)


###########################################################################
# 12. Produce a diagnostic plot if specified
###########################################################################
if 

###########################################################################
# 13. Write the values to a binary file
###########################################################################














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






