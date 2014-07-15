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
from scipy import signal
import time
###########################################################################
# 2. Include sub-functions.
###########################################################################


from noise_BaseLineDrift import noise_BaseLineDriftSmooth,noise_BaseLineDriftPower
from noise_Impulse import noise_ImpulseSmooth,noise_ImpulsePower
from noise_Narrowband import noise_NarrowbandSmooth,noise_NarrowbandPower
import quantizationOfSignalValues as QZ

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
    "   --tobs,-T           Total observation time, s (def=300)\n"\
    "   --tsamp,-t          Sampling interval, us (def=1000)\n"\
    "   --mjd,-m            MJD of first sample (def=56000.0)\n"\
    "   --fch1,-F           Frequency of channel 1, MHz (def=1550.0)\n"\
    "   --foff,-f           Channel bandwidth, MHz (def=-0.078125)\n"\
    "   --nchans,-c         Output number of channels (def=1024)\n"\
    "   --noiseInput,-n     File name: file containing noise specifications\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "   --plot,-p           Diagnostic plot (Yes or No) (def=Yes)\n"\
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

    return 0

###########################################################################
# 6. Start of the main function: fake_noise.py
###########################################################################

if __name__ == '__main__':
    t=time.time()
###########################################################################
# 7. Set the default values for generating fake_noise.
###########################################################################
    telescope_id    = 4
    nchans          = 1024            #Number of frequency channels across observed band
    obstime         = 300             #Observation time in seconds
    tsamp           = 1000            #Sampling period microseconds
    fch1            = 1550            #Frequency of highest recorded channel in MHz
    foff            = -0.078125       #Bandwidth of each channel as negative value in MHz
    nbits           = 8               #Number of bits {1,2,8,16}
    nifs            = 1
    nbeams          = 1
    ibeam           = 1
    tstart          = 56000.0
    seed            = np.random.seed()  #Sets the seed value for the number generator by using the current time
    source_name     = "FAKE"
    diagnosticplot  = "Yes"
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
# 8.Parse the arguments passed to fake_noise from the terminal.
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
# 9. Write the Header part of the output file
###########################################################################


    with open(outputFile, 'wb') as f:
        f.write(struct.pack('<I', 12))
        f.write(bytes("HEADER_START", 'UTF-8'))
        f.write(struct.pack('<I', 11))
        f.write(bytes("source_name", 'UTF-8'))
        f.write(struct.pack('<I', 10))
        f.write(bytes(source_name, 'UTF-8'))
        f.write(struct.pack('<I', 10))
        f.write(bytes("machine_id",'UTF-8'))
        f.write(struct.pack('<I', 10))
        f.write(struct.pack('<I', 12))
        f.write(bytes("telescope_id",'UTF-8'))
        f.write(struct.pack('<I', telescope_id))
        f.write(struct.pack('<I', 9))
        f.write(bytes("data_type", 'UTF-8'))
        f.write(struct.pack('<II',1,4))
        f.write(bytes("fch1", 'UTF-8'))
        f.write(struct.pack('d',fch1))
        f.write(struct.pack('<I', 4))
        f.write(bytes("foff", 'UTF-8'))
        f.write(struct.pack('d', foff))
        f.write(struct.pack('<I', 6))
        f.write(bytes("nchans", 'UTF-8'))
        f.write(struct.pack('i', nchans))
        f.write(struct.pack('<I', 5))
        f.write(bytes("nbits", 'UTF-8'))
        f.write(struct.pack('i', nbits))
        f.write(struct.pack('<I', 6))
        f.write(bytes("nbeams", 'UTF-8'))
        f.write(struct.pack('i', nbeams))
        f.write(struct.pack('<I', 5))
        f.write(bytes("ibeam", 'UTF-8'))
        f.write(struct.pack('i', ibeam))
        f.write(struct.pack('<I', 6))
        f.write(bytes("tstart", 'UTF-8'))
        f.write(struct.pack('d', tstart))
        f.write(struct.pack('<I', 5))
        f.write(bytes("tsamp", 'UTF-8'))
        f.write(struct.pack('d', tsamp))
        f.write(struct.pack('<I', 4))
        f.write(bytes("nifs", 'UTF-8'))
        f.write(struct.pack('i', nifs))
        f.write(struct.pack('<I', 6))
        f.write(bytes("signed", 'UTF-8'))
        f.write(struct.pack('B', 1))
        f.write(struct.pack('<I', 10))
        f.write(bytes("HEADER_END", 'UTF-8'))

        print("Finished writing Header to binary file")
###########################################################################
# 10. Generate Baseline drift noise, Impulse noise and Narrowband noise
###########################################################################

    z=[]

    numberOfSamples=np.uint32(obstime/(tsamp*1e-06))

    out3=[]
    listOfListsI=[]
    listOfListsN=[]

    out1=noise_BaseLineDriftSmooth(1, lamda, numberOfSamples, obstime, seed)
    print("Finished generating the smooth baseline drift")
    for m in range(0,np.uint8(I_Occurrences)):
        TimeDuration=(I_tEnd[m]-I_tStart[m])
        nrOfSamplesI=np.floor(TimeDuration/(tsamp*1e-06))
        seedValueForImpulseNoise=np.random.seed()
        out3.append(noise_ImpulseSmooth(1, np.uint16(nrOfSamplesI),TimeDuration, seedValueForImpulseNoise))
        listOfListsI.append(out3[:])
        out3.clear()
    print("Finished generating the Impulse noise occurrences")
    for n in range(0,np.uint8(N_Occurrences)):
        TimeDuration=(N_tEnd[n]-N_tStart[n])
        nrOfSamplesN=np.floor(TimeDuration/(tsamp*1e-06))
        seedValueForNarrowbandNoise=np.random.seed()
        out3.append(noise_NarrowbandSmooth(1, np.uint16(nrOfSamplesN),TimeDuration, seedValueForNarrowbandNoise))
        listOfListsN.append(out3[:])
        out3.clear()
    print("Finished generating the Narrowband noise occurrences")

    for k in range(0, nchans):
        print(k)
        print(time.time()-t)
        # 10.1 Generate Baseline drift noise
        out=noise_BaseLineDriftPower(out1, numberOfSamples)
        out2=QZ.quantizationOfBaseLineSignal(out)
        # 10.2 Generate impulse noise
        for m in range(0,np.uint8(I_Occurrences)):
            out7=noise_ImpulsePower(listOfListsI[m],np.array(np.shape(listOfListsI[m][0])))
            out4=QZ.quantizationOfImpulseNoise(I_Magnitude[m],out7)
            out2[np.floor((I_tStart[m])/(tsamp*1e-06)):(np.floor((I_tStart[m])/(tsamp*1e-06))+np.array(np.shape(listOfListsI[m][0])))]=out4[:,0]
        # 10.3 Generate narrowband noise
        for n in range(0,np.uint8(N_Occurrences)):
            diff1=fch1-N_FEnd
            diff2=fch1-N_FStart
            if ((np.ceil(diff1/np.abs(foff)) <= k) and( k <= np.floor(diff2/np.abs(foff)))):
                print(k)
                out8=noise_NarrowbandPower(listOfListsN[n],np.array(np.shape(listOfListsN[n][0])))
                out9=QZ.quantizationOfNarrowbandNoise(N_Magnitude[n],out8)
                out2[np.floor((N_tStart[n])/(tsamp*1e-06)):(np.floor((N_tStart[n])/(tsamp*1e-06))+np.array(np.shape(listOfListsN[n][0])))]=out9[:,0]

        # 10.4 Generate the diagnostic plot if diagnosticplot=="yes"
        if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
            y=signal.decimate((out2[:].T),3)
            z=np.hstack((z,y))
        # 10.5 Write the values out to the binary file
        f = open(outputFile, 'ab')
        z2=np.uint8(out2)
        f.write(bytes(z2))
        f.close()

    if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
        z=np.reshape(z, (nchans,len(y)))
        z=z/255*6.67
        print(z.shape)
        plt.imshow(z, vmin=0, vmax=6.67, origin='upper',extent=[0,obstime,(fch1+nchans*foff),fch1])
        plt.xlabel('Time(sec)')
        plt.ylabel('Frequency channels (MHz)')
        cbar=plt.colorbar()
        cbar.set_label('$\sigma$ from expected value', rotation=270, labelpad=20, y=0.5)
        plt.show()

    print(time.time()-t)
    plt.plot(out2)
    plt.show()












