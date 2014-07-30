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
from array import array
import argparse
import sys
import os.path
import re
from scipy import signal
import time
from scipy.interpolate import interp1d
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
    "   --tobs,-T           Total observation time, s (def=30)\n"\
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

###########################################################################
# 7. Set the default values for generating fake_noise.
###########################################################################
    telescope_id    = 4
    nchans          = 1024            #Number of frequency channels across observed band
    obstime         = 10             #Observation time in seconds
    tsamp           = 10000             #Sampling period microseconds
    fch1            = 1550            #Frequency of highest recorded channel in MHz
    foff            = -0.078125       #Bandwidth of each channel as negative value in MHz
    nbits           = 8               #Number of bits {1,2,8,16}
    nifs            = 1
    nbeams          = 1
    ibeam           = 1
    tstart          = 56000.0
    seed            = np.random.seed()  #Sets the seed value for the number generator by using the current time
    source_name     = "FAKE"
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
    tsamp=tsamp*1e-06

    with open(outputFile, 'wb') as f:
        f.write(struct.pack('<I', 12))
        unicode_array=array('b',b'HEADER_START')
        unicode_array.tofile(f)
        f.write(struct.pack('<I', 11))
        unicode_array=array('b',b'source_name')
        unicode_array.tofile(f)
        f.write(struct.pack('<I', 10))
        w = bytearray(source_name)
        unicode_array=array('b',w)
        unicode_array.tofile(f)
        f.write(struct.pack('<I', 10))
        unicode_array=array('b',b'machine_id')
        unicode_array.tofile(f)
        f.write(struct.pack('<I', 10))
        f.write(struct.pack('<I', 12))
        unicode_array=array('b',b'telescope_id')
        unicode_array.tofile(f)
        f.write(struct.pack('<I', telescope_id))
        f.write(struct.pack('<I', 9))
        unicode_array=array('b',b'data_type')
        unicode_array.tofile(f)
        f.write(struct.pack('<II',1,4))
        unicode_array=array('b',b'fch1')
        unicode_array.tofile(f)
        f.write(struct.pack('d',fch1))
        f.write(struct.pack('<I', 4))
        unicode_array=array('b',b'foff')
        unicode_array.tofile(f)
        f.write(struct.pack('d', foff))
        f.write(struct.pack('<I', 6))
        unicode_array=array('b',b'nchans')
        unicode_array.tofile(f)
        f.write(struct.pack('i', nchans))
        f.write(struct.pack('<I', 5))
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
        f.write(struct.pack('<I', 6))
        unicode_array=array('b',b'signed')
        unicode_array.tofile(f)
        f.write(struct.pack('B', 1))
        f.write(struct.pack('<I', 10))
        unicode_array=array('b',b'HEADER_END')
        unicode_array.tofile(f)
        f.close()

        print("Finished writing Header to binary file")

###########################################################################
# 10. Generate Baseline drift noise, Impulse noise and Narrowband noise
###########################################################################
    t=time.time()
    cut=0


    numberOfSamples=np.uint32(obstime/(tsamp))

    listOfListsI=[]
    listOfListsN=[]
    z=np.zeros(nchans).T
    averagePerSampleChannel=[]

# Generate the smooth baseline drift
    out1=noise_BaseLineDriftSmooth(1, lamda, 10000, 100, seed)

    x = np.linspace(0, 10000, 10000)
    xi = np.linspace(0, 10000, numberOfSamples)

    if (numberOfSamples>=10000):
        # use linear interpolation method
        linear = interp1d(x, out1)
        y = linear(xi)
        y = y -np.min(y)

    else:
        # use the decimate method
        y=signal.resample(out1, numberOfSamples)
        y = y -np.min(y)


    for m in range(0,np.uint8(I_Occurrences)):
        out3=[]
        TimeDuration=np.float64(I_tEnd[m]-I_tStart[m])
        nrOfSamplesI=np.floor(TimeDuration/(tsamp))

        seedValueForImpulseNoise=np.random.seed()
        out3.append(noise_ImpulseSmooth(1, np.uint16(nrOfSamplesI),TimeDuration, seedValueForImpulseNoise))
        listOfListsI.append(out3[:])
        del out3
    print("Finished generating the Impulse noise occurrences")
    for n in range(0,np.uint8(N_Occurrences)):
        out3=[]
        TimeDuration=(N_tEnd[n]-N_tStart[n])
        nrOfSamplesN=np.floor(TimeDuration/(tsamp))
        seedValueForNarrowbandNoise=np.random.seed()
        out3.append(noise_NarrowbandSmooth(1, np.uint16(nrOfSamplesN),TimeDuration, seedValueForNarrowbandNoise))

        listOfListsN.append(out3[:])
        del out3
    print("Finished generating the Narrowband noise occurrences")

    # Define the function specifying the standard deviation of the noise
    mask=np.copy(np.abs(y))
    sigma=(np.sqrt(mask)+ np.abs(np.mean(y)))*1.0
    del mask


    NormalizationValueBaseline=np.power((np.mean(np.abs(y))+(np.mean(np.abs(sigma))*2.2)),2)*4
    f = open(outputFile, 'ab')
    for k in range(0, numberOfSamples):

        # 10.1 Generate Baseline drift noise
        out=noise_BaseLineDriftPower(y[k], sigma[k], nchans)
        out2=QZ.quantizationOfBaseLineSignal(out,NormalizationValueBaseline)
        averagePerSampleChannel.append(np.mean(out2))

        del out

        # 10.2 Write the values out to the binary file
        z2=np.uint8(out2)
        z2.tofile(f)

        if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
            z=np.column_stack((z,z2))

        del z2
    f.close()


    print("Finished generating and writing BaselineDrift to output file.")

    outImage=[]
    outImage=np.copy(out2)
    toets=np.copy(out2)


    f = open(outputFile, 'rb+')
    for m in range(0,np.uint8(I_Occurrences)):

            z1=np.copy(listOfListsI[m])
            mask=np.copy(np.abs(z1))


            maksimum=np.max(mask)
            mask = -1*(mask-maksimum)
            average=np.abs(np.mean(z1))
            mask[mask<average]=1*average
            sigma=np.sqrt(mask) + np.abs(np.mean(z1))
            z1=mask
            NormalizationValue=np.power((np.mean(np.abs(z1))+(np.mean(np.abs(sigma))*2.4)),2)*4
            for k in range(0, np.shape(mask)[1]):

                out7=noise_ImpulsePower(z1[0,k],sigma[0,k],nchans)
                out4=QZ.quantizationOfImpulseNoise(I_Magnitude[m],out7,NormalizationValue)

                out4[out4[:]<averagePerSampleChannel[np.uint64(np.floor((I_tStart[m])/(tsamp)) + k)]]=averagePerSampleChannel[np.uint64(np.floor((I_tStart[m])/(tsamp)) + k)]
                del out7

                position=np.uint64(258+ np.floor((I_tStart[m])/(tsamp)) + k*nchans )
                f.seek(position)
                z2=np.uint8(out4[:])
                z2.tofile(f)
                del z2

                # 10.5 Generate the diagnostic plot if diagnosticplot=="yes"
                if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
                    place=(np.floor((I_tStart[m])/(tsamp)) + k)
                    cut=np.uint64(place-1)

                    z[:,place]=out4
    f.close()


    f = open(outputFile, 'rb+')

    # 10.4 Generate narrowband noise
    for m in range(0,np.uint8(N_Occurrences)):

            z1=np.copy(listOfListsN[m])

            mask=np.copy(np.abs(z1))


            maksimum=np.max(mask)
            mask = -1*(mask-maksimum)
            average=np.abs(np.mean(z1))
            mask[mask<average]=1*average
            sigma=np.sqrt(mask) + np.abs(np.mean(z1))
            z1=mask


            NormalizationValue=np.power((np.mean(np.abs(z1))+(np.mean(np.abs(sigma))*2.4)),2)*4
            for k in range(0, np.shape(mask)[1]):

                diff1=np.ceil((fch1-N_FEnd)/np.abs(foff))
                diff2=np.floor((fch1-N_FStart)/np.abs(foff))
                cut=np.uint64(diff1+2)
                numberOfchans=diff2-diff1
                out7=noise_NarrowbandPower(z1[0,k],sigma[0,k],numberOfchans)
                out4=QZ.quantizationOfNarrowbandNoise(N_Magnitude[m],out7,NormalizationValue)

                del out7

                position=np.uint64(258+ np.floor((I_tStart[m])/(tsamp)) + k*nchans +  diff1)
                f.seek(position)
                z2=np.uint8(out4[:])
                z2.tofile(f)
                del z2
                if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
                    placeX=(np.floor((N_tStart[m])/(tsamp)) + k)
                    placeY=diff1
                    z[placeY:(placeY+len(out4)), placeX ] = out4
    f.close()







    f.close()
    print(time.time()-t)

    if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
        z =np.copy(z[:,1::])
        plt.close('all')
        fig = plt.figure()
        ax1 = plt.subplot2grid((6,6), (0,0), colspan=5, rowspan=5)
        ax3 = plt.subplot2grid((6,6), (0,5), colspan=1, rowspan=5)
        ax2 = plt.subplot2grid((6,6), (5,0), colspan=5, rowspan=1)

        pos=np.uint16(cut)
        trace=np.copy(z[pos,:])
        ax2.plot(np.linspace(0, 10, numberOfSamples),trace.T,'k')
        ax2.set_xlim((0, obstime))
        ax2.set_ylim((0, 255))
        plt.tight_layout()

        z=z/255*6.67
        im= ax1.imshow(z, vmin=0, vmax=6.67,aspect='auto',extent=[0,obstime,(fch1+nchans*foff),fch1])
        ax1.set_xlabel('Time(sec)')
        ax1.set_ylabel('Frequency channels (MHz)')
        cbar=plt.colorbar(im ,cax=ax3)
        cbar.set_label('$\sigma$ from expected value', rotation=270, labelpad=20, y=0.5)
        plt.show()





#     infile1 = open('temp.fil', 'rb')
#     infile1.seek(258)
#     x1 = infile1.read(256)#.decode("utf-8")
#     infile2 = open('output.fil', 'rb')
#     infile2.seek(258)
#     x2 = infile2.read(128)#.decode("utf-8")
#     for k in range(0,256):
#         print(struct.unpack('B', x1[k])[0])#,struct.unpack('B', x2[k])[0])

#     infile = open('output.fil', 'rb')
#     outfile = open('output2.fil', 'wb')
#     outfile.write(infile.read(258))
#     outfile.close()
#
#     outfile = open('output2.fil', 'ab')
#     for l in range(0,numberOfSamples):
#         for k in range(0, nchans):
#             infile.seek((258+l+k*numberOfSamples))
#             value=infile.read(1)
#             outfile.write(value)
#     infile.close()
#     outfile.close()





