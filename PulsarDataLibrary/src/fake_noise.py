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
import os
import re
from scipy import signal
import time
from scipy.interpolate import interp1d
###########################################################################
# 2. Include sub-functions.
###########################################################################


from noise_BaseLineDrift import noise_BaseLineDriftSmooth,noise_BaseLineDriftPower, noise_BaseLineDriftPowerPlot
from noise_Impulse import noise_ImpulseSmooth,noise_ImpulsePower,noise_ImpulsePowerPlot
from noise_Narrowband import noise_NarrowbandSmooth,noise_NarrowbandPower,noise_NarrowbandPowerPlot
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
    "   --nbits,-b          Number of bits (8, 16) (def=8)\n"\
    "   --noiseInput,-n     File name: file containing noise specifications\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "   --plot,-p           Diagnostic plot (Yes/No) (def=No)\n"\
    "   --header,-H         Write header to output file (Yes/No) (def=Yes)\n"\
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
    global lamda, I_Occurrences, I_tStart, I_tEnd, I_Magnitude, N_Occurrences, N_FStart, N_FEnd, N_tStart, N_tEnd, N_Magnitude, filename

# Read from a text file the sources of noise
    filename = inputFile
    with open(filename, 'rt') as f:
        for row in f:
            if ('>B' and 'Lambda') in row:
                lamda = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('Occurrence' in row)):
                I_Occurrences = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('t_start' in row)):
                I_tStart = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('t_end' in row)):
                I_tEnd = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('I' in row) and ('Magnitude' in row)):
                I_Magnitude = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('Occurrence' in row)):
                N_Occurrences = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('F_start' in row)):
                N_FStart = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('F_end' in row)):
                N_FEnd = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('t_start' in row)):
                N_tStart = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('t_end' in row)):
                N_tEnd = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))
            elif (('N' in row) and ('Magnitude' in row)):
                N_Magnitude = np.float64(re.findall(r"[-+]?\d*\.\d+|\d+", row))

    return 0

###########################################################################
# 6. Start of the main function: fake_noise.py
###########################################################################

if __name__ == '__main__':

###########################################################################
# 6. Set the default values for generating fake_noise.
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
    seed            = np.random.seed()  #Sets the seed value for the number generator by using the current time
    source_name     = "Fake"
    diagnosticplot  = "Yes"
    outputFile      = "output.fil"
    header          = "No"

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
    cut             = 0
    positionHeader  = 0

###########################################################################
# 7.Parse the arguments passed to fake_noise from the terminal.
###########################################################################


    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--info", help="Information regarding the function",action="store_true")
    parser.add_argument("-T","--tobs", help="Total observation time, s (def=270)", action="store", type=float )
    parser.add_argument("-t","--tsamp", help="Sampling interval, us (def=64)", action="store", type=float )
    parser.add_argument("-m","--mjd", help="MJD of first sample (def=56000.0)", action="store", type=float )
    parser.add_argument("-F","--fch1", help="Frequency of channel 1, MHz (def=1581.804688)", action="store", type=float )
    parser.add_argument("-f","--foff", help=" Channel bandwidth, MHz (def=-0.390625)", action="store", type=float )
    parser.add_argument("-c","--nchans", help="Output number of channels (def=1024)", type=int)
    parser.add_argument("-b","--nbits", help="Number of bits (8,16) (def=8)", type=int)
    parser.add_argument("-n","--noiseInput", help="File name: file containing noise specifications",
                        metavar="FILE", type=lambda x: is_valid_file(parser,x))
    parser.add_argument("-S","--seed", help="Random seed (def=time())", action="store", type=float )
    parser.add_argument("-s","--name", help="Source name for header (def=FAKE)", action="store")
    parser.add_argument("-o","--outFile", help="Output file name (def=output.fil)", action="store")
    parser.add_argument("-p","--plot", help="Diagnostic plot (Yes or No) (def=No)", action="store")
    parser.add_argument("-H","--header", help="Write header to output file (Yes/No) (def=Yes)", action="store")

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
        seed=args.seed
    if args.name:
        source_name=args.name
    if args.outFile:
        outputFile=args.outFile
    if args.plot:
        diagnosticplot=args.plot
    if args.header:
        header=args.plot

    tsamp=tsamp*1e-06
###########################################################################
# 8. Generate a summary of the parameters of the observation to be printed to the screen
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
    "   seed        ", seed             ,"\n"\
    "   name        ", source_name      ,"\n"\
    "   plot        ", diagnosticplot   ,"\n"\
    "   header      ", header           ,"\n"\
    "   outputFile  ", outputFile           ,"\n"\
    "\n"
    print('##################################################################')
    print('          Parameter values from the input file:')
    print('##################################################################')
    print"\n"\
    "   Lambda                    ", lamda          ,"\n\n"\
    "   Impulse Occurrences       ", I_Occurrences  ,"\n"\
    "   Impulse t_start           ", I_tStart       ,"\n"\
    "   Impulse t_end             ", I_tEnd         ,"\n\n"\
    "   Narrowband Occurrences    ", N_Occurrences  ,"\n"\
    "   Narrowband F_start        ", N_FStart       ,"MHz\n"\
    "   Narrowband F_start        ", N_FEnd         ,"MHz\n"\
    "   Narrowband t_start        ", N_tStart       ,"\n"\
    "   Narrowband t_end          ", N_tEnd         ,\
    "\n"
    print('##################################################################')

####################################################################if __name__ == '__main__':######
# 9. Check that all the parameters are kosher
###########################################################################
    if (np.any(I_tStart>=obstime) or np.any(I_tEnd>=obstime)):
        clear = lambda: os.system('cls')
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
    if (np.any(N_tStart>=obstime) or np.any(N_tEnd>=obstime)):
        clear = lambda: os.system('cls')
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
        clear = lambda: os.system('cls')
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


###########################################################################
# 10. Write the Header part of the output file
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
            f.write(struct.pack('<I', 6))
            unicode_array=array('b',b'signed')
            unicode_array.tofile(f)
            f.write(struct.pack('B', 1))
            f.write(struct.pack('<I', 10))
            unicode_array=array('b',b'HEADER_END')
            unicode_array.tofile(f)
            f.close()
            positionHeader = os.path.getsize(outputFile)

            print("Finished writing Header to binary file")



##########################################################################
# 11. Generate Baseline drift noise profile
###########################################################################

    t=time.time()

    numberOfSamples=np.uint64(obstime/(tsamp))

    listOfListsI=[]
    listOfListsN=[]
    z=np.zeros(nchans).T
    averagePerSampleChannel=[]
    diagnosticPlot=[]

    # Generate the smooth baseline drift
    print
    print('Start: Smooth baseline drift profile')
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
    print('End: Baseline drift profile')
    print
###########################################################################
# 12. Generate Impulse noise profile(s)
###########################################################################
    if (np.uint8(I_Occurrences)!=0):
        print('Start: Impulse noise profile(s)')
    for m in range(0,np.uint8(I_Occurrences)):
        out3=[]
        TimeDuration=np.float64(I_tEnd[m]-I_tStart[m])
        nrOfSamplesI=np.floor(TimeDuration/(tsamp))

        seedValueForImpulseNoise=np.random.seed()
        out3.append(noise_ImpulseSmooth(1, np.uint16(nrOfSamplesI),TimeDuration, seedValueForImpulseNoise))
        listOfListsI.append(out3[:])
        del out3
    if (np.uint8(I_Occurrences)!=0):
        print('End: Impulse noise profile(s)')
    
###########################################################################
# 13. Generate Narrowband noise profile(s)
###########################################################################
    if (np.uint8(N_Occurrences)!=0):
        print('Start: Narrowband noise profile(s)')
    for n in range(0,np.uint8(N_Occurrences)):
        out3=[]
        TimeDuration=(N_tEnd[n]-N_tStart[n])
        nrOfSamplesN=np.floor(TimeDuration/(tsamp))
        seedValueForNarrowbandNoise=np.random.seed()
        out3.append(noise_NarrowbandSmooth(1, np.uint16(nrOfSamplesN),TimeDuration, seedValueForNarrowbandNoise))

        listOfListsN.append(out3[:])
        del out3
    if (np.uint8(N_Occurrences)!=0):
        print("End: Narrowband noise profile(s)")
    print

    c=[]
###########################################################################
# 14. DIAGNOSTIC PLOT
###########################################################################

    if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):

        # Define the function specifying the standard deviation of the noise
#        y = (y-np.min(y))/(np.max(y)-np.min(y))*24+72
#        sigma = 2.4*np.sqrt(y)
        std = np.round(((np.power(2,nbits))/10.67),0)
        PosiveOffset= 3.5*std
        y = (y-np.min(y))/(np.max(y)-np.min(y))*std+PosiveOffset
        sigma = 2.2*np.sqrt(y)
        # 10.1 Generate Baseline drift noise
#        out=noise_BaseLineDriftPowerPlot(y, sigma, numberOfSamples)
#        out2=QZ.quantizationOfBaseLineSignal(out,nbits)
        t1 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
        t2 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
        t3 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
        t4 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
        out2= (t1+t2+t3+t4)/(4*np.power(156,2))*256
                
        t5 = 4*np.power((y+0.15*np.sqrt(y)),2)/(4*np.power(156,2))*256
        t6 = 96*np.ones(1,numberOfSamples)+np.multiply(np.random.normal(0,1,numberOfSamples),24)
        plt.plot(t6,'b')
        plt.plot(out2,'k')
        plt.plot(t5,"white")
        plt.title('Power versus time with $\sigma_{noise} = 2.2\sqrt{mean}$')
        plt.show()
        
        outImage=[]
        outImage=np.copy(out2)

        for m in range(0,np.uint8(I_Occurrences)):
            out7=noise_ImpulsePowerPlot(listOfListsI[m],np.array(np.shape(listOfListsI[m][0])))
            out4=QZ.quantizationOfImpulseNoisePlot(I_Magnitude[m],out7,nbits)
            del out7
            outImage[np.floor((I_tStart[m])/(tsamp)):(np.floor((I_tStart[m])/(tsamp))+np.array(np.shape(listOfListsI[m][0])))]=out4[:,0]
            del out4

        toets=np.copy(outImage)

        for k in range(0, nchans):
            # 10.4 Generate narrowband noise
            for n in range(0,np.uint8(N_Occurrences)):
                diff1=fch1-N_FEnd[n]
                diff2=fch1-N_FStart[n]
                if ((np.ceil(diff1/np.abs(foff)) <= k) and ( k <= np.floor(diff2/np.abs(foff)))):
                    out8=noise_NarrowbandPowerPlot(listOfListsN[n],np.array(np.shape(listOfListsN[n][0])))
                    out9=QZ.quantizationOfNarrowbandNoisePlot(N_Magnitude[n],out8,nbits)
                    del out8
                    outImage[np.floor((N_tStart[n])/(tsamp)):(np.floor((N_tStart[n])/(tsamp))+np.array(np.shape(listOfListsN[n][0])))]=out9[:,0]
                    del out9
            c=signal.decimate((outImage[:].T),2)
            diagnosticPlot=np.hstack((diagnosticPlot,c))
            del outImage
            outImage=np.copy(toets)

        plt.close('all')
        fig = plt.figure()
        ax1 = plt.subplot2grid((6,6), (0,0), colspan=5, rowspan=5)
        ax3 = plt.subplot2grid((6,6), (0,5), colspan=1, rowspan=5)
        ax2 = plt.subplot2grid((6,6), (5,0), colspan=5, rowspan=1)
        width=np.uint64(len(diagnosticPlot)/nchans)
        diagnosticPlot=np.reshape(diagnosticPlot, (nchans,width))

        pos=np.uint16(cut)
        trace=np.copy(diagnosticPlot[pos,:])
        ax2.plot(np.linspace(0, obstime, width ),trace,'k')
        ax2.set_xlim((0, obstime))
        ax2.set_ylim((0, (np.power(2,nbits)-1)))
        plt.tight_layout()

        #diagnosticPlot=diagnosticPlot/(np.power(2,nbits)-1)*6.67
        im= ax1.imshow(diagnosticPlot, vmin=0, vmax=255,aspect='auto',extent=[0,obstime,(fch1+nchans*foff),fch1])
        ax1.set_xlabel('Time(sec)')
        ax1.set_ylabel('Frequency channels (MHz)')
        cbar=plt.colorbar(im ,cax=ax3)
        cbar.set_label('$\sigma$ from expected value', rotation=270, labelpad=20, y=0.5)
#     plt.draw()
        plt.show(block=False)


###########################################################################
# 15. Write the Baseline Drift noise to the output file.
###########################################################################
    # Determine how much free memory is available:
    

    # Define the function specifying the standard deviation of the baseline drift  noise
    std = np.round(((np.power(2,nbits))/10.67),0)
    PosiveOffset= 3.5*std
    y = (y-np.min(y))/(np.max(y)-np.min(y))*std+PosiveOffset
    sigma = 2.2*np.sqrt(y)

#    t1 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
#    t2 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
#    t3 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
#    t4 = np.power((y + np.multiply(np.random.normal(0,1,numberOfSamples),sigma)),2)
#    t= (t1+t2+t3+t4)/(4*np.power(156,2))*256
#    t7 = np.array(t, dtype=np.int8)
#    t5 = 4*np.power((y+0.3*np.sqrt(y)),2)/(4*np.power(156,2))*256
#    t6 = 96*np.ones(1,numberOfSamples)+np.multiply(np.random.normal(0,1,numberOfSamples),24)
#    plt.plot(t6,'b')
#    plt.plot(t,'k')
#    plt.plot(t5,'r')

    if ((header=="Yes") or (header=="yes")):
        f = open(outputFile, 'ab')
    else:
        f = open(outputFile, 'wb')
    
    np.random.seed()
    noise= np.random.normal(0,1,(nchans*numberOfSamples*4))
    Normalize = np.max(y)+(np.sqrt(np.max(y))*2.2*3) 

    for k in range(0, numberOfSamples):

        # 10.1 Generate Baseline drift noise
        z1_noise= sigma[k]*noise[(k*4*nchans):(k*4*nchans+nchans)] + y[k]
        z2_noise= sigma[k]*noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)] + y[k]
        z3_noise= sigma[k]*noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)] + y[k]
        z4_noise= sigma[k]*noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)] + y[k]
        z_pow=np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2)
        z_pow=(z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits)
       
#        out=noise_BaseLineDriftPower(y[k], sigma[k], nchans)
#        out2=QZ.quantizationOfBaseLineSignal(out,nbits)
        averagePerSampleChannel.append(np.mean(z_pow))

        # 10.2 Write the values out to the binary file
        if (nbits==8):
            z2=np.uint8(z_pow)
        else:
            z2=np.uint16(z_pow)
        z2.tofile(f)

        del z2
    f.close()
    #print(z_pow)

    print("Finished generating and writing BaselineDrift to output file.")


###########################################################################
# 16. Write the Impulse noise to the output file.
###########################################################################

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
                out4=QZ.quantizationOfImpulseNoise(I_Magnitude[m],out7,NormalizationValue,nbits)

                out4[out4[:]<averagePerSampleChannel[np.uint64(np.floor((I_tStart[m])/(tsamp)) + k)]]=averagePerSampleChannel[np.uint64(np.floor((I_tStart[m])/(tsamp)) + k)]
                del out7

                position=np.uint64(positionHeader + (np.floor((I_tStart[m])/(tsamp)) + k)*nchans )
                f.seek(position)
                if (nbits==8):
                    z2=np.uint8(out4[:])
                else:
                    z2=np.uint16(out4[:])
                z2.tofile(f)
                del z2

    f.close()

###########################################################################
# 17. Write the Narrowband noise to the output file.
###########################################################################

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

                diff1=np.ceil((fch1-N_FEnd[m])/np.abs(foff))
                diff2=np.floor((fch1-N_FStart[m])/np.abs(foff))
                cut=np.uint64(900)
                numberOfchans=diff2-diff1
                out7=noise_NarrowbandPower(z1[0,k],sigma[0,k],numberOfchans)
                out4=QZ.quantizationOfNarrowbandNoise(N_Magnitude[m],out7,NormalizationValue,nbits)

                del out7
                position=np.uint64(positionHeader+ (np.floor((N_tStart[m])/(tsamp)) + k)*nchans +  diff1)
                f.seek(position)
                if (nbits==8):
                    z2=np.uint8(out4[:])
                else:
                    z2=np.uint16(out4[:])
                z2.tofile(f)
                del z2
#                 if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
#                     placeX=(np.floor((N_tStart[m])/(tsamp)) + k)
#                     placeY=diff1
#                     z[placeY:(placeY+len(out4)), placeX ] = out4
    f.close()
    plt.show()
    print(time.time()-t)

#     if ((diagnosticplot=="Yes") or (diagnosticplot=="yes")):
#         z =np.copy(z[:,1::])
#         plt.close('all')
#         fig = plt.figure()
#         ax1 = plt.subplot2grid((6,6), (0,0), colspan=5, rowspan=5)
#         ax3 = plt.subplot2grid((6,6), (0,5), colspan=1, rowspan=5)
#         ax2 = plt.subplot2grid((6,6), (5,0), colspan=5, rowspan=1)
#
#         pos=np.uint16(cut)
#         trace=np.copy(z[pos,:])
#         ax2.plot(np.linspace(0, obstime, numberOfSamples),trace.T,'k')
#         ax2.set_xlim((0, obstime))
#         ax2.set_ylim((0, 255))
#         plt.tight_layout()
#
#         z=z/255*6.67
#         im= ax1.imshow(z, vmin=0, vmax=6.67,aspect='auto',extent=[0,obstime,(fch1+nchans*foff),fch1])
#         ax1.set_xlabel('Time(sec)')
#         ax1.set_ylabel('Frequency channels (MHz)')
#         cbar=plt.colorbar(im ,cax=ax3)
#         cbar.set_label('$\sigma$ from expected value', rotation=270, labelpad=20, y=0.5)
#         plt.show()




#
#     infile1 = open('temp.fil', 'rb')
#     infile1.seek(0)
#     x1 = infile1.read(258)#.decode("utf-8")
#     infile2 = open('BLDLambda5.fil', 'rb')
#     infile2.seek(0)
#     x2 = infile2.read(258)#.decode("utf-8")
#     for k in range(0,258):
#         print(struct.unpack('B', x1[k])[0],struct.unpack('B', x2[k])[0])

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





