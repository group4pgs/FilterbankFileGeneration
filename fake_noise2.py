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


from noise_BaseLineDrift import noise_BaseLineDriftSmooth
from noise_Impulse import noise_ImpulsePowerPlot
from noise_Narrowband import noise_NarrowbandPowerPlot

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
    "   --nbits,-b          Number of bits (8, 16, 32) (def=8)\n"\
    "   --noiseInput,-n     File name: file containing noise specifications\n"\
    "   --seed,-S           Random seed (def=time())\n"\
    "   --name,-s           Source name for header (def=FAKE)\n"\
    "   --plot,-p           Diagnostic plot (Yes/No) (def=No)\n"\
    "   --header,-H         Write header to output file (Yes/No) (def=Yes)\n"\
    "   --stationary,-y     Produce stationary noise (Yes/No) (def=No)\n"\
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
    seedValue	    = np.uint64(np.random.randint(100000,size=1))
    source_name     = "Fake"
    diagnosticplot  = "No"
    outputFile      = "output.fil"
    header          = "Yes"
    stationary      = "Yes"

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
    parser.add_argument("-b","--nbits", help="Number of bits (8,16,32) (def=8)", type=int)
    parser.add_argument("-n","--noiseInput", help="File name: file containing noise specifications",
                        metavar="FILE", type=lambda x: is_valid_file(parser,x))
    parser.add_argument("-S","--seed", help="Random seed (def=time())", action="store", type=float )
    parser.add_argument("-s","--name", help="Source name for header (def=FAKE)", action="store")
    parser.add_argument("-o","--outFile", help="Output file name (def=output.fil)", action="store")
    parser.add_argument("-p","--plot", help="Diagnostic plot (Yes or No) (def=No)", action="store")
    parser.add_argument("-H","--header", help="Write header to output file (Yes/No) (def=Yes)", action="store")
    parser.add_argument("-y","--stationary", help="Produce stationary noise (Yes/No) (def=Yes)", action="store")

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
    if args.plot:
        diagnosticplot=args.plot
    if args.header:
        header=args.header
    if args.stationary:
        stationary=args.stationary

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
    "   seed        ", seedValue        ,"\n"\
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
    "   Narrowband F_end          ", N_FEnd         ,"MHz\n"\
    "   Narrowband t_start        ", N_tStart       ,"\n"\
    "   Narrowband t_end          ", N_tEnd         ,\
    "\n"
    print('##################################################################')

####################################################################if __name__ == '__main__':######
# 9. Check that all the parameters are kosher
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
# 11. Generate Baseline drift noise profile
###########################################################################
    if (nbits==32):
        nbits=8
        mbits=32
    else:
        mbits=nbits
        
    t=time.time()
    
    numberOfSamples=np.uint64(obstime/(tsamp))

    listOfListsI=[]
    listOfListsN=[]
    z=np.zeros(nchans).T
    diagnosticPlot=[]

    # Generate the smooth baseline drift
    print('Start: Smooth baseline drift profile')
    if (stationary=='No'):
        seedValue=np.uint64(np.random.randint(100000,size=1))
        points = obstime/(lamda*2/10)
        out1=noise_BaseLineDriftSmooth(1, lamda, points , obstime , seedValue)
        x = np.linspace(0, obstime,points)
        xi = np.linspace(0, obstime, numberOfSamples)
        linear = interp1d(x, out1)
        y = linear(xi)
        y = y -np.min(y)

	'''
        seedValue=np.int64(time.time()/10000000)
        out1=noise_BaseLineDriftSmooth(1, lamda, 10000, 100, seedValue)
    
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
        '''
        
        std = np.round(((np.power(2,nbits))/21.34),0)
        PosiveOffset= 3.5*std
        dynamicRangeOfY=std/np.sqrt(nchans)/5*4.0
        y = (y-np.min(y))/(np.max(y)-np.min(y))*dynamicRangeOfY+PosiveOffset-dynamicRangeOfY/2
        Normalize = np.max(y)+(np.sqrt(np.max(y))*1.0*3.1)  
    else:
        seedValue=np.uint64(np.random.randint(100000,size=1))
        std = np.round(((np.power(2,nbits))/21.34),0)
        PosiveOffset= 3.5*std 
        dynamicRangeOfY=std/np.sqrt(nchans)/5*4.0
        y=np.ones(numberOfSamples)*(PosiveOffset) # + 0.5*dynamicRangeOfY)
        Normalize = np.max(y)+(np.sqrt(np.max(y))*1.0*3.1)
        print('End: Baseline drift profile')
###########################################################################
# 12. Generate Impulse noise profile(s)
###########################################################################
    if (np.uint64(I_Occurrences)!=0):
        print('Start: Impulse noise profile(s)')
#        x1=np.arange(-3,3.000072,0.001)
#        y1=np.exp(-np.power(x1,2)/2)/np.sqrt(2*np.pi)/0.4
#        del x1
    for m in range(0,np.uint64(I_Occurrences)):
        
        TimeDuration=np.float64(I_tEnd[m]-I_tStart[m])
        nrOfSamplesI=np.floor(TimeDuration/(tsamp))
        out3=I_Magnitude[m]*np.ones(nrOfSamplesI)#np.sqrt(PosiveOffset)*signal.resample(y1,nrOfSamplesI)
        location=np.uint64(np.round(I_tStart[m]/tsamp))
        y[location:(location+nrOfSamplesI)]=y[location:(location+nrOfSamplesI)]+out3[:]
        del out3
    if (np.uint64(I_Occurrences)!=0):
        print('End: Impulse noise profile(s)')
    
 
    
###########################################################################
# 13. Generate Narrowband noise profile(s)
###########################################################################
    if (np.uint64(N_Occurrences)!=0):
        print('Start: Narrowband noise profile(s)')
#        x1=np.arange(-3,3.000072,0.001)
#        y1=np.exp(-np.power(x1,2)/2)/np.sqrt(2*np.pi)/0.4
#        del x1
        y2=np.multiply(y,np.ones(nchans).reshape((nchans,1)))
        del y
        y=y2
    for n in range(0,np.uint64(N_Occurrences)):
        firstchannel=int((fch1-N_FEnd[n])/np.abs(foff))
        channelsaffected=np.abs(int((N_FEnd[n]-N_FStart[n])/np.abs(foff)))

        x1=np.arange(1,(np.ceil(channelsaffected/2.0+1)))     
        x2=np.power(x1,(-1/np.sqrt(channelsaffected)))
        if np.mod(channelsaffected,2)==0:
            y1=np.concatenate((x2[::-1],x2[:]),axis=1)
        else:
            y1=np.concatenate((x2[:-1],x2[:]),axis=1)

        del x1,x2
        TimeDuration=np.float64(N_tEnd[n]-N_tStart[n])
        nrOfSamplesI=np.floor(TimeDuration/(tsamp))
        out3=N_Magnitude[n]*np.ones(nrOfSamplesI)#*np.sqrt(PosiveOffset)*signal.resample(y1,nrOfSamplesI)
        location=np.uint64(np.round(N_tStart[n]/tsamp))
        for indx in range(0,channelsaffected):
            #y[firstchannel:(firstchannel+channelsaffected),location:(location+nrOfSamplesI)]=y[firstchannel:(firstchannel+channelsaffected),location:(location+nrOfSamplesI)]+out3[:]
            y[(firstchannel+indx),location:(location+nrOfSamplesI)]=y[(firstchannel+indx),location:(location+nrOfSamplesI)]+out3[:]*y1[indx]
        del out3
    if (np.uint64(N_Occurrences)!=0):
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
        sigma = 1*np.sqrt(y)
        # 10.1 Generate Baseline drift noise
#        out=noise_BaseLineDriftPowerPlot(y, sigma, numberOfSamples)
#        out2=QZ.quantizationOfBaseLineSignal(out,nbits)


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
    # Define the function specifying the standard deviation of the baseline drift  noise

    #sigma = 1.0*np.sqrt(y)
    sigma=y
  
   
    if ((header=="Yes") or (header=="yes")):
        f = open(outputFile, 'ab')
    else:
        f = open(outputFile, 'wb')
    
    seedValue=np.uint64(np.random.randint(100000,size=1))
    np.random.seed(seedValue)
    noise= np.random.normal(0,1,(nchans*numberOfSamples*4))
    if (np.uint64(N_Occurrences)!=0):
        for k in range(0, numberOfSamples):
            # 10.1 Generate Baseline drift noise
            z1_noise= np.multiply(np.sqrt(sigma[:,k]),noise[(k*4*nchans):(k*4*nchans+nchans)]) + (y[:,k]).reshape((1,nchans))
            z2_noise= np.multiply(np.sqrt(sigma[:,k]),noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)]) + (y[:,k]).reshape((1,nchans))
            z3_noise= np.multiply(np.sqrt(sigma[:,k]),noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)]) + (y[:,k]).reshape((1,nchans))
            z4_noise= np.multiply(np.sqrt(sigma[:,k]),noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)]) + (y[:,k]).reshape((1,nchans))
            z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
            z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))              
        
            # 10.2 Write the values out to the binary file
            if (mbits==8):
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
        for k in range(0, numberOfSamples):
            # 10.1 Generate Baseline drift noise
            z1_noise= (np.sqrt(sigma[k])*noise[(k*4*nchans):(k*4*nchans+nchans)] + y[k])
            z2_noise= (np.sqrt(sigma[k])*noise[(k*4*nchans+nchans):(k*4*nchans+ 2*nchans)] + y[k])
            z3_noise= (np.sqrt(sigma[k])*noise[(k*4*nchans+ 2*nchans):(k*4*nchans+ 3*nchans)] + y[k])
            z4_noise= (np.sqrt(sigma[k])*noise[(k*4*nchans+ 3*nchans):(k*4*nchans+ 4*nchans)] + y[k])
            z_pow=(np.power(z1_noise,2)+np.power(z2_noise,2)+np.power(z3_noise,2)+np.power(z4_noise,2))
            z_pow=np.float32((z_pow.T)/(4*np.power((Normalize),2))*np.power(2,nbits))              
        
            # 10.2 Write the values out to the binary file
            if (mbits==8):
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
















