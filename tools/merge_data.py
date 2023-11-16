import numpy as np
from sigpyproc.readers import FilReader
import argparse
import os.path
import sys
from sigpyproc.block import FilterbankBlock


psr,rfi,fil,tobs = '','','',0
parser = argparse.ArgumentParser(description='Merges two Filterbanks and writes out an output file')

parser.add_argument('-a','--psr', type=str, help='Input Fitlerbank File 1',required=True)
parser.add_argument('-b','--rfi', type=str, help='Input Filterbank File 2',required=True)
#parser.add_argument('-h','--header',type=int, help='Choose the header of output file to match with File 1 or 2, Type: 1 or 2',default=1)
parser.add_argument('-o','--output',type=str, help='Name of the file to write the data to',default='psr_rfi.fil')
parser.add_argument('-s','--start',type=int, help='Start sample after which data will be merged',default=0)
parser.add_argument('-t','--tobs',type=int, help='Time of observation (tsamp is taken from filterbank)',default=-1)

args = parser.parse_args()
if args.psr:
    psr = args.psr
if args.rfi:
    rfi = args.rfi
if args.output:
    fil = args.output
    
if not os.path.isfile(psr):
    print('Pulsar File does not exist')
    sys.exit()
if not os.path.isfile(rfi):
    print('RFI File doesn not exist')
    sys.exit()

psr = FilReader(psr)
rfi = FilReader(rfi)

if args.output:
    fil = args.output
if args.tobs:
    if args.tobs != -1:
        tobs = psr.header.tobs if psr.header.tobs<rfi.header.tobs else rfi.header.tobs
    else:
        tobs = args.tobs

psr_data = psr.read_block(args.start,int(tobs/psr.header.tsamp))
rfi_data = rfi.read_block(args.start,int(tobs/psr.header.tsamp))

out_data = FilterbankBlock((psr_data+rfi_data)/2,header=psr.header)
print('Data has been written to the file -> ',out_data.to_file(fil))

