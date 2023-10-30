import numpy as np
import matplotlib.pyplot as plt
from rich.pretty import Pretty
from sigpyproc.readers import FilReader
from sys import argv
import argparse

filname = ''
end_data = -1
start_data = 0

parser = argparse.ArgumentParser(description='Plotting the filterbank')

parser.add_argument('-i','--input', type=str, help='Input Fitlerbank File',required=True)
parser.add_argument('-s','--start', type=int, help='Start of the timeseries',default=0)
parser.add_argument('-e','--end',type=int, help='End of timeseries (sample number [-1 for entire file])',default=-1)
args = parser.parse_args()

if args.input:
    filname = args.input
if args.start:
    start_data = args.start
if args.end:
    end_data = args.end

fil = FilReader(filname)
Pretty(fil.header)

if end_data==-1:
    end_data = fil.header.nsamples
data = fil.read_block(start_data,end_data)

plt.imshow(data,aspect='auto')
plt.title("Filterbank Data")
plt.show()