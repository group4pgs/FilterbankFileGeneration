import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.readers import FilReader
import argparse

filname = ''
end_data = -1
start_data = 0

parser = argparse.ArgumentParser(description='Plotting the filterbank')

parser.add_argument('input', type=str, help='Input Fitlerbank File')
parser.add_argument('-s','--start', type=int, help='Start of the timeseries',default=0)
parser.add_argument('-e','--end',type=int, help='End of timeseries (sample number [-1 for entire file])',default=-1)
parser.add_argument('-d','--display',help='Setting this flag will show the image, and not print it',action='store_true')
parser.add_argument('-f','--filename',type=str, help='Specific filename to save it in (default=img.png)')
args = parser.parse_args()

if args.input:
    filname = args.input
if args.start:
    start_data = args.start
if args.end:
    end_data = args.end

fil = FilReader(filname)
y_plt = np.linspace(0,fil.header.nchans,6)
y_ticks = np.round(np.linspace(fil.header.fch1+(fil.header.nchans*fil.header.foff),fil.header.fch1,6),0)

if end_data==-1:
    end_data = fil.header.nsamples
data = fil.read_block(start_data,end_data)

time = end_data*(fil.header.tsamp)
x_ticks = np.round(np.linspace(0,time,10),2)

plt.matshow(data,aspect='auto')
plt.yticks(y_plt,y_ticks[::-1])
plt.xticks(np.linspace(0,end_data,10),x_ticks)
plt.xlabel("Time $(s)$")
plt.ylabel("Frequency $(MHz)$")
plt.title("Filterbank Data")
if args.filename == 'show':
    plt.show()
elif args.filename:
    plt.savefig(args.filename)
else:
    plt.savefig('img.png')