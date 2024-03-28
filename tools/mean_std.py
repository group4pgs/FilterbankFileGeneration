import numpy as np
import matplotlib.pyplot as plt
from sigpyproc.readers import FilReader
import argparse
from rich.pretty import Pretty

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

if end_data==-1:
    end_data = fil.header.nsamples
data = fil.read_block(start_data,end_data)

mean = np.mean(data)
std = np.std(data)
iqr = (np.percentile(data,0.75)-np.percentile(data,0.25))/1.349

print('Mean =',mean)
print('Std =',std)
print('Std (IQR) =',iqr)
print('Median =',np.median(data))

data = data.reshape((data.shape[0]*data.shape[1]))
plt.hist(data,bins=50)
plt.title('Histogram of filterbank data')
plt.grid()
plt.xlabel('Value')
plt.ylabel('Count')
if args.filename == 'show':
    plt.show()
elif args.filename:
    plt.savefig(args.filename)
else:
    plt.savefig('img.png')