import numpy as np
import matplotlib.pyplot as plt
from rich.pretty import Pretty
from sigpyproc.readers import FilReader
import argparse
from pathlib import Path
from sys import exit

fil = ''
period, dm = 0.0,0.0

parser = argparse.ArgumentParser(description='Displays pulse profiles and trying for SNR too from Filterbank data')

parser.add_argument('input', type=str, help='Input Fitlerbank File')
parser.add_argument('-p','--period',type=float, help='Set the period for folding the pulse (s)',required=True)
parser.add_argument('-d','--dm',type=float,help='Set the DM to dedisperse the time-series (kpc/cm^3)',required=True)
parser.add_argument('-n','--normalise',type=bool,default=False)
parser.add_argument('-f','--filename',type=str, help='Specific filename to save it in (default=img.png)')

args = parser.parse_args()
if args.input:
    fil = args.input

if not (Path(fil).is_file()):
    print('ERROR: File not found!')
    exit()
if args.period:
    period = args.period
if args.dm:
    dm = args.dm

fil = FilReader(fil)
print('**************\n',
      'filename\t:',fil.header.filename,'\n',
      'f_stop\t\t:',fil.header.fch1,'MHz\n',
      'f_res\t\t:',fil.header.foff,'MHz\n',
      'nchans\t\t:',fil.header.nchans,'\n',
      'nbits\t\t:',fil.header.nbits,'\n',
      'tsample\t:',fil.header.tsamp*1e6,'us\n'
      '**************\n',
      'Parameters chosen to fold\n',
      'DM\t:',dm,'kpc/cm^3\n',
      'Period\t:',period,'s\n****************'
      )
fil_fold = fil.fold(period=period, dm=dm, nints = 1, nbands=1, nbins=128)
#fil_fold = fil_fold.centre()

plt.plot(np.linspace(0,1,fil_fold.shape[2]),fil_fold[0][0])
plt.title("Folded Pulse Profile")
plt.grid()
if args.filename == 'show':
    plt.show()
elif args.filename:
    plt.savefig(args.filename)
else:
    plt.savefig('img.png')