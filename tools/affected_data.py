#!/usr/bin/env python

"""
*****************************************************
|       Testing RFI Injected Filterbank             |
*****************************************************

Description: Verifies if a filterbank file contains RFI or not using several method

Method 1: Kolmogorov-Smirnov Test (KS-Test), where the underlying distribution of the filterbank is compared to Gaussian Noise with chosen Mean and Std
Method 2: Expected Mean is calculated using the RFI configuration (.yaml) file and compared with the actual mean of the Filterbank. The test is passed if 
            the difference of means is not more than 5% of the mean. 
            Note: User can provide the YAML file in the arguments or the location of the folder of rfi_config files, where the program searches for the 
            YAML file using the RFI_ID from filename of the test-vector.
"""


import argparse
import os

import numpy as np
import yaml
from scipy.stats import kstest
from sigpyproc.readers import FilReader

# pylint: disable=W1202
# pylint: disable=C0301
# pylint: disable=W0611
# pylint: disable=R0914

def fexists(this_file):
    """
    Checks if a file exists. True if yes, else False
    """
    if os.path.isfile(this_file):
        print("File {} found".format(this_file))
        return True
    else:
        print("File {} not found".format(this_file))
        return False

def load_yaml(yamlfile):
    """
    Loads in the YAML file. Checks it's valid YAML.
    Exception raised if not.
    """
    with open(yamlfile, 'r') as this_yaml:
        try:
            pars = yaml.load(this_yaml, Loader=yaml.FullLoader)
        except yaml.parser.ParserError:
            print("{} not valid YAML".format(yamlfile))
            raise Exception("YamlError")
    return pars

def calc_affect(noise_filename,rfi_yaml):
    print('Reading Filterbank file')
    fil = FilReader(noise_filename)
    
    tsamp = fil.header.tsamp
    nchans = fil.header.nchans
    fch1 = fil.header.fch1
    foff = fil.header.foff
    total_samples = nchans*fil.header.nsamples

    print('Reading RFI YAML Configuration File')
    rfi = load_yaml(rfi_yaml)
    affected_cells = 0

    print ('Calculating the number of affected cells')
    for i in rfi:
        tstart = int(np.floor(i['tstart']/tsamp))
        tstop = int(np.floor(i['tstop']/tsamp))
        fstart = int((fch1-i['fstart']+foff*nchans)/foff)
        fstop = int((fch1-i['fstop']+foff*nchans)/foff)
        period = int(np.floor(i['period']/tsamp))

        if period==0:
            affected_cells+= (tstop-tstart)*(fstop-fstart)
        else:
            on_time = period*i['duty']*int((tstop-tstart)/period)
            affected_cells+= on_time*(fstop-fstart)

    print('Total number of affected cells => ',affected_cells)
    print('Fraction of noise file affected => ',round(100*affected_cells/total_samples,3),'%')

def main():
    """
    Main method
    """
    parser = argparse.ArgumentParser(description="Calculating the fraction of Data affected")
    parser.add_argument('noise', help='Path to the noise or rfi file, to get header information',type=str)
    parser.add_argument('--rfi','-r',help='Path the RFI YAML Configuration file',type=str)
    args = parser.parse_args()

    rfi_yaml_file = None
    if args.rfi:
        args.rfi_loc=None           #Setting the RFI_loc to None in order to avoid searching for the file again
        if not fexists(args.rfi):
            exit(1)
        rfi_yaml_file = args.rfi

    #Double chekcing the existence of RFI YAML File
    if not fexists(args.noise):
        exit(1)

    calc_affect(args.noise,rfi_yaml_file)
if __name__=='__main__':
    main()
