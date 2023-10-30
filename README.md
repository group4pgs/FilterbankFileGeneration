# Filterbank File Generator

## About
The code Ersatz.py generates non-stationary Gaussian noise and RFI. 
The output is saved as a SIGPROC filterbank file. The correlation 
length and the amplitude of the non-stationary noise along with the 
specifics for the injected RFI are specified in a .txt file that is 
passed to the function.

SIGPROC is a pulsar search processing software: See the documentation at [http://sigproc.sf.net/](SIGPROC)

**Dependencies**:
- Numpy
- Scipy

## How to use
Read the file __Readme_For_Ersartz.txt__