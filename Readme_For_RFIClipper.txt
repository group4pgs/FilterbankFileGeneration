About
==============
This is the readme for a function that similtaneously mitigates RFI and removes any varying baseline from SIGPROC filterbank files.

SIGPROC is a pulsar search processing software: See the documentation at http://sigproc.sf.net/

How to compile .cpp file
=========================

g++ -std=c++11 RFI_Clipper_Modified.cpp -o RFI_Clipper_Modified


How to use
==================
1. Download and compile the function.

2. Call the function as follows:
./RFI_Clipper_Modified -inputFile "SNInject.fil" -outputFile "SNInject_normalised.fil"  





