About
==============
This is the readme for generating non-Gaussian non-Stationary noise to be used in Michael Keith's release of Duncan Lorimer's
SIGPROC.

SIGPROC is a pulsar search processing software: See the documentation at
http://sigproc.sf.net/

Software needed
==================
python

Libraries needed
==================
Scipy
Numpy
Matplotlib


How to use
==================
1. Download the python project from the github repository at https://github.com/EllieVanH/PulsarDetectionLibrary.git

2. Create an input text file having the exact layout as shown between the starred border (The values may be customised to suit individual scenarios):
***************************
B Lambda 3
I Occurrences 3
I t_start 4,5,6
I t_end 4.5,5.2,6.1
I Magnitude 1,2,6.6
N Occurrences 1
N F_start 1480
N F_end 1500
N t_start 6.7
N t_end 7.0
N Magnitude 3
****************************

where:
’B’ precedes the variables associated with the baseline drift noise
’I’ precedes the variables associated with impulse RFI
’N’ precedes the variables associated with narrowband RFI

Furthermore:
Lambda - Regularity of upward baseline drifts during an observatoin (Lambda should lie between 1 and 10, 
where 1=highly regular,...,10=least regular.)
Occurrences - Number of RFI artefacts in an observation
Magnitude- The Magnitude of the noise i.t.o. the number of standard deviations (std) from the expected value.
The magnitude specified should be between 
0 and 6.67 (i.e. 0 < std < 6:67)
F_start - First frequency affected by RFI artefact (MHz)
F_end - Ending frequency affected by RFI artefact (MHz)
t_start - The start time of the RFI artefact (s)
t_end - The end time of the RFI artefact (s)

NB. If any of the noise phenomena should be ignored just set the ’Occurrences’ variable to zero.

3. In the terminal/command prompt go to the direcotry where the following pyhton files are saved: 
fake_noise.py
noise_BaseLineDrift.py
noise_Narrowband.py
noise_Impulse.py
quantizationOfSignalValues.py

4. To see how to run the code, type the following in the terminal/command prompt:
>python fake_noise.py --info

5. An example of how to run the code with its default values then type the following in the terminal/command prompt, the output will be written to the output file called 'output.fil':
>python fake_noise.py --noiseInput test.txt



