About
==============
The code Ersatz.py generates non-stationary Gaussian noise and RFI. 
The output is saved as a SIGPROC filterbank file. The correlation 
length and the amplitude of the non-stationary noise along with the 
specifics for the injected RFI are specified in a .txt file that is 
passed to the function.

SIGPROC is a pulsar search processing software: See the documentation at http://sigproc.sf.net/

Software needed
==================
python

Libraries needed
==================
Scipy
Numpy


How to use
==================
1. Download the stand-alone python function ersatz.py.

2. Create an input text file having the exact layout as shown between 
the starred border (The values may be customised to suit individual scenarios):
*******************************************
Baseline Lambda 2
Baseline Amplitude 4
Broadband Occurrences 0
Broadband t_start 0
Broadband t_end 0
Broadband Magnitude 0
Narrowband Occurrences 0
Narrowband F_start 0
Narrowband F_end 0
Narrowband t_start 0
Narrowband t_end 0
Narrowband Magnitude 0
Periodic Broadband Occurrences 0
Periodic Broadband Period 0
Periodic Broadband Duty cycle 0
Periodic Broadband t_start 0
Periodic Broadband t_end 0
Periodic Broadband Magnitude 0
Periodic Narrowband Occurrences 0
Periodic Narrowband Period 0
Periodic Narrowband Duty cycle 0
Periodic Narrowband F_start 0
Periodic Narrowband F_end 0
Periodic Narrowband t_start 0
Periodic Narrowband t_end 0
Periodic Narrowband Magnitude 0
BandPass ramp-up 0.3
BandPass ramp-down 0.2
BandPass Amplitude 0.5
*******************************************
            
Identifier ||	      Description
========================================================
Baseline 	        Precedes the variables associated with the non-stationary baseline.
Broadband 	        Precedes the variables associated with impulse RFI.
Narrowband 	        Precedes the variables associated with narrowband RFI.
Bandpass 	        Precedes the variables associated with giving the bandpass a shape.
Baseline Lambda         Non-stationary baseline correlation length in seconds.
Baseline Amplitude      Amplitude of the non-stationary noise (0< value ≤ 6).
Occurrences             Number of RFI instances.
t_start 	        The start time of the RFI instances ( s )
t_end                   The end time of the RFI instances  ( s )
Magnitude 	        The magnitude of the noise i.t.o. the number of standard deviations from  
                        the expected value (0< value ≤ 6).
F_start                 First frequency affected by RFI artefact (MHz)
F_end 	                        Ending frequency affected by RFI artefact (MHz)
Period 	                        Period of the periodic noise (s)
Duty cycle 	        Specify how long the RFI should be 'on' as a percentage of its period (0<value<99)
Ramp-up  	        Specify the protion of highest frequency channels to be modulated (0<value<1).
Ramp-down 	        Specify the portion of lowest frequency channels to be modulated (0<value<1).
Amplitude 	        Specifies the height of the slope (0<value<1).


NB. If any of the noise phenomena should be ignored just set the 'Occurrences' variable to zero.
NB. For the bandpass shape to be enforced the flag '"--bandPass" should be set when ersatz.py is calle


4. To see how to run the code, type the following in the terminal/command prompt:
>> python fake_noise.py --info

5. An example of how to run the code then type the following in the terminal/command prompt, the output will be written to the output file called 'output.fil':
>> python ersatz.py --tobs 20.0 --tsamp 64 --fch1 1536 --foff -0.62890625 --nchans 512 --nbits 8 --noiseInput FakeNoiseParameters.txt --stationary "No" --bandPass



