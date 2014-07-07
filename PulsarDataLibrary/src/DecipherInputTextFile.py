#!/usr/bin/python
'''
@author: Elmarie
'''

def DiagnosticPlot():
# Read from a text file the sources of noise
    inputFile='test.txt'
    filename = inputFile
    with open(filename, 'rt') as f:
    #         reader = csv.reader(f)
        for row in f:
            if '#Hello' in row:
                print('sukses')
#             if row == ('#Hello'):
#                 print (4)
#             elif row == 1 or n == 9 or n == 4:
#                 print "n is a perfect square\n"
#             elif row == 2:
#                 print "n is an even number\n"
#             elif  row == 3 or n == 5 or n == 7:
#                 print "n is a prime number\n"
            
    return 0


DiagnosticPlot()
