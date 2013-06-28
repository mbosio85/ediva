#!/usr/bin/python
from sys import exit
from sys import argv
import re

###############
# subroutines #
###############
def error(reason):
    if reason == 'fileformat':
        print 'input file does not have the correct format'
    elif reason == 'outfileproblem':
        print 'output file could not be opened'
    exit()

def usage(reason):
    print """usage: extend_bed.py locations.bed extension
    
    This script takes a tab delimited bed file
    chr1	12	25
    chr1	155	899
    chr1	720	1800
    ...
    Overlapping regions will be fused and afterwards
    extended by the 'extension' given.
    File name will get a suffix stating the extension applied.
    
    This script needs at least one file name to work on.
    and a number to extend the probeset [default 0]"""
    if reason == 'number':
        print '\n\n second argument is not a number - set to default = 0 \n'
    elif reason == 'filenotfound':
        print '\n\n first argument is not a valid file name \n'
        exit()
    else:
        exit()

####################
# code starts here #
####################

#initialize variables
chr = 'NA'
chr_old = 'NA'
probestart = 0
probeend = 0
extension = 0
extendedprobestart = 0
extendedprobeend = 0
extendedprobestart_old = 0
extendedprobeend_old = 0

#if the second argument is not given set extension to 0
try:
    int(argv[2])
except IndexError:
    usage('number')
    extension = 0
except ValueError:
    usage('number')
    extension = 0
else:
    extension = int(argv[2])

#open and read input file
try:
    fh = open(argv[1]) #open file
except IndexError:
    usage('filenotfound')
except IOError:
    usage('filenotfound')
else:
    #read first line
    line = fh.readline()

#open output file
try:
    outname = argv[1] + '_plus_' + str(extension)
    outfh = open(outname, 'w')
except IOError:
    error('outfileproblem')

#work through the file
while (line != ''):
    #skip header
    if not line.startswith('chr'):
        line = fh.readline()
        continue
    
    #check line elements
    line_elements = line.split()
    try:
        chr = line_elements[0]
        probestart = int(line_elements[1])
        probeend = int(line_elements[2])
        str(chr)
    except IndexError:
        error('fileformat')
    except TypeError:
        error('fileformat')
    
    #reset _old values if we encounter a new chromosome, but before write old data to outfile
    if chr != chr_old and chr_old != 'NA':
        line_new = chr_old + '\t' + str(extendedprobestart_old) + '\t' + str(extendedprobeend_old) + '\n'
        outfh.write(line_new)
        extendedprobestart_old = 0
        extendedprobeend_old = 0
        chr_old = chr
    
    #extend probe set
    extendedprobestart = probestart - extension
    extendedprobeend = probeend + extension
    if extendedprobestart < 0:
        extendedprobestart = 1
    
    #test for overlap
    if extendedprobeend_old >= extendedprobestart and chr_old == chr:
        extendedprobeend_old = extendedprobeend
    #write to file, if no overlap and memorize current line for next round
    else:
        # skip the first round and memorize first line
        if extendedprobestart_old == 0 and chr_old == 'NA':
            extendedprobestart_old = extendedprobestart
            extendedprobeend_old = extendedprobeend
            chr_old = chr
        #skip chromosome switches    
        elif extendedprobestart_old == 0 and chr_old != 'NA':
            extendedprobestart_old = extendedprobestart
            extendedprobeend_old = extendedprobeend
        #write to file
        else:
            line_new = chr_old + '\t' + str(extendedprobestart_old) + '\t' + str(extendedprobeend_old) + '\n'
            outfh.write(line_new)
            chr_old = chr
            extendedprobestart_old = extendedprobestart
            extendedprobeend_old = extendedprobeend
    
    
    #read next line
    line = fh.readline()

#write out last line
line_new = chr_old + '\t' + str(extendedprobestart_old) + '\t' + str(extendedprobeend_old) + '\n'
outfh.write(line_new)
