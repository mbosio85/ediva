#!/usr/bin/python
from sys import exit
from sys import argv
import re
import argparse

parser = argparse.ArgumentParser(description='Fuse target regions given in an Exome Target file and extend by the given length. If target regions overlap after extension, they will get fused!')

parser.add_argument('--infile', type=argparse.FileType('r'), required=True, dest='infile', help="The Bed file to work on.")
parser.add_argument('--outfile', type=argparse.FileType('w'), required=True, dest='outfile', help='The name of the output file.')
parser.add_argument('--extension', type=int, required=False, dest='extension', default=0, help='The amount of nucleotides by which the probe regions should get extended. [default=0]')

args = parser.parse_args()

####################
# code starts here #
####################

#initialize variables
chr = 'NA'
chr_old = 'NA'
probestart = 0
probeend = 0
extension = args.extension
extendedprobestart = 0
extendedprobeend = 0
extendedprobestart_old = 0
extendedprobeend_old = 0

## write header --- no header defined in bed file format
#header = '\t'.join(["#chr","start","end"])
#args.outfile.write(header + "\n")

#work through the file
for line in args.infile:
    
    #header (should not start with chr or any number)
    if not line.startswith('chr') and not line[:1].isdigit():
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
        args.outfile.write(line_new)
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
    #elif extendedprobeend_old != 0:
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
            args.outfile.write(line_new)
            chr_old = chr
            extendedprobestart_old = extendedprobestart
            extendedprobeend_old = extendedprobeend

else:
    #write out last line
    line_new = chr_old + '\t' + str(extendedprobestart_old) + '\t' + str(extendedprobeend_old) + '\n'
    args.outfile.write(line_new)
