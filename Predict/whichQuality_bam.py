####!/usr/bin/python

# samtools view sam/bam | python whichQuality_bam.py --- prints out a value of 33 or 64

import sys

counter = 0
max = 0
min = 127

try:
    for line in sys.stdin:
        if line.startswith('@'): continue
        
        stripped = line.rstrip()
        strippedSplit = stripped.split('\t')
        for letter in strippedSplit[10]:
                if ord(letter) > max:
                    max = ord(letter)
                if ord(letter) < min:
                    min = ord(letter)
            #sys.stdout.write(line)
        counter += 1
        
        if counter == 10000:
            if 30 <= min <= 40 and 70 <= max <= 80:
                sys.stdout.write('33')
            elif 60 <= min <= 75 and 100 <= max <= 110:
                sys.stdout.write('64')
            else:
                sys.stdout.write('unknown')
            sys.exit()
            
except IOError:
    sys.exit('end of file')

