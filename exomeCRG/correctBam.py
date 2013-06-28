
###!/soft/bin/python
###!/usr/bin/python

import sys
import argparse
import pysam

def correctChromsome(subfirstread, subsecondread):
    # return directly, if there's no mapping to compare...
    if subfirstread.is_unmapped or subsecondread.is_unmapped:
        #print firstread.qname, ' or ', secondread.qname, ' is unmapped.'
        return(subfirstread, subsecondread)
    
    if subfirstread.qname != subsecondread.qname:
        print subfirstread.qname, ' and ', subsecondread.qname, ' produced a problem.'
    
    # tid&pos is the chromosome&nucleotide location of the read, whereas rnext&pnext is the chromosome&nucleotide location of it's mate
    if subfirstread.tid != subsecondread.tid:
        subfirstread.rnext = subsecondread.tid
        subsecondread.rnext = subfirstread.tid
        
    return(subfirstread, subsecondread)
    # end correctChromsome



parser = argparse.ArgumentParser(description = 'Correct bam files, that are spit out from gem2sam. Input can be sam or bam, Output will always be bam')
parser.add_argument('--input', type=str, required=True, dest='infile', help = 'bam/sam file produced by GEM')
parser.add_argument('--output', type=str, required=True, dest='outfile', help = 'bam file corrected for wrong chromosome entries')

args = parser.parse_args()

try:
    #open file in r (read) b (bam) mode
    #bam_infile = pysam.Samfile(args.infile, "rb")
    bam_infile = pysam.Samfile(args.infile) # could be SAM or BAM file 
except:
    sys.exit('could not open input bam file')

try:
    # open/create a new (b)am file for (w)rinting output and use the same header as in the input file
    bam_outfile = pysam.Samfile( args.outfile, "wb", template=bam_infile )
except:
    sys.exit('could not open output bam file')

counter = 0
linecounter = 0

firstread = pysam.AlignedRead()
secondread = pysam.AlignedRead()



for read in bam_infile.fetch(until_eof=True): # fetch all reads including unmapped reads ###until_eof=True
    
    ## produce a shortened file for testing
    #if linecounter == 1000000:
    #    break 
    
    # read even and uneven read lines
    if linecounter % 2 == 0:
        firstread = read
        
    if linecounter % 2 == 1:
        secondread = read
        
        # try to fix off paires
        if firstread.qname != secondread.qname:
            # leave one out and do NOT increase linecounter so that the second read is read again.
            print firstread , ' was dropped.'
            firstread = secondread
            continue
        
        # every second read is the read partner
        (firstread, secondread) = correctChromsome(firstread, secondread)
        
        # write corrected read information to output file
        bam_outfile.write(firstread)
        bam_outfile.write(secondread)

        # re-set variables
        firstread = pysam.AlignedRead()
        secondread = pysam.AlignedRead()
        
    linecounter += 1
    
    # end reading file

#finalize stuff
print linecounter, 'lines read and processed'

bam_infile.close()
bam_outfile.close()



#@SQ	SN:1	LN:249250621
#@SQ	SN:2	LN:243199373
#@SQ	SN:3	LN:198022430
#@SQ	SN:4	LN:191154276
#@SQ	SN:5	LN:180915260
#@SQ	SN:6	LN:171115067
#@SQ	SN:7	LN:159138663
#@SQ	SN:8	LN:146364022
#@SQ	SN:9	LN:141213431
#@SQ	SN:10	LN:135534747
#@SQ	SN:11	LN:135006516
#@SQ	SN:12	LN:133851895
#@SQ	SN:13	LN:115169878
#@SQ	SN:14	LN:107349540
#@SQ	SN:15	LN:102531392
#@SQ	SN:16	LN:90354753
#@SQ	SN:17	LN:81195210
#@SQ	SN:18	LN:78077248
#@SQ	SN:19	LN:59128983
#@SQ	SN:20	LN:63025520
#@SQ	SN:21	LN:48129895
#@SQ	SN:22	LN:51304566
#@SQ	SN:X	LN:155270560
#@SQ	SN:Y	LN:59373566
#@SQ	SN:MT	LN:16571
#@RG	ID:VH026	SM:VH026
#@PG	ID:bwa	PN:bwa	VN:0.5.9-r26-dev
#HWI-ST227:243:C1D97ACXX:4:2305:6497:14539	163	1	10004	17	1S99M	=	10112	208	CCCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCTACCCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTCACCCTCACCCTCACCC@@CFFFFFHGHGDIIGGIIJJJIIHEDHGEHHGGGFE?FHCH9;FD=BFG=FCDGHGGFCHHD?6?C>C?A;;;=;;<CA?<ABD<1<22?8(8A<?###	RG:Z:VH026	XT:A:M	NM:i:6	SM:i:17	AM:i:17	XM:i:6	XO:i:0	XG:i:0	MD:Z:35A5A5A34A5A5A4	XA:Z:12,-133841511,100M,6;X,-155260270,100M,6;X,-155260276,100M,6;12,-133841517,100M,6;X,-155260264,100M,6;
#HWI-ST227:243:C1D97ACXX:4:1306:20035:86291	99	1	10010	0	100M	=	10125	215	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCA@@@FFFDF?FHDBFGIBGGIIIIIEIJJGGHFEH;GHGJIGIIJGJJJ(=FGHEGIHGGIEEEADDB;?@>;;=??;3=?<A?CB###############	RG:Z:VH026	XT:A:U	NM:i:1	SM:i:0	AM:i:0	X0:i:1	XM:i:1	XO:i:0	XG:i:0	MD:Z:93T6
#HWI-ST227:243:C1D97ACXX:4:1307:12454:134599	99	1	10010	29	100M	=	10112	202	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCA@BCFFFFFHGHHHJJJJJJIJJIJJJJIFEIIJJJGIEBDGBHIIGIEEGIGGHIIJIICHF=;C@@BD>>>A=A;==?BA@D<B<BD<(2888<A?BD#	RG:Z:VH026	XT:A:U	NM:i:0	SM:i:29	AM:i:29	X0:i:1	XM:i:0	XO:i:0	XG:i:0	MD:Z:100
