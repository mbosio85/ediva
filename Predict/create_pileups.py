#! /usr/bin/python

import argparse
import os
#import pprint

parser = argparse.ArgumentParser(description = 'Create a merged vcf file and accordingly read bam files to gather coverage information at the variant positions.')

parser.add_argument('--vcffile', type=str, dest='infile', required=True, action='append', help='VCF file(s) to merge and investigate in the format --vcffile SampleID:file . Can be called multiple times. [required]')
parser.add_argument('--bamfile', type=str, dest='bamfile', required=False, action='append', help='BAM file(s) to retrieve coverage information. Can be called multiple times and should correlate to the given VCF files.')
parser.add_argument('--outfolder', type=str, dest='outfolder', required=True, help='Folder to write the script to run and mpileup result to. (some temporary files are also written to that folder)')
parser.add_argument('--qsub_name', type=str, dest='qsub_name', required=True, help='Name of the script you should start once created.' )

args = parser.parse_args()

## DEBUG
#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(args.infile)
## /DEBUG

if not os.path.exists(args.outfolder):
    print "Creating output folder %s" % args.outfolder
    os.makedirs(args.outfolder)

# element.split(':') for element in args.infile -- creates a list of tuples [ [ id, path ] [ id, path] ]
# zip ( *list ) splits the whole thing into two lists
ids, paths = zip( *[ element.split(':') for element in args.infile ])

# create header elements
id_elements = list()
for element in ids:
    id_elements.append(element + "count")
    id_elements.append(element + "obs")
    id_elements.append(element + "qual")
header = '\\t'.join(id_elements)


# will contain vcf --variant vcf --variant vcf ..., i.e. a --variant call in GATK is still necessary
vcf_files = ' --variant '.join(paths)

try:
    bam_file_list = '\n'.join(args.bamfile)
    
    # save the locations of the bam files to a temporary file
    try:
        bam_list = args.outfolder + '/bam.list'
        FH = open(bam_list, 'w')
    except:
        print "Could not open output folder for writing."
        exit(1)

    FH.write(bam_file_list)
    FH.close()
except:
    print "No bam files given. Only producing a merge vcf script."

# qsub doesn't like jobs starting with a digit
if args.qsub_name[0].isdigit():
    qsub_name = 's' + args.qsub_name
else:
    qsub_name = args.qsub_name


# basic merge vcf script
out = """

#$ -N %s
#$ -e %s
#$ -o %s

GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
SAMTOOLS=$(which samtools)

REF=/users/GD/resource/human/hg19/bwa7/hg19.fasta
OUTF=%s

# merge vcf files
java -Xmx2g -jar $GATK -T CombineVariants -R $REF --variant %s -o $OUTF/combined.variants.vcf

""" % (qsub_name, args.outfolder, args.outfolder, args.outfolder, vcf_files)

# if enough (vcf=bam) bam files were given:
if args.bamfile is None:
    #already checked and printed a message
    pass

# if there's an uneven number (vcf=bam) of files given:
elif not len(args.bamfile) == len(args.infile):
    print """ You did not provide the same amount of bam files and vcf files.
    Just creating a merge vcf script.
    """

# if everything is fine add the mpileup part of the script
elif len(args.bamfile) == len(args.infile):
    out = out + """

# pileup call
sed -e 's/chr//' $OUTF/combined.variants.vcf | awk '{OFS=\"\\t\"; if (!/^#/){print $1,$2}}'  > $OUTF/variant.position.bed

# header
echo -e \"chr\\tstart\\tref\\t%s\" > $OUTF/variant.positions.mpileup

# pileup
$SAMTOOLS mpileup -l $OUTF/variant.position.bed -f $REF -b $OUTF/bam.list >> $OUTF/variant.positions.mpileup



""" % (header)

# output
try:
    script_location = args.outfolder + args.qsub_name
    SCRIPT = open( script_location, 'w')
except:
    print "Could not open script file for writing."
    exit(1)

SCRIPT.write(out)
SCRIPT.close()

exit(0)