import pprint
import argparse
import os

pp = pprint.PrettyPrinter(indent = 5)

parser = argparse.ArgumentParser(description = 'Create a file containing all necessary information to run the eDiVa prioritization.')

parser.add_argument('--outfile',  type=argparse.FileType('w'), dest='config',  required=True, help='The output file will contain all sample names, the affection status and the file locations of the vcf and bam file.')
parser.add_argument('--sampleinfo', type=str, dest='sampleinfo', required=False, action='append', help="""Instead of going through the wizzard, you can also provide all information in one line.
Please follow the pattern --sampleinfo ID:affected:vcf:bam
affected is either 0 (unaffected) or 1 (affected)
vcf is the full vcf file path
bam is the full bam file path [optional]

This option can be called multiple times.
If GATK multi-sample calls are supplied, just enter the same vcf file name for each sample and toggle the --multisample option creating the pipeline.

""")

args = parser.parse_args()

# write header
args.config.write('\t'.join(['ID', 'status', 'vcf', 'bam']))
args.config.write('\n')

if args.sampleinfo:
    for line in args.sampleinfo:
        splitline = line.split(':')
        name = splitline[0]
        affect = splitline[1]
        try:
            vcf = splitline[2]
            vcf = os.path.expanduser(vcf)
        except Exception,e:
            print("There was an error reading the folder.")
            print(e)
        try:
            bam = splitline[3]
            try:
                bam = os.path.expanduser(bam)
            except Exception,e:
                print("There was an error reading the folder.")
                print(e)
        except:
            bam = str()
        
        args.config.write('\t'.join([name, affect, vcf, bam]))
        args.config.write('\n')


else:
    
    while True:
        print("If you supply GATK multi-sample calls, then just mention the same file name, if you're asked for a vcf file.")
        print("Remember. If you want to stop entering data, just leave a field empty and hit enter.")
        name   = raw_input("Please enter the sample ID. The format should be same as the vcf file, you will provide later. : ")
        if name == '':
            print("Field empty. Quitting...")
            break
        
        affect =  raw_input("Please enter if the sample is affected or not. [0 - unaffected, 1 - affected]: ")
        if affect == '':
            print("Field empty. Quitting...")
            break
        
        vcf    = raw_input("Please enter the location of the respective vcf file. : ")
        if vcf == '':
            print("Field empty. Quitting...")
            break
        
        bam    = raw_input("Please enter the location of the respective bam file. : ")
        if bam == '':
            print("Field empty. Quitting...")
            break
        
        args.config.write('\t'.join([name, affect, vcf, bam]))
        args.config.write('\n')