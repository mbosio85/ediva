import argparse
import pprint



parser = argparse.ArgumentParser(description='filter out variants, that show a . as an alternative allele')

parser.add_argument('--infile',  dest='infile',  type=argparse.FileType('r'), help='vcf input file')
parser.add_argument('--outfile', dest='outfile', type=argparse.FileType('w'), help='vcf output file')

args = parser.parse_args()

for line in args.infile:
    
    # header
    if line.startswith('#'):
        args.outfile.write(line)
        continue
    
    # body
    line = line.rstrip("\n")
    splitline = line.split("\t")
    
    # alt allele only .?
    if splitline[4] == '.':
        continue
    # otherwise write
    else:
        args.outfile.write(line + "\n")