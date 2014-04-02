import argparse
import csv
import pprint
import sys

parser = argparse.ArgumentParser(description = 'check SNPs, if they fit to the inheritance and make sense. Are the parents really the parents? Currently this tool only works in trios.')
    
parser.add_argument('--infile', type=argparse.FileType('r'), dest='infile', required=True, help='comma separated list of SNPs annotated with mutation impact data. [required]')
parser.add_argument('--family',  type=argparse.FileType('r'), dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')

args = parser.parse_args()

def main (args):
    
    pp = pprint.PrettyPrinter(indent=4)
    
    # read family relationships
    ####
    # expected input:
    # sample	affected
    # VH017	1
    # VH018	0
    # VH019	0
    ####
    
    child = 2
    parent1 = 2
    parent2 = 2
    counter = 0
    
    for line in args.famfile:
        # skip header
        if line.startswith('sample'):
            continue
        
        counter += 1
        if counter > 3:
            sys.exit("More than 3 samples found in family. Unlikely to be trio.")
        
        line = line.rstrip('\n')
        splitline = line.split('\t')
        
        # assign parent and child sample ids to varibles
        if int(splitline[1]) == 1:
            child = splitline[0]
        elif int(splitline[1]) == 0 and parent1 == 2:
            parent1 = splitline[0]
        else:
            parent2 = splitline[0]
    
    pp.pprint([child, parent1, parent2])
    
    # read all data
    alldata = list(csv.reader(args.infile))
    
    header = alldata.pop(0)
    
    index_sample      = identifycolumns(header, 'samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)')
    
    all_fail = 0
    all_ok   = 0
    
    # start reading data
    for line in alldata:
        try:
            samples = line[index_sample]
        except:
            sys.exit("A severe misunderstanding happened, while reading the data. Seemingly some data lines have less fields than the header.")
        
        # input looks like
        # VH011>0/1>2>7>0.777;VH012>./.>1>13>0.928;VH013>1/1>0>20>1
        genotypes = dict()
        
        split_samples = samples.split(';')
        
        for sample in split_samples:
            split_information = sample.split('>')
            sample_id = split_information[0]
            sample_gt = split_information[1]
            
            genotypes[sample_id] = sample_gt
        
        
        if genotypes[child] == '0/0' and ( genotypes[parent1] == '1/1' or genotypes[parent2] == '1/1' ):
            all_fail += 1
        elif genotypes[child] == '0/1' and ( (genotypes[parent1] == '1/1' and genotypes[parent2] == '1/1') or (genotypes[parent1] == '0/0' and genotypes[parent2] == '0/0') ):
            all_fail += 1
        elif genotypes[child] == '1/1' and ( genotypes[parent1] == '0/0' or genotypes[parent2] == '0/0' ):
            all_fail += 1
        else:
            all_ok += 1
    
    failed = all_fail / float(all_fail + all_ok)
    
    pp.pprint([all_fail, all_ok, failed])
    
    if failed > 0.1:
        print "WARNING: The trio might not have the correct relationship. Assuming that less than 10% of all variants violate a logical inheritance pattern."
    else:
        print "REPORT: Checked the relationship in the supplied trio and it looked OK. Assuming that less than 10% of all variants violate a logical inheritance pattern."
    



###############
# SUBROUTINES #
###############


def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        exit("ERROR: %s column could not be identified in the annotated data" % question)
    return( int(index) )


main(args)