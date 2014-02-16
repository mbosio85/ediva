#! /usr/bin/python

import argparse
import os
import re
import collections
import pprint

parser = argparse.ArgumentParser(description = 'Create a merged vcf file and accordingly read bam files to gather coverage information at the variant positions.')

parser.add_argument('--vcffile',  type=argparse.FileType('r'), dest='vcffile',  required=True, help='VCF file to supplement read data to. [required]')
parser.add_argument('--readfile', type=argparse.FileType('r'), dest='readfile', required=True, help='Mpileup file to read data from. [required]')

parser.add_argument('--outfile',  type=argparse.FileType('w'), dest='outfile',  required=True, help='Output vcf file which will contain read depth data. [required]')

args = parser.parse_args()

# DEBUG
pp = pprint.PrettyPrinter(indent=5)
# /DEBUG

# read header of readfile
sample_columns_pile = dict()
for pileline in args.readfile:
    pileline = pileline.rstrip('\n')
    split_pileline = pileline.split('\t')
    
    # identify columns with pileup information
    if split_pileline[0] == 'chr':
        endindex = len(split_pileline)
        for i in range(3, endindex, 3):
            sample = split_pileline[i].rstrip('count')
            sample_columns_pile[sample] = i
        
        break

# # DEBUG
# pp.pprint(sample_columns_pile)
# exit(0)
# # /DEBUG

# iterate through vcf file
sample_columns_vcf = dict()
for line in args.vcffile:
    
    writetofile = False
    
    # header
    if line.startswith('##'):
        args.outfile.write(line)
        continue
    
    elif line.startswith('#'):
        args.outfile.write(line)
        
        line = line.rstrip('\n')
        splitline = line.split('\t')
        
        endindex = len(splitline)
        
        # identify column positions
        for i in range(9, endindex):
            sample_columns_vcf[splitline[i]] = i
        continue

    line = line.rstrip('\n')
    splitline = line.split('\t')
    
    chromosome  = str(splitline[0])
    position    = int(splitline[1])
    reference   = splitline[3].lower()
    observed    = splitline[4].lower()
    format_info = splitline[8]
    
    # iterate through pile up file to find the respective position of the variant
    for pileline in args.readfile:
        pileline = pileline.rstrip('\n')
        split_pileline = pileline.split('\t')
        
        pile_chromosome = str(split_pileline[0])
        pile_position   = int(split_pileline[1])
        pile_reference  = split_pileline[2].lower()

        
        # position found?
        if pile_chromosome == chromosome and pile_position == position and pile_reference == reference:
            
            new_vcf_field = str()
            
            for sample_vcf in sample_columns_vcf.keys():
                # find the linked information in vcf and pileup file, using samples in vcf file as reference
                sample_pile_count_column = sample_columns_pile[sample_vcf]
                sample_pile_obs_column   = sample_pile_count_column + 1
                
                count_ref = 0
                count_alt = 0
                for letter in split_pileline[sample_pile_obs_column]:
                    if letter == ',' or letter == '.':
                        count_ref += 1
                    elif letter.lower() == observed:
                        count_alt += 1
                
                # now fill up the vcf information
                genotype_info = collections.OrderedDict() #, because the order is important here.
                for fo in format_info.split(':'):
                    # now this dict either has 5 or 6 keys [GT:AD:DP:GQ:PL] [GT:AD:DP:FT:GQ:PL]
                    genotype_info[fo] = None
                
                # split vcf field
                column_vcf = sample_columns_vcf[sample_vcf]
                vcf_field  = splitline[column_vcf]
                
                # supplement the info only if information is missing here
                if vcf_field.startswith('./.'):
                    split_vcf_field = vcf_field.split(':')
                    
                    # get genotype information, sometimes the vcf entries lack values for fields
                    for i in range( len(genotype_info.keys()) ):
                        current_field = genotype_info.keys()[i]
                        try:
                            genotype_info[current_field] = split_vcf_field[i]
                        except IndexError:
                            genotype_info[current_field] = '.'
                    if int(position) == 117803588:
                        
                        pp.pprint( genotype_info.keys() )
                        pp.pprint( genotype_info['FT'] )
                        genotype_info['FT'] = 'FILTER'
                        pp.pprint( genotype_info['FT'] )
                    
                    # check if there is a FT tag and if it was set to PASS, which is not quite logical (since genotype was ./.)
                    try:
                        if genotype_info['FT'] == 'PASS':
                            genotype_info['FT'] = 'FILTER'
                    except:
                        # FT tag does not exist, well so be it...
                        pass
                    
                    # AD field
                    genotype_info['AD'] = "%i,%i" % (count_ref, count_alt)
                    # DP field
                    genotype_info['DP'] = count_ref + count_alt
                    # GQ field - set to 0 since we do not produce any Genotype call
                    genotype_info['GQ'] = 0
                    
                    # stich together & update field
                    new_vcf_field = ':'.join(map( str,genotype_info.values() ))
                    splitline[column_vcf] = new_vcf_field
                    
            # stop iterating through pileup file
            # if quality fits to criteria, write to file. Afterwards continue with next line in vcf
            break
    
    # check if everything has a filter tag
    # no filter tags or 'PASS' is the only thing should be printed to vcf file
    sample_count = 0
    rejected_count = 0
    
    for i in range(9, len(splitline)):
        sample_count += 1
        split_vcf_field = splitline[i].split(':')
        
        # does a 3rd field exist?
        try:
            temp = split_vcf_field[3]
        except:
            # next sample
            continue
        
        # is the 3rd field a number? So no Filter tag exists?
        if split_vcf_field[3].isdigit():
            # fine. next sample
            continue
        # if the filter exist and is 'PASS' things are fine
        elif split_vcf_field[3] == 'PASS':
            # fine. next sample
            continue
        
        # in all other cases reject the vcf field
        rejected_count += 1
        
    # did all fields get rejected?
    if not sample_count == rejected_count:
        # no? OK write to file
        writetofile = True
    
    # write to output file
    if writetofile is True:
        new_vcf_line = '\t'.join(map(str,splitline))
        args.outfile.write(new_vcf_line + '\n')