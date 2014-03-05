from subprocess import call
import pprint
import sys
import os
import re
import argparse
pp = pprint.PrettyPrinter(indent = 5)

parser = argparse.ArgumentParser(description = 'Create a file containing all necessary paths for running ediva pipeline.')

parser.add_argument('--oldconfig',  type=argparse.FileType('r'), dest='oldconfig',  required=False, help='If you have already an old file you might want to use it as a template.')
parser.add_argument('--newconfig',  type=argparse.FileType('w'), dest='newconfig',  required=True, help='Location of the ediva config file.')

args = parser.parse_args()

variables = ["EDIVA", "REFERENCE", "SHORE_REFERENCE", "DBINDEL", "DBSNP", "BWA", "GATK", "SAMTOOLS", "NOVOSORT", "PICARD", "BEDTOOLS", "CLINDEL", "EXOME", "EXOME_SHORE"]
installdir  = str()
ref_genome  = str()
shore_ref   = str()
dbindel     = str()
dbsnp       = str()
bwa         = str()
gatk        = str()
samtools    = str()
novosort    = str()
clindel     = str()
picard      = str()
bedtools    = str()
exome       = str()
exome_shore = str()


### read a template config file ###

if args.oldconfig:
    for line in args.oldconfig:
        line = line.rstrip()
        splitline = line.split('=')
        
        if splitline[0] == 'EDIVA':
            installdir = splitline[1].rstrip()
            
        elif splitline[0] == 'REFERENCE':
            ref_genome = splitline[1].rstrip()
            
        elif splitline[0] == 'SHORE_REFERENCE':
            shore_ref  = splitline[1].rstrip()
            
        elif splitline[0] == 'DBINDEL':
            dbindel    = splitline[1].rstrip()
            
        elif splitline[0] == 'DBSNP':
            dbsnp      = splitline[1].rstrip()
            
        elif splitline[0] == 'BWA':
            bwa        = splitline[1].rstrip()
            
        elif splitline[0] == 'GATK':
            gatk       = splitline[1].rstrip()
            
        elif splitline[0] == 'SAMTOOLS':
            samtools     = splitline[1].rstrip()
            
        elif splitline[0] == 'NOVOSORT':
            novosort    = splitline[1].rstrip()
        
        elif splitline[0] == 'CLINDEL':
            clindel     = splitline[1].rstrip()
            
        elif splitline[0] == 'PICARD':
            picard       = splitline[1].rstrip()
            
        elif splitline[0] == 'BEDTOOLS':
            bedtools     = splitline[1].rstrip()
            
        elif splitline[0] == 'EXOME':
            exome        = splitline[1].rstrip()
            
        elif splitline[0] == 'EXOME_SHORE':
            exome_shore  = splitline[1].rstrip()


### create a novel config file ###

for var in variables:
    
    if var == 'EDIVA':
        if installdir == '':
            # guess install folder from called script location
            guessed    = os.path.dirname(os.path.realpath(__file__))
            guessed    = '/'.join(guessed.split('/')[0:-1]) # chomp off last directory level
            installdir = raw_input("Please enter the path eDiVa was installed to. (path only, if you want to use %s just hit enter)" % guessed)
            if installdir == '':
                installdir = guessed
            
        else:
            # ask for intput, but suggest location that has been read from old file
            installdir_manual = raw_input("Please enter the path eDiVa was installed to. (path only, if you want to use %s just hit enter)" % installdir)
            if not installdir_manual == '':
                installdir = installdir_manual
    
    elif var == 'REFERENCE':
        # human genome reference file
        if ref_genome == '':
            ref_genome  = raw_input("Please enter the path of the human genome reference (including file name): ")
        else:
            ref_genome_manual = raw_input("Please enter the path of the human genome reference (including file name, to use %s hit enter): " % ref_genome)
            if not ref_genome_manual == '':
                ref_genome = ref_genome_manual
        ref_genome = re.sub('\'', '', ref_genome)
    
    elif var == 'SHORE_REFERENCE':
        # human genome reference file in SHORE format
        if shore_ref == '':
            shore_ref   = raw_input("Please enter the path of the human genome reference in SHORE format (including file name): ")
        else:
            shore_ref_manual = raw_input("Please enter the path of the human genome reference in SHORE format (including file name, to use %s hit enter): " % shore_ref)
            if not shore_ref_manual == '':
                shore_ref = shore_ref_manual
        shore_ref = re.sub('\'', '', shore_ref)
        
    elif var == 'DBINDEL':
        # dbSNP indels
        if dbindel == '':
            dbindel     = raw_input("Please enter the path of the dbSNP data restricted to InDels (including file name): ")
        else:
            dbindel_manual = raw_input("Please enter the path of the dbSNP data restricted to InDels (including file name, to use %s hit enter): " % dbindel)
            if not dbindel_manual == '':
                dbindel = dbindel_manual
        dbindel = re.sub('\'', '', dbindel)
    
    elif var == 'DBSNP':
        # dbSNP snp
        if dbsnp == '':
            dbsnp     = raw_input("Please enter the path of the dbSNP data restricted to SNPs (including file name): ")
        else:
            dbsnp_manual = raw_input("Please enter the path of the dbSNP data restricted to SNPs (including file name, to use %s hit enter): " % dbsnp)
            if not dbsnp_manual == '':
                dbsnp = dbsnp_manual
        dbsnp = re.sub('\'', '', dbsnp)
    
    elif var == 'BWA':
        if bwa == '':
            bwa         = raw_input("Please enter the path of BWA (including file name): ")
        else:
            bwa_manual  = raw_input("Please enter the path of BWA (including file name, to use %s hit enter): " % bwa)
            if not bwa_manual == '':
                bwa = bwa_manual
        bwa = re.sub('\'', '', bwa)
    
    elif var == 'GATK':
        if gatk == '':
            gatk        = raw_input("Please enter the path of GATK (including file name): ")
        else:
            gatk_manual = raw_input("Please enter the path of GATK (including file name, to use %s hit enter): " % gatk)
            if not gatk_manual == '':
                gatk = gatk_manual
        gatk = re.sub('\'', '', gatk)
    
    elif var == 'SAMTOOLS':
        if samtools == '':
            samtools        = raw_input("Please enter the path of SAMTOOLS (including file name): ")
        else:
            samtools_manual = raw_input("Please enter the path of SAMTOOLS (including file name, to use %s hit enter): " % samtools)
            if not samtools_manual == '':
                samtools = samtools_manual
        samtools = re.sub('\'', '', samtools)
    
    elif var == 'NOVOSORT':
        if novosort == '':
            novosort        = raw_input("Please enter the path of NOVOSORT (including file name): ")
        else:
            novosort_manual = raw_input("Please enter the path of NOVOSORT (including file name, to use %s hit enter): " % novosort)
            if not novosort_manual == '':
                novosort = novosort_manual
        novosort = re.sub('\'', '', novosort)
    
    elif var == 'CLINDEL':
        if clindel == '':
            clindel         = raw_input("Please enter the path of CLINDEL (including file name): ")
        else:
            clindel_manual  = raw_input("Please enter the path of CLINDEL (including file name, to use %s hit enter): " % clindel)
            if not clindel_manual == '':
                clindel = clindel_manual
        clindel = re.sub('\'', '', clindel)
    
    elif var == 'PICARD':
        if picard == '':
            picard        = raw_input("Please enter the path of PICARD (path only): ")
        else:
            picard_manual = raw_input("Please enter the path of PICARD (path only, to use %s hit enter): " % picard)
            if not picard_manual == '':
                picard = picard_manual
        picard = re.sub('\'', '', picard)
    
    elif var == 'BEDTOOLS':
        if bedtools == '':
            bedtools        = raw_input("Please enter the path of BEDTOOLS (path only): ")
        else:
            bedtools_manual = raw_input("Please enter the path of BEDTOOLS (path only, to use %s hit enter): " % bedtools)
            if not bedtools_manual == '':
                bedtools = bedtools_manual
        bedtools = re.sub('\'', '', bedtools)
    
    elif var == 'EXOME':
        if exome == '':
            exome        = raw_input("Please enter the path of EXOME target region file (including file name in BED format): ")
        else:
            exome_manual = raw_input("Please enter the path of EXOME target region file (including file name in BED format, to use %s hit enter): " % exome)
            if not exome_manual == '':
                exome = exome_manual
        exome = re.sub('\'', '', exome)
    
    #elif var == 'EXOME_SHORE':
    #    exome_shore = raw_input("Please enter the path of EXOME bed file (in SHORE format, including file name): ")
    

args.newconfig.write("EDIVA='%s'\n" % installdir)
args.newconfig.write("REFERENCE='%s'\n" % ref_genome)
args.newconfig.write("SHORE_REFERENCE='%s'\n" % shore_ref)
args.newconfig.write("DBINDEL='%s'\n" % dbindel)
args.newconfig.write("DBSNP='%s'\n" % dbsnp)
args.newconfig.write("BWA='%s'\n" % bwa)
args.newconfig.write("GATK='%s'\n" % gatk)
args.newconfig.write("SAMTOOLS='%s'\n" % samtools)
args.newconfig.write("NOVOSORT='%s'\n" % novosort)
args.newconfig.write("CLINDEL='%s'\n" % clindel)
args.newconfig.write("PICARD='%s'\n" % picard)
args.newconfig.write("BEDTOOLS='%s'\n" % bedtools)
args.newconfig.write("EXOME='%s'\n" % exome)
#args.newconfig.write("EXOME_SHORE='%s'\n" % exome_shore)

args.newconfig.close()