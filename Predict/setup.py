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

# if you want to configure in one go on the command line
parser.add_argument('--edivapath', type=str, dest='ediva', required=False, help='The path where eDiVa is located, e.g. ~/ediva/ [optional]')
parser.add_argument('--reference', type=str, dest='reference', required=False, help='The location of your reference genome. (fasta format and indexed with bwa) [optional]')
parser.add_argument('--shoreref', type=str, dest='shore_reference', required=False, help='The location of your reference genome. (fasta format required to be in SHORE format) [optional]')
parser.add_argument('--dbindel', type=str,  dest='dbindel', required=False, help='The location of dbSNP variants, filtered for InDels. [optional]')
parser.add_argument('--dbsnp',   type=str,  dest='dbsnp',   required=False, help='The location of dbSNP variants, filtered for SNPs [optional]')
parser.add_argument('--bwa',     type=str,  dest='bwa',     required=False, help='The location of the BWA program. [optional]')
parser.add_argument('--gatk',    type=str,  dest='gatk',    required=False, help='The location of GATK GenomeAnalysisTK.jar [optional]')
parser.add_argument('--samtools', type=str, dest='samtools', required=False, help='The location of the SAMTOOLS program. [optional]')
parser.add_argument('--novosort', type=str, dest='novosort', required=False, help='The location of the NOVOSORT program. [optional]')
parser.add_argument('--picard',  type=str,  dest='picard',  required=False, help='The path of the PICARD tools. [optional]')
parser.add_argument('--bedtools', type=str, dest='bedtools', required=False, help='The path of BEDTOOLS. [optional]')
parser.add_argument('--clindel', type=str,  dest='clindel', required=False, help='The location of the CLINDEL program. [optional]')
parser.add_argument('--exome',   type=str,  dest='exome',   required=False, help='The location of the EXOME enrichment target files. [optional]')

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
            
        #elif splitline[0] == 'EXOME_SHORE':
        #    exome_shore  = splitline[1].rstrip()

elif args.ediva != None and args.reference != None and args.dbindel != None and args.dbsnp != None and args.bwa != None and args.gatk != None and args.samtools != None and args.novosort != None and args.picard != None and args.bedtools != None and args.exome != None:
    installdir  = os.path.expanduser(args.ediva)
    ref_genome  = os.path.expanduser(args.reference)
    shore_ref   = os.path.expanduser(args.shore_reference)
    dbindel     = os.path.expanduser(args.dbindel)
    dbsnp       = os.path.expanduser(args.dbsnp)
    bwa         = os.path.expanduser(args.bwa)
    gatk        = os.path.expanduser(args.gatk)
    samtools    = os.path.expanduser(args.samtools)
    novosort    = os.path.expanduser(args.novosort)
    clindel     = os.path.expanduser(args.clindel)
    picard      = os.path.expanduser(args.picard)
    bedtools    = os.path.expanduser(args.bedtools)
    exome       = os.path.expanduser(args.exome)

### create a novel config file ###

else:
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

### check for compatibility between reference genome and target bed file
try:
    REF = open(ref_genome, 'r')
except Exception,e:
    print("Could not check for conconrdance between Reference genome and Exome bed file, because:")
    print(e)

try:
    EXO = open(exome, 'r')
except Exception,e:
    print("Could not check for conconrdance between Reference genome and Exome bed file, because:")
    print(e)

# get the chromosome designator of the reference
ref_chrom = str()
for line in REF:
    if line.startswith('>'):
        line = line.lstrip('>')
        splitline = line.split(' ')
        ref_chrom = str(splitline[0])
        break

# get the chromosome designator of the Exome 
counter = 0
exo_chrom = str()
for line in EXO:
    if counter == 2:
        splitline = re.split('\W+', line)
        exo_chrom = str(splitline[0])
    counter += 1

if not len(ref_chrom) == len(exo_chrom):
    print("WARNING: Your Reference file and Exome file seem to have different naming patterns. This will make eDiVa fail later on.")