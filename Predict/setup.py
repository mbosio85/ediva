from subprocess import call
import pprint
import sys
import os
import re
import argparse
import readline


def parse_argument(variable,msg,sub=True):
    token = 1
    while token:
        if variable == '':
            readline.parse_and_bind("tab: complete")
            readline.set_completer_delims(' \t\n;')
            out_value = raw_input(msg + "\n>")
        else:
            out_value = variable
            # ask for intput, but suggest location that has been read from old file
            readline.parse_and_bind("tab: complete")
            readline.set_completer_delims(' \t\n;')
            out_value_manual = raw_input(msg + "Just hit enter to use %s\n>"%out_value)
            if not out_value_manual == '':
                out_value = out_value_manual
        if len(out_value)>0:
            if out_value[0] =='~':
                out_value = os.path.expanduser("~") + "/"+ out_value[1:] 
            else :
                pass        
        if os.path.exists(out_value):
            token = 0
        else:
            print("please insert a valid path")
    
    if sub:
        ref_genome = re.sub('\'', '', ref_genome)
    return out_value

def main(newconf = False, newfile = ""):
    print("\nWelcome to setup.py\n")
    pp = pprint.PrettyPrinter(indent = 5)
    
    parser = argparse.ArgumentParser(description = 'Create a file containing all necessary paths for running ediva pipeline.')
    
    parser.add_argument('--oldconfig',  type=argparse.FileType('r'), dest='oldconfig',  required=False, help='If you have already an old file you might want to use it as a template.')
    if newconf:
        parser.add_argument('--newconfig',  type=argparse.FileType('w'), dest='newconfig',  required=False, help='Location of the ediva config file.')
    else:
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
    if newconf:
        args.newconfig = open(newfile,'w')
    
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
                installdir = parse_argument(installdir, "Please enter the path eDiVa was installed to. (path only)",False)
            elif var == 'REFERENCE':
                ref_genome = parse_argument(ref_genome, "Please enter the path of the human genome reference (including file name): ")
            elif var == 'SHORE_REFERENCE':
                shore_ref = parse_argument(shore_ref, "Please enter the path of the human genome reference in SHORE format (including file name): ")
            elif var == 'DBINDEL':
                # dbSNP indels
                dbindel = parse_argument(dbindel,"Please enter the path of the dbSNP data restricted to InDels (including file name): ")
            elif var == 'DBSNP':
                # dbSNP snp
                dbsnp = parse_argument(dbsnp,"Please enter the path of the dbSNP data restricted to SNPs (including file name): ")
            elif var == 'BWA':
                bwa = parse_argument(bwa,"Please enter the path of BWA (including file name): ")
            elif var == 'GATK':
                gatk = parse_argument(gatk,"Please enter the path of GATK (including file name): ")
            elif var == 'SAMTOOLS':
                samtools = parse_argument(samtools,"Please enter the path of SAMTOOLS (including file name): ")
            elif var == 'NOVOSORT':
                novosort = parse_argument(novosort,"Please enter the path of NOVOSORT (including file name): ")
            elif var == 'CLINDEL':
                clindel = parse_argument(clindel, "Please enter the path of CLINDEL (including file name): ")
            elif var == 'PICARD':
                picard = parse_argument(picard,"Please enter the path of PICARD (path only): ")        
            elif var == 'BEDTOOLS':
                bedtools = parse_argument(bedtools,"Please enter the path of BEDTOOLS (path only): ")
            elif var == 'EXOME':
                exome = parse_argument(exome,"Please enter the path of EXOME target region file (including file name in BED format): ")
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
        print("Could not check for concordance between Reference genome and Exome bed file, because:")
        print(e)
        exit(0)
    
    try:
        EXO = open(exome, 'r')
    except Exception,e:
        print("Could not check for concordance between Reference genome and Exome bed file, because:")
        print(e)
        exit(0)
    
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
    
    # check if the given files and folders might be correct and contain expected files.
    
    if not os.path.exists(installdir):
        print("WARNING: Your given ediva folder might not be correct.")
    if not os.path.exists(ref_genome):
        print("WARNING: Your given reference genome file might not exist.")
    if not os.path.exists(dbindel):
        print("WARNING: Your given dbSNP Indel file might not exist.")
    if not os.path.exists(dbsnp):
        print("WARNING: Your given dbSNP SNP file might not exist.")
    if not os.path.exists(bwa):
        print("WARNING: Your given bwa file might not exist.")
    if not os.path.exists(gatk):
        print("WARNING: Your given GATK file might not exist.")
    if not os.path.exists(samtools):
        print("WARNING: Your given SAMTOOLS file might not exist.")
    if not os.path.exists(novosort):
        print("WARNING: Your given NovoSort file might not exist.")
    if not os.path.exists(exome):
        print("WARNING: Your given Exome file might not exist.")
    if not os.path.exists( "%s/MarkDuplicates.jar" % (picard) ):
        print("WARNING: Your given Picard folder might not be correct.")
    if not os.path.exists( "%s/intersectBed" % (bedtools) ):
        print("WARNING: Your given BedTools folder might not be correct.")
        
    return None

if __name__ == "__main__":
    main()