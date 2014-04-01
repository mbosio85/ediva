import pprint
import argparse
import os.path
import re

pp = pprint.PrettyPrinter(indent = 5)

parser = argparse.ArgumentParser(description = 'Create a file that will be run in the cluster to fully run the prioritization in eDiVa.')

parser.add_argument('--config',  type=argparse.FileType('r'), dest='config',  required=True, help='The config file you created using setup.py, containing all paths of GATK, SAMTOOLS, etc.')
parser.add_argument('--family',  type=argparse.FileType('r'), dest='family',  required=True, help='A file containing all paths to vcf and bam files, including the affection status. Can be created with collect_family_info.py')
parser.add_argument('--outfolder', type=str, dest='outfolder', required=True, help='The folder where the output data should be written to.')
parser.add_argument('--inheritance', choices=['dominant_denovo', 'dominant_inherited', 'recessive', 'Xlinked', 'compound'], dest='inheritance', required=True, action = 'append', help="""choose a inheritance model [required]
This option can be called multiple times.
dominant_inherited: used for families
dominant_denovo: apply to novel variants seen in the affected individuals

recessive: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
Xlinked: used for X linked recessive variants in trios only
compound: detect compound heterozygous recessive variants
""")
parser.add_argument('--familytype', choices=['trio', 'family'], dest='familytype', required=True, help="choose if the data you provide is a trio or a larger family")
parser.add_argument('--multisample', required=False, action='store_true', help="if your input variants have been called using GATK multi sample calling, you should toggle this option.")
parser.add_argument('--qsubname', type=str, dest='qsub_name', required=True, help='name of the script you want to start later on')
parser.add_argument('--jobname', type=str, dest='job_name', required=False, default='prioritize', help='name of the job that will be run on your cluster [default: prioritize]')
parser.add_argument('--force', action = 'store_true', help='Enable eDiVa to overwrite old output')

args = parser.parse_args()

# does the output folder exist?
if not os.path.exists(args.outfolder):
    print "Creating output folder %s" % args.outfolder
    os.makedirs(args.outfolder)

# check if pipeline was run in that folder already
if os.path.isfile("%s/combined.variants.vcf" % args.outfolder) and not args.force:
    print("The output files of prior run still exist in the specified folder. Please use the --force option or choose another folder.")

# read the general config file for paths of GATK, eDiVa and so on

for line in args.config:
    line = line.rstrip('\n')
    splitline = line.split('=')
    
    if splitline[0] == 'EDIVA':
        ediva = splitline[1]
        
    elif splitline[0] == 'REFERENCE':
        ref = splitline[1]
        
    elif splitline[0] == 'SHORE_REFERENCE':
        shore_ref = splitline[1]
    
    elif splitline[0] == 'DBINDEL':
        dbindel = splitline[1]
        
    elif splitline[0] == 'DBSNP':
        dbsnp = splitline[1];

    elif splitline[0] == 'BWA':
        bwa = splitline[1];

    elif splitline[0] == 'GATK':
        gatk = splitline[1];

    elif splitline[0] == 'SAMTOOLS':
        samtools = splitline[1];

    elif splitline[0] == 'NOVOSORT':
        novosort = splitline[1];

    elif splitline[0] == 'PICARD':
        picard = splitline[1];

    elif splitline[0] == 'BEDTOOLS':
        bedtools = splitline[1];

    elif splitline[0] == 'CLINDEL':
        clindel = splitline[1];

    elif splitline[0] == 'EXOME':
        exome = splitline[1];


# read a family config file, that gives the ID, affection status, vcf location, bam location
variant_string = list()
sample_string  = list()
bam_list    = list()
vcf_list    = list()
sample_list = list()
affection   = dict()

no_bams_given = False

for line in args.family:
    
    # find the header and skip it
    m = re.search("ID.*status.*vcf.*bam", line)
    
    if m:
        continue
    
    
    line = line.rstrip('\n')
    splitline = line.split('\t')
    
    sample = splitline[0]
    affect = splitline[1]
    vcf    = splitline[2]
    try: # only try to append the bam name, because this is optional
        bam    = splitline[3]
        bam_list.append(bam)
    except:
        pass
    
    # save all sample IDs
    sample_list.append(sample)
    # collect samples in a string, that can be used to select samples from multi sample calls
    sample_string.append( "-sn %s" % (sample) )
    # save all vcf names
    vcf_list.append(vcf)
    # collect elements for a string, that can be used to merge the vcf files
    variant_string.append("--variant:%s %s" % (sample, vcf))
    # create a dictionary linking sample and disease status
    affection[sample] = affect
    
# finally - create a string, that can be used to merge the vcf files
variant_joint_string = ' '.join(variant_string)

# create a string, that can be used to select samples out of a multi sample call file
sample_joint_string  = ' '.join(sample_string)


# create header
id_elements = list()
for element in sample_list:
    id_elements.append(element + "count")
    id_elements.append(element + "obs")
    id_elements.append(element + "qual")
header = '\\t'.join(id_elements)

# write a pedigree file
try:
    pedigree_file = args.outfolder + '/pedigree.tree'
    FH = open(pedigree_file, 'w')
    
    # write header
    outline = '\t'.join(["sample", "affected"])
    FH.write(outline + "\n")
    
    for sample in affection.keys():
        affect = affection[sample]
        outline = '\t'.join([sample, affect])
        FH.write(outline + "\n")
    
    FH.close()
except Exception,e:
    print("There was something wrong writing the pedigree file, better aborting.")
    print(str(e))
    exit(1)

# write out a list of bam files
try:
    bam_file_list = '\n'.join(bam_list)
    
    if len(bam_list) < len(vcf_list):
        no_bams_given = True
    
    if no_bams_given is False:
        # save the locations of the bam files to a temporary file
        try:
            list_file = args.outfolder + '/bam.list'
            FH = open(list_file, 'w')
        except Exception,e:
            print "Could not open output folder for writing."
            print(str(e))
            exit(1)
    
        FH.write(bam_file_list)
        FH.close()
except:
    print "No bam files given. Only producing a merge vcf script." # not giving a bam file is not yet implemented


# qsub doesn't like jobs starting with a digit
if args.job_name[0].isdigit():
    job_name = 'p' + args.job_name
else:
    job_name = args.job_name

# create a qsub script, that:
script_content = str()

# bash script header
script_content = ("""
#!/bin/bash

#$ -N %s
#$ -e %s
#$ -o %s

source /etc/profile
export _JAVA_OPTIONS=\"-Djava.io.tmpdir=$TMPDIR $_JAVA_OPTIONS\"

OUTF=%s
GATK=%s
EDIVA=%s
SAMTOOLS=%s
REF=%s

""" % (job_name, args.outfolder, args.outfolder, args.outfolder, gatk, ediva, samtools, ref))

######
# if a multisample call was given...
######

#  extract information from the multi sample call file and save to combined.variants.supplement.vcf

if args.multisample is True:
    script_content = script_content + """

# select samples from multisample call file
java -Xmx2g -jar $GATK -R $REF -T SelectVariants --variant %s -o $OUTF/combined.variants.temp.vcf %s -env -ef

# filter out variants, where no sample has more support than 5 reads
perl $EDIVA/Prioritize/post_gatkms_filter.pl --infile $OUTF/combined.variants.temp.vcf --outfile $OUTF/combined.variants.vcf

# clean up
rm $OUTF/combined.variants.temp.vcf

""" % (vcf_list[0], sample_joint_string)


#######
# merges the vcfs (if no multisample call was given)
#######

if args.multisample is False:

    # produces a line like this:
    # java -Xmx4g -jar /users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T CombineVariants -R /users/GD/resource/human/hg19/hg19.fasta --variant:40ACVi 40ACVi_indel.vcf --variant:40ACVm 40ACVm_indel.vcf --variant:40ACVp 40ACVp_indel.vcf -o 40ACV/combined.variants.indel.vcf --unsafe LENIENT_VCF_PROCESSING
    script_content = script_content + ("""
# merge vcf files
java -jar $GATK -T CombineVariants -R $REF %s -o $OUTF/combined.variants.vcf --unsafe LENIENT_VCF_PROCESSING
    
    """ % (variant_joint_string))




######
# produces the mpileups (if no multisample call was given)
######

# if not enough (vcf=bam) bam files were given:
# or if there's an uneven number (vcf=bam) of files given:
if ( len(bam_list) == 0 or not len(bam_list) == len(vcf_list) ):
    print """ You did not provide the same amount of bam files and vcf files.
    Just creating a merge vcf script.
    """
    script_content = script_content + """

# only renaming combined.vcf to supplemented vcf, because supplementing will not be done
cp $OUTF/combined.variants.vcf $OUTF/combined.variants.supplement.vcf

    """
    

# if everything is fine, add the mpileup part of the script
elif len(bam_list) == len(vcf_list):
    script_content = script_content + """

# pileup call
# create position list
sed -e 's/chr//' $OUTF/combined.variants.vcf | awk '{OFS=\"\\t\"; if (!/^#/){print $1,$2}}'  > $OUTF/variant.position.bed

# header
echo -e \"chr\\tstart\\tref\\t%s\" > $OUTF/variant.position.mpileup

# pileup
$SAMTOOLS mpileup -l $OUTF/variant.position.bed -f $REF -b $OUTF/bam.list >> $OUTF/variant.position.mpileup


# do the supplementing
python $EDIVA/Predict/supplement_vcf.py --vcffile $OUTF/combined.variants.vcf --readfile $OUTF/variant.position.mpileup --outfile $OUTF/combined.variants.supplement.vcf

""" % (header)


######
# do the annotation
######

script_content = script_content + """


# annotation
perl $EDIVA/Annotate/annotate.pl --input $OUTF/combined.variants.supplement.vcf --sampleGenotypeMode complete -f

"""
######
# do the ranking
######

script_content = script_content + """

# rank the variants given
python $EDIVA/Prioritize/rankSNP.py --infile $OUTF/combined.variants.supplement.sorted.annotated --outfile $OUTF/combined.variants.supplement.ranked

"""

######
# inheritance pattern
######

# take the inheritance pattern(s) to be interogated
# create a subfolder for each inheritance mode wished and
# runs the inheritance script for the pattern

for inhet_mode in args.inheritance:
    
    # create a subfolder for each inheritance mode wished and
    inhet_folder = args.outfolder + "/" + inhet_mode
    if not os.path.exists(inhet_folder):
        print "Creating output folder %s" % inhet_folder
        os.makedirs(inhet_folder)
    
    script_content = script_content + """
    
# run inheritance mode: %s
python $EDIVA/Prioritize/familySNP.py --infile $OUTF/combined.variants.supplement.ranked  --outfile $OUTF/%s/combined.variants.supplement.%s --filteredoutfile $OUTF/%s/combined.variants.supplement.filtered%s --family $OUTF/pedigree.tree --inheritance %s --familytype %s --geneexclusion $EDIVA/Resource/gene_exclusion_list.txt
    
    """ % ((inhet_mode, inhet_mode, inhet_mode, inhet_mode, inhet_mode, inhet_mode, args.familytype))
    pass

# write out the script file
try:
    if os.path.isdir( os.path.dirname(args.qsub_name) ):
        script_location = os.path.expanduser(args.qsub_name)
    else:
        script_location = args.outfolder + args.qsub_name
    SCRIPT = open( script_location, 'w')
except Exception,e:
    print "Could not open script file for writing."
    print(e)
    exit(1)

SCRIPT.write(script_content)
SCRIPT.close()