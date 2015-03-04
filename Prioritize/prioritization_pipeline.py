import prioritization_support_functions
import pprint
import argparse
import os
import re
if not os.environ.get('DRMAA_LIBRARY_PATH'):
    print "Adding the DRMAA library path to the environment"
    os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
import pickle
import subprocess
import imp
import datetime

curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
curpath = '/'.join(curpath[:-2])
curpath += '/'

qsubclass = imp.load_source('qsubclass', os.path.abspath(curpath + 'pipeline_control/qsubclass.py'))
pipeline_element = imp.load_source('pipeline_element', os.path.abspath(curpath +'pipeline_control/pipeline_element.py'))

python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'
pipe_script  = os.path.abspath(curpath +'/pipeline_control/pipeline_element.py')
qsub_script = os.path.abspath(curpath +'pipeline_control/qsubclass.py')



def multisample_check(multi,vcf):
    predicted = False
    base = vcf[0]
    for i in vcf:
        predicted = i==base            
    
    
    if predicted != multi:
        msg= "The multisample variable is set to %s but the vcf files suggest it should be changed to %s"%(multi,predicted)
        while True:
            print msg
            a =raw_input("Should I change it ? please write 'yes' or 'no'")
            if a == 'yes':
                multi = predicted
                break
            elif a == 'no':
                break
    return multi


########## main routine starts here ##################
pp = pprint.PrettyPrinter(indent = 5)

args = prioritization_support_functions.parse_args()

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
if args.geneexclusion == None:
    print('Warning: Using the default gene exclusion list')
    gene_exclusion_list='$EDIVA/Resource/gene_exclusion_list.txt'
else:
    gene_exclusion_list=args.geneexclusion


# read a family config file, that gives the ID, affection status, vcf location, bam location
variant_string = list()
sample_string  = list()
bam_list    = list()
vcf_list    = list()
sample_list = list()
affection   = dict()
no_bams_given = False
qoptions_def  = False




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
 
# #Automatic check if multisample is set to the correct value   
#args.multisample = multisample_check(args.multisample,vcf_list)

##########################
## NEW PART WITH PIPELINE CONTROL AND AUTORUN ...
##########################
#qclass = qsubclass.disclaimer()
qclass = True

if qclass:
    qlist =list()
    logfile = os.path.abspath(args.outfolder) + '/'+ datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + "qsub_log.log"
    if args.qoptions ==None:
        qoptions,logfile = qsubclass.getOptions()
        qoptions = dict()
        qoptions,logfile = qsubclass.getOptions()
        qoptions['-e'] = list()
        qoptions['-e'].append(str(args.outfolder))
        qoptions['-o'] = list()
        qoptions['-o'].append(str(args.outfolder))
        qoptions['-N'] = [job_name]
        if qoptions.get('-l',False):
            qoptions['-l'].append("virtual_free=20G"%mem)
        else:
            qoptions['-l'] = list()
            qoptions['-l'].append("virtual_free=20G,h_rt=51600")
        qoptions['-pe'] = list() 
    else:
        #print "Automatic Logfile:%s"%logfile
        qoptions = qsubclass.parse_command_options(args.qoptions,args.outfolder,args.qsub_name)
    #print qoptions
    pipe = list()
    
##########################

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

export OUTF=%s
export GATK=%s
export EDIVA=%s
export SAMTOOLS=%s
export REF=%s

DBSNP=%s

""" % (job_name, args.outfolder, args.outfolder, args.outfolder, gatk, ediva, samtools, ref, dbsnp))
env_var = script_content
######
# if a multisample call was given...
######

#  extract information from the multi sample call file and save to combined.variants.supplement.vcf

if args.multisample:
    text= """

# select samples from multisample call file
java -Xmx2g -jar $GATK -R $REF -T SelectVariants --variant %s -o $OUTF/combined.variants.temp.vcf %s -env -ef

# filter out variants, where no sample has more support than 5 reads
%s $EDIVA/Prioritize/post_gatkms_filter.py --infile $OUTF/combined.variants.temp.vcf --outfile $OUTF/combined.variants.vcf

# clean up
rm $OUTF/combined.variants.temp.vcf

""" % (vcf_list[0], sample_joint_string,python_path)
    script_content += text
    p_element = pipeline_element.pipeline_element(env_var+text,"Multisample + filtering")
    p_element.set_error("Error in multisample+filtering executions Please refer to SGE job error file")
    p_element.set_alias("Select multisample variants")
    pipe.append(p_element)


#######
# merges the vcfs (if no multisample call was given)
#######
else:

    # produces a line like this:
    # java -Xmx4g -jar /users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T CombineVariants -R /users/GD/resource/human/hg19/hg19.fasta --variant:40ACVi 40ACVi_indel.vcf --variant:40ACVm 40ACVm_indel.vcf --variant:40ACVp 40ACVp_indel.vcf -o 40ACV/combined.variants.indel.vcf --unsafe LENIENT_VCF_PROCESSING
    text = ("""
# merge vcf files
java -jar $GATK -T CombineVariants -R $REF %s -o $OUTF/combined.variants.vcf --unsafe LENIENT_VCF_PROCESSING
    
    """ % (variant_joint_string))
    script_content += text
    p_element = pipeline_element.pipeline_element(env_var+text,"Non multisample + vcf merging")
    p_element.set_error("Error in merging executions Please refer to SGE job error file")
    filename = args.outfolder+'/combined.variants.vcf'
    p_element.set_condition([filename,None])
    p_element.set_alias("Combine variants")
    pipe.append(p_element)



######
# produces the mpileups (if no multisample call was given)
######

# if not enough (vcf=bam) bam files were given:
# or if there's an uneven number (vcf=bam) of files given:
if ( len(bam_list) == 0 or not len(bam_list) == len(vcf_list) ):
    print """ You did not provide the same amount of bam files and vcf files.
    Just creating a merge vcf script.
    """
    text =  """

# only renaming combined.vcf to supplemented vcf, because supplementing will not be done
cp $OUTF/combined.variants.vcf $OUTF/combined.variants.supplement.vcf

    """
    script_content += text
    p_element = pipeline_element.pipeline_element(env_var+text,"Renaming VCF")
    p_element.set_error("Error in renaming VCF executions Please refer to SGE job error file")
    filename = args.outfolder+'/combined.variants.supplement.vcf'
    p_element.set_condition([filename,None])
    pipe.append(p_element)

# if everything is fine, do genotyping in all family members
elif len(bam_list) == len(vcf_list):
    text= """

# do Genotyping in all family members
java -jar $GATK -T UnifiedGenotyper -R $REF -I %s --dbsnp $DBSNP -o $OUTF/combined.variants.supplement.temp.vcf -alleles $OUTF/combined.variants.vcf --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -glm BOTH

%s $EDIVA/Prioritize/vcf_filter.py --infile $OUTF/combined.variants.supplement.temp.vcf --outfile $OUTF/combined.variants.supplement.vcf

rm $OUTF/combined.variants.supplement.temp.vcf

""" % (list_file,python_path)
    script_content += text
    p_element = pipeline_element.pipeline_element(env_var+text,"Genotyping in all family members")
    p_element.set_error("Error in Genotyping in all family members execution Please refer to SGE job error file")
    filename = args.outfolder+'/combined.variants.supplement.vcf'
    p_element.set_condition([filename,None])
    pipe.append(p_element)

######
# do the annotation
######

text= """


# annotation
#perl $EDIVA/Annotate/annotate.pl --input $OUTF/combined.variants.supplement.vcf --sampleGenotypeMode complete -f
%s $EDIVA/Annotate/annotate.py --input $OUTF/combined.variants.supplement.vcf --sampleGenotypeMode complete -f


"""%(python_path)
script_content += text
p_element = pipeline_element.pipeline_element(env_var+text,"Annotation in all family members")
p_element.set_error("Error in Annotation execution Please refer to SGE job error file")
filename = args.outfolder+'/combined.variants.supplement.sorted.annotated'
p_element.set_condition([filename,None])
pipe.append(p_element)

######
# do the ranking
######

text= """

# rank the variants given
%s $EDIVA/Prioritize/rankSNP.py --infile $OUTF/combined.variants.supplement.sorted.annotated --outfile $OUTF/combined.variants.supplement.ranked

"""%(python_path)
script_content += text
p_element = pipeline_element.pipeline_element(env_var+text,"Ranking in all family members")
p_element.set_error("Error in Ranking execution Please refer to SGE job error file")
filename = args.outfolder+'/combined.variants.supplement.ranked'
p_element.set_condition([filename,None])
pipe.append(p_element)

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
    
    text ="""
    
# run inheritance mode: %s
%s $EDIVA/Prioritize/familySNP.py --infile $OUTF/combined.variants.supplement.ranked  --outfile $OUTF/%s/combined.variants.supplement.%s --filteredoutfile $OUTF/%s/combined.variants.supplement.filtered%s --family $OUTF/pedigree.tree --inheritance %s --familytype %s --geneexclusion %s
    
    """ % ((inhet_mode,python_path, inhet_mode, inhet_mode, inhet_mode, inhet_mode, inhet_mode, args.familytype,gene_exclusion_list))
    script_content += text
    p_element = pipeline_element.pipeline_element(env_var+text,"Inheritance mode %s in all family members"%inhet_mode)
    p_element.set_error("Error in Inheritance %s in all family members execution Please refer to SGE job error file"%inhet_mode)
    filename = args.outfolder+'/%s/combined.variants.supplement.%s'%(inhet_mode,inhet_mode)
    p_element.set_condition([filename,None])
    pipe.append(p_element)
    

# write out the script file
try:
    if os.path.isdir( os.path.dirname(args.qsub_name) ):
        script_location = os.path.expanduser(args.qsub_name)
    else:
        script_location = args.outfolder + '/' + args.job_name
    SCRIPT = open( script_location, 'w')
except Exception,e:
    print "Could not open script file for writing."
    print(e)
    exit(1)

SCRIPT.write(script_content)
SCRIPT.close()

##############
# NEW SECTION FOR PIPELINE CONTROL
###########

pipeline_element.save_pipeline(pipe,script_location+'.pipe')

if qclass:
    command = "%s %s %s %s"%(python_path,pipe_script,script_location+'.pipe','log.log')
    qlist.append(qsubclass.qsubCall(command,qoptions,list(),logfile))
    #print qoptions
    print "\n\ncommand\n"
    command = str.replace(command,'//','/')
    print command
    
    
    
 
if qclass:
    try:
        print "Now I'm saving the queue list in %s file"%(script_location+'.qlist')
        with open(script_location +'.qlist','wb') as outfile:
            pickle.dump(qlist, outfile)
    except :
        print('error saving pipeline list to file')
    print """
    
    Now I'm running the pipeline files:
    
    ############################################################################
    Please remind that the terminal should not be closed for a correct execution
    ############################################################################
        UNLESS:
    If you want to simply let me run in background follow the proposed instructions:
    
        >Ctrl+Z
        >bg
        >disown -h %jid # Where jid is the currend job id. Usually it is 1 but
                        # Please check it before placing a wrong number
    
    
    """    
    qsubclass.run_qlist(script_location+'.qlist')
    try:
        subprocess.call('rm %s',logfile)
    except:
        pass
    print "Finished, have a nice day :)"