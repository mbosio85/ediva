#!/usr/bin/python -u
# Filename: prediction_support_functions.py
import readline
import os
import argparse
import subprocess
curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
predictpath = '/'.join(curpath[:-2])
predictpath+='/Predict/'
curpath = '/'.join(curpath[:-1])
curpath += '/'

python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'


cursor ="\n>"
def usage() :
    print ("\n usage:\n"
           "--infolder \t folder containing the raw sequence data in fastq format \n"
            "--outfolder \t folder to write the analysis to \n"
            "--qsubname \t name of the script you want to start lateron \n"
            "--config \t name of the config file to be used \n"
            "--max_coverage \t used in SNP filtering with samtools [default = 400] \n"
            "--namestart \t  start of a substring in the read file name (first letter is numbered 1) \n"
            "--namelength \t length of the part of the filenames which should be taken as sample (and folder) name \n"
            "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n"
            "--secondreadextension \t describes the name of the first read file (e.g. 2.fastq.gz or 2_sequence.txt.gz)\n"
            "--indelcaller \t sets the indel caller to be used. Choose between gatk or clindel [default = clindel]\n"
            "--fusevariants \t call this parameter to write all variants (SNPs and InDels) into one file [sample_name.vcf]\n"
            "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 4]\n"
            "--mem \t amount of Memory dedicated to your job in Gb. The amount of memory must not be bigger than available in the machine. [default = 12]\n"
            "--sample_list \t comma separated list of sample filenames, it can be used to select just few samples in the folder"
            "--help \t\t show help \n")
    return None   
    
    
def parse_commandline():
    parser_ = dict()
    commandLine = True
    try:
        parser = argparse.ArgumentParser(description = '.',usage='below is explained')
        parser.add_argument('--infolder'  ,             type=str, dest ="infolder",     required=True, 			    	 help= "--infolder \t folder containing the raw sequence data in fastq format \n")
        parser.add_argument('--outfolder',              type=str, dest ="outfolder",    required=True,				 help= "--outfolder \t folder to write the analysis to \n")
        parser.add_argument('--qsubname'  ,             type=str, dest ="qsubname",     required=True, 				 help= "--qsubname \t name of the script you want to start lateron \n")
        parser.add_argument('--config',		               dest ="config",       required=True,                              help= "--config \t name of the config file to be used \n")
        parser.add_argument('--max_coverage', 	    type =int,dest ="max_coverage",     required=False,                              help= "--max_coverage \t used in SNP filtering with samtools [default = 400] \n")
        parser.add_argument('--namestart'  ,            type=int, dest ="namestart",     required=True, 			 help= "--namestart \t  start of a substring in the read file name (first letter is numbered 1) \n")
        parser.add_argument('--namelength',             type=int, dest ="namelength",    required=True,         	         help= "--namelength \t length of the part of the filenames which should be taken as sample (and folder) name \n")
        parser.add_argument('--firstreadextension',     type=str, dest ="firstreadextension",      required=True,	         help= "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n")
        parser.add_argument('--secondreadextension',    type=str, dest ="secondreadextension",      required=True, 		 help= "--secondreadextension \t describes the name of the second read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n")  
        parser.add_argument('--indelcaller',            type=str, dest ="indelcaller",  required=False,      		        help= "--indelcaller \t sets the indel caller to be used. Choose between gatk or clindel [default = clindel]\n")  
        parser.add_argument('--fusevariants',                     dest ="fusevariants", action='store_true' , required=False,	 help= "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n")  
        parser.add_argument('--cpu',                    type=int, dest ="cpu",          required=False,      		        help= "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 4]\n")  
        parser.add_argument('--mem',                    type=int, dest ="mem",          required=False,     		        help= "--mem \t amount of Memory dedicated to your job in Gb. The amount of memory must not be bigger than available in the machine. [default = 12]\n")
        parser.add_argument('--qoptions',               type=str, dest ="qoptions",     required=False,     		        help= "--qoptions \t qsub options with this format: -X val,-Y val,-Z, -XX value\n")
        parser.add_argument('--sample_list',            type=str, dest ="sample_list",  required=False,     		        help= "--sample_list \t comma separated list of sample filenames, it can be used to select just few samples in the folder")
        Err_message = "not enough input arguments for command line launch: "    
        args = parser.parse_args()
        if args.infolder:
            if os.path.isdir(args.infolder):
                parser_["infolder"] = args.infolder
            else:
                Err_message  = "ERROR: --infolder is not a recognized directory."
                raise SystemExit
        if args.outfolder:
            if os.path.isdir(args.outfolder):
                parser_["outfolder"]=args.outfolder
            else:
                try:
                    os.makedirs(args.outfolder)
                except OSError as exeption:
                    if exception.errno != os.errno.EEXIST:
                        raise SystemExit           
        if args.qsubname:
            parser_["qsubname"] = args.qsubname
        if args.config:
            if os.path.isfile(args.config):
                parser_["config"] =args.config
            else:
                Err_message  = "ERROR: --config is not a recognized file."
                raise SystemExit
        if args.max_coverage:
            parser_["max_coverage"]=args.max_coverage
        if args.namestart :
            parser_["namestart"] =  args.namestart        
        if args.namelength:
            parser_["namelength"] = args.namelength
        if args.firstreadextension:
            parser_["firstreadextension"] = args.firstreadextension
        if args.secondreadextension:
            parser_["secondreadextension"] = args.secondreadextension
        if args.indelcaller:
            possibilities = ["gatk", "clindel", "both"]
            if args.indelcaller in possibilities:
                parser_["indelcaller"] = args.indelcaller
            else:
                Err_message += "ERROR: the indel caller should be gatk , clined, both"
                raise SystemExit
        if args.cpu:
            parser_["cpu"] = args.cpu
        if args.mem:
            parser_["mem"] = args.mem
        if args.fusevariants:
            parser_["fusevariants"] = 'yes'
        else:
            parser_["fusevariants"] = 'no'
        if args.qoptions:
            parser_["qoptions"] = args.qoptions
        else:
            Err_message+= "ERROR: somthing wrong with qoptions"
        if args.sample_list:
            parser_["sample_list"] = args.sample_list
            
            
        ## now that we have the variables - check for the required ones if they exist or not:
        
        
    except SystemExit :
        commandLine = False
        usage()
        print Err_message+"\nResetting and starting with the interactive mode"
        

      
    
    
    return (commandLine,parser_)


def process_config_file(text):
    installdir =""
    ref_genome=""
    shore_ref=""
    dbindel=""
    dbsnp=""
    bwa =""
    gatk =""
    samtools=""
    novosort=""
    clindel=""
    picard=""
    bedtools=""
    exome=""
    ##now actually parse them from config file
    with open(text) as oldconfig:  
        for line in oldconfig:
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
    return{"installdir":installdir,"ref_genome":ref_genome,"shore_ref":shore_ref,"dbindel":dbindel,"dbsnp":dbsnp, "bwa":bwa,
           "gatk":gatk,"samtools":samtools,"novosort":novosort,"clindel":clindel,"picard":picard,"bedtools":bedtools,"exome":exome}


def parse_config_file():
    token =1
    while token:
        answer = raw_input("Do you already have a config file ? please print 'y' or 'n'\n")
        if answer == 'y':
            token = 0
        elif answer =='n':
            while token:
                readline.parse_and_bind("tab: complete")
                readline.set_completer_delims(' \t\n;')
                text = raw_input("please insert the full path to your new config file\n >")
                if len(text)>0:
                    if text[0] =='~':
                        text = os.path.expanduser("~") + "/"+ text[1:] 
                    else :
                        pass                                            
                #if (os.path.isfile(text)):
                #cmd =  ("python ./setup.py --newconfig %s"%text)
                cmd =  ("%s %s/setup.py --newconfig %s"%(python_path,predictpath,text))
                #print cmd
                p=subprocess.call(cmd,shell=True)
                stdout,stderr = p.communicate()
                p.wait()
                token = 0
                
    
    
    token =1
    prompt_0 = "Please insert the config file \n >"
    prompt = prompt_0
    while token : 
        #outputs a prompt, parse it to a full path and  checks that it works
        readline.parse_and_bind("tab: complete")
        readline.set_completer_delims(' \t\n;')
        text  = raw_input(prompt)
        if len(text)>0:
            if text[0] =='~':
                text = os.path.expanduser("~") + "/"+ text[1:] 
            else :
                pass        
                            
            if (os.path.isfile(text)):
                token =0
            else:
                prompt = "Error: the chosen path does not exist. Please correct it" + prompt_0
        else :
            prompt = "Please insert something: \n " + prompt_0
            
    print (">The chosen configuration file is %s" % (text))
    #Here we load the file and check for all the parameters:
    #init default variables:
    return process_config_file(text)

#---------------------------------------------------
#This one gets the path from raw_input and veryfies it exists
#used for infolder and outfolder [that case it checks if the folder exists and if not, it creates it]
def get_path(prompt,outfolder):
    prompt_0 = prompt
    token =1
    while token : 
        #outputs a prompt, parse it to a full path and  checks that it works
        readline.parse_and_bind("tab: complete")
        readline.set_completer_delims(' \t\n;')
        text  = raw_input(prompt)
        if len(text)>0:
            if text[0] =='~':
                text = os.path.expanduser("~") + "/"+ text[1:] 
            else :
                pass        
                
            if text[-1] != "/":
                text = text + "/"                
            else:
                pass
            
            if (os.path.exists(text)):
                token =0
            else:
                if not(outfolder):
                    prompt = "Error: the chosen path does not exist. Please correct it" + prompt_0
                else:
                    try:
                        os.makedirs(text)
                    except OSError as exeption:
                        if exception.errno != os.errno.EEXIST:
                            raise
        else :
            prompt = "Please insert something: \n " + prompt_0
            
    #missing verifying that the directory exists and if it's the case to create it
    print (">The chosen path is %s" % (text))
    return text
#---------------------------------------------------
    
def parse_input(prompt,min_v, max_v,def_v):
    token = 1
    s = None
    while token:
        tmp = raw_input(prompt)
        if len(tmp)<1:
            token = 0 
            return def_v
        else:
            try:
                s = int(tmp)
                if ((s<min_v) or (s>max_v)):
                    raise Exception                    
                else:
                    token = 0
                    
            except :
                tmp =  "Wrong input please put a number between %d and %d %s" %(min_v,max_v,prompt)
                prompt = tmp

    return s

## MANAGES THE INDEL CHOICE

def complete_indel(text,state):
    possibilities = ["gatk", "clindel", "both"]
    results = [x for x in possibilities if x.startswith(text)]
    return results[state]

def indel_parser(def_value):
    possibilities = ["gatk", "clindel", "both"]
    readline.parse_and_bind("tab: complete")
    readline.set_completer(complete_indel)
    token = 1
    out_str = None
    while token:
        out_str = raw_input("please choose the indel caller. \n Hit return for default = %s  %s"%(def_value,cursor))
        if (out_str in possibilities):
            token = 0
        elif len(out_str)==0:
            token = 0
            out_str = def_value
    print (">Indel parser: %s"%(out_str))
    return out_str

## MANAGES THE FUSEPARSER CHOICE

def complete_fusevariants(text,state):
    possibilities = ["yes", "no"]
    results = [x for x in possibilities if x.startswith(text)]
    return results[state]

def fuse_parser(fuse_def):
    print "Want to fuse the variants ? Default yes. Type no to change it. Otherwise hit return"
    token = 1
    out_str = None
    dictionary =["yes","no"]
    readline.parse_and_bind("tab: complete")
    readline.set_completer(complete_fusevariants)
    while token:
        out_str = raw_input("please choose if fusevariants is 'yes' or 'no' Default is %s %s"%(fuse_def, cursor))
        if len(out_str)<1:
            out_str = fuse_def
            token =0
        elif (out_str in dictionary):
            token = 0
    print (">Fusevariants: %s"%(out_str))
    return out_str

#Predicts and gathers the start_name and namelength form infolder
def get_namelength(infolder):
    namestart = 1
    namelength = 6
    filenames = next(os.walk(infolder))[2]
    tmp = list()
    for i in range(0,len(filenames)):
        tmp_var = filenames[i].split('.')
        tmp.append(tmp_var[0])

    namelength = len(tmp[0])
    print("Inferred nameStart: %d and NameLength: %d  "%(namestart,namelength))
    print("The result is %s"%(tmp[0]))
    token = 1
    token2= 1
    while token:
        print("If you agree, just hit return")
        print("Otherwise, write no")
        tmp_text = raw_input()
        if len(tmp_text)==0:
            token = 0
            token2 =0
        elif tmp_text =="no":
            token =0
        else:
            pass
    if token2:
        namestart = parse_input("Please insert the namestart. Current value= %d %s"%(namestart,cursor),1,1000000,namestart)
        namelength = parse_input("Please insert the namelength. Current value= %d %s"%(namelength,cursor),1,1000000,namelength);
    tmp = tmp[0]
    print("The result is: \n"
          "Namestart = %d \n Namelength : %d \n Result: %s"
          %(namestart,namelength,tmp[namestart-1:namelength]))
    return(namestart,namelength)


#Predicts and gathers name extension from the infolder
def get_extensions(infolder):
    tmp = list()
    filenames = next(os.walk(infolder))[2]
    for i in range(0,len(filenames)):
        extension = ""
        tmp_var = filenames[i].split('.')
        for j in range(1,len(tmp_var)):
            extension = extension + '.' + tmp_var[j]
        tmp.append(extension)
    first_e = tmp[0]
    second_e = tmp[1]
    print("Inferred firstExtension: %s \n and secondExtension: %s  "%(first_e,second_e))
    token = 1
    token2=1
    while token:
        print("If you agree, just hit return")
        print("Otherwise, write no")
        tmp_text = raw_input()
        if len(tmp_text)==0:
            token = 0
            token2 = 0
        elif tmp_text =="no":
            token =0
        else:
            pass
    if token2:
        first_e  = raw_input("Please insert the first read extension \n >")
        second_e = raw_input("Please insert the second read extension \n >")
    print("The result is \n"
          "First read extension : %s \n Second read extension : %s"
          %(first_e,second_e))
    return (first_e,second_e)


#######################################################################
########################################################################
# PIPELINE SCRIPT WRITING
#######################################################################
#######################################################################
def find_samples(infolder, namestart,namelength,first_e, second_e,sample_list):
    samples =  list()
    if sample_list == None:
        first_read  = [each for each in os.listdir(infolder) if each.endswith(first_e)]
        second_read = [each for each in os.listdir(infolder) if each.endswith(second_e)]
        for i in range(0,len(first_read)):
            tmp = first_read[i]
            tmp = tmp[:-len(first_e)]
            first_read[i] = tmp
        for i in range(0,len(second_read)):
            tmp = second_read[i]
            tmp = tmp[:-len(second_e)]
            second_read[i] = tmp
        samples = list(set(first_read) & set(second_read))
    else:
        #samples are already selected:
        samples = sample_list.split(',')
        
    return(samples)

def make_dirs(outfolder,sample):
    print outfolder +'/'+ sample
    try:
        os.makedirs(outfolder+'/'+sample)
    except OSError as exx:
        if exx.errno != os.errno.EEXIST:
            raise
    try:
        os.makedirs(outfolder+'/'+sample+'/intermediate_files/')
    except OSError as exx:
        if exx.errno != os.errno.EEXIST:
            raise
    return(0)
############################
#in_paths fields [will be used in the header]
#Parsed_config fields:
#installdir
#ref_genome
#shore_ref
#dbindel
#dbsnp
#bwa
#gatk
#samtools
#novosort
#clindel
#picard
#bedtools
#exome

##########################
# in_vars fields [will be used in all following methods to write the pipeline to the script file]
# sname
# qsubname
# first_e
# second_e
# outfolder
# cpu
# mem
# indel


def write_header(script,in_paths,in_vars):
    '''    # Header and references ^^^^^   '''
    value_dictionary = dict(in_paths.items() + in_vars.items())

    header = ("""
            #!/bin/bash
            set -e
            #$ -N {0[sdir]}_{0[qsubname]}
            #$ -e {0[outfolder]}/{0[sdir]}
            #$ -o {0[outfolder]}/{0[sdir]}
            #$ -pe smp {0[cpu]}
            #$ -l virtual_free={0[mem]}G
            
            ### source /etc/profile
            """.format(value_dictionary))
    env_var =("""
        export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMPDIR $_JAVA_OPTIONS"
       
        NAME={0[sdir]}
        READ1={0[infolder]}/{0[sname]}{0[first_e]}
        READ2={0[infolder]}/{0[sname]}{0[second_e]}
        OUTF={0[outfolder]}/{0[sdir]}
       
       
       
        REF={0[ref_genome]}
        SHOREREF={0[shore_ref]}
        DBINDEL={0[dbindel]}
        DBSNP={0[dbsnp]}
        BWA={0[bwa]}
        EDIVA={0[installdir]}
        GATK={0[gatk]}
        PICARD={0[picard]}
        SAMTOOLS={0[samtools]}
        NOVOSORT={0[novosort]}
        BEDTOOLS={0[bedtools]}
        CLINDEL={0[clindel]}
        EXOME={0[exome]}
        #TMPDIR={0[outfolder]}/{0[sdir]}       
    """.format(value_dictionary))

    script.write(header)
    script.write(env_var)

    return header,env_var

def seq_pipeline(script,in_paths,in_vars):
    '''
    #   BWA alignment    ^^^^^^^^^^^^
    #   Quality check and 64-33 transform    
    #   Sort BAM file + cleanup    
    #   Local Re alignment    
    #   Duplicate marking    
    #   Base quality recalibration    
    #   Cleanup    
    #   GATK SNP and INDEL with GATK Haplotype Caller    
    #   Filter and compare SNP    
    #   Correct Sample names in VCF file    
    #   Intersecting    
    #   Evaluation    
    #   Annotate Enrichment    
    #   Borrow Header from GATK vcf File '''
    textlist = list()
    condition_list = list()
    textlist.append("""### Align reads with bwa
            $BWA mem -M -t {0[cpu]}  -R "@RG\\tID:$NAME\\tSM:$NAME" $REF $READ1 $READ2 | time $SAMTOOLS view -h -b -S -F 0x900 -  > $TMPDIR/$NAME.noChimeric.bam
            
            ### check for Quality encoding and transform to 33 if 64 encoding is encountered
            OFFSET=$($SAMTOOLS view $TMPDIR/$NAME.noChimeric.bam | python $EDIVA/Predict/whichQuality_bam.py)
            if [[ $OFFSET == 64 ]];
            then
                echo "fixing 64 quality encoding"
                $SAMTOOLS view -h $TMPDIR/$NAME.noChimeric.bam | python $EDIVA/Predict/bam_rescale_quals.py - | $SAMTOOLS view -bS - > $TMPDIR/$NAME.transformed.bam
                rm $TMPDIR/$NAME.noChimeric.bam
                mv $TMPDIR/$NAME.transformed.bam $TMPDIR/$NAME.noChimeric.bam
            fi
            
            cp $TMPDIR/$NAME.noChimeric.bam $OUTF/intermediate_files/$NAME.noChimeric.bam ##
    """.format(in_vars))
    condition_list.append(None)
    
    textlist.append("""### Sort BAM file
            if ! [ -s $TMPDIR/$NAME.noChimeric.bam ] && [ -s $OUTF/intermediate_files/$NAME.noChimeric.bam ];
            then
                #copy from intermediate folder
                cp $OUTF/intermediate_files/$NAME.noChimeric.bam $TMPDIR/$NAME.noChimeric.bam
            fi
                    
            if [ -s $TMPDIR/$NAME.noChimeric.bam ];
            then
                echo Sort BAM
                $NOVOSORT --threads {0[cpu]} --tmpdir $TMPDIR --forcesort --output $TMPDIR/$NAME.sort.bam -i -m {0[mem]}G $TMPDIR/$NAME.noChimeric.bam
                cp $TMPDIR/$NAME.sort.bam* $OUTF/
                # clean up
                rm $TMPDIR/$NAME.noChimeric.bam
            else
                echo $TMPDIR/$NAME.noChimeric.bam not found >&2
                exit 1
            fi
    """.format(in_vars))
    #condition_list.append([in_vars["outfolder"]+'/'+in_vars["sdir"]+".sort.bam" ,None, in_vars["outfolder"]+'/'+in_vars["sdir"]+".sort.bam.bai" ,None])
    condition_list.append(None)
    
    textlist.append("""### Local Re-alignment
            if ! [ -s  $TMPDIR/$NAME.sort.bam.bai ] && [ -s $OUTF/intermediate_files/$NAME.sort.bam.bai ];
                then
                    cp   $OUTF/intermediate_files/$NAME.sort.bam.bai $TMPDIR/$NAME.sort.bam.bai
            fi
            
            if ! [ -s  $TMPDIR/$NAME.sort.bam ] && [ -s $OUTF/intermediate_files/$NAME.sort.bam ];
                then
                    cp   $OUTF/intermediate_files/$NAME.sort.bam $TMPDIR/$NAME.sort.bam
            fi
            
            
            if [ -s $TMPDIR/$NAME.sort.bam.bai ] && [ -s $TMPDIR/$NAME.sort.bam ] ;
            then
               echo Local Re-alignment
               echo -e "\\n #### doing Local Realignment: java -jar $GATK -nt {0[cpu]} -T RealignerTargetCreator -R $REF -I $TMPDIR/$NAME.sort.bam -o $OUTF/$NAME.intervals -known $DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE \\n\"
               java -jar $GATK -nt {0[cpu]} -T RealignerTargetCreator -R $REF -I $TMPDIR/$NAME.sort.bam -o $OUTF/$NAME.intervals -known $DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE # -L $EXOME
               java -jar $GATK  -T IndelRealigner -R $REF -I $TMPDIR/$NAME.sort.bam -targetIntervals $OUTF/$NAME.intervals -o $TMPDIR/$NAME.realigned.bam -known $DBINDEL --maxReadsForRealignment 10000 --consensusDeterminationModel USE_SW --downsampling_type NONE # -L $EXOME
            else
               echo $NAME.sort.bam.bai not found >&2
               exit 1
            fi
            
            cp $TMPDIR/$NAME.realigned.bam*  $OUTF/intermediate_files/
    """.format(in_vars))
    #condition_list.append([in_vars["outfolder"]+'/'+in_vars["sdir"]+".intervals", None, in_vars["outfolder"]+'/'+in_vars["sdir"]+".realigned.bam", None])
    condition_list.append(None)
    
    textlist.append("""### Duplicate marking
            if ! [ -s $TMPDIR/$NAME.realigned.bam ] && [-s $OUTF/intermediate_files/$NAME.realigned.bam ];
            then
                cp $OUTF/intermediate_files/$NAME.realigned.bam $TMPDIR/$NAME.realigned.bam
            fi
            
            if ! [ -s $TMPDIR/$NAME.realigned.bam.bai ] && [-s $OUTF/intermediate_files/$NAME.realigned.bam.bai ];
            then
                cp $OUTF/intermediate_files/$NAME.realigned.bam.bai $TMPDIR/$NAME.realigned.bam.bai
            fi
            
                    
            if [ -s $TMPDIR/$NAME.realigned.bam ];
            then
               # clean up
               rm $TMPDIR/$NAME.sort.bam*
            
               echo Duplicate marking
               echo -e \"\\n #### duplicate marking using: java -jar $PICARD/MarkDuplicates.jar INPUT=$TMPDIR/$NAME.realigned.bam OUTPUT=$TMPDIR/$NAME.realigned.dm.bam METRICS_FILE=$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true \\n\"
               java -jar $PICARD/MarkDuplicates.jar INPUT=$TMPDIR/$NAME.realigned.bam OUTPUT=$TMPDIR/$NAME.realigned.dm.bam METRICS_FILE=$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
               cp $TMPDIR/$NAME.realigned.dm.* $OUTF/intermediate_files/
            else
               echo $TMPDIR/$NAME.realigned.bam not found >&2
               exit 1
            fi
            
            
    """.format(in_vars))
    #condition_list.append([ in_vars["outfolder"]+'/'+in_vars["sdir"]+".realigned.dm.bam", None])
    condition_list.append(None)
    
    textlist.append("""### Base quality recalibration
            if ! [ -s $TMPDIR/$NAME.realigned.dm.bam ] && [ -s $OUTF/intermediate_files/$NAME.realigned.dm.bam ];
            then
                cp $OUTF/intermediate_files/$NAME.realigned.dm* $TMPDIR/
            fi
                    
            if [ -s $TMPDIR/$NAME.realigned.dm.bam ];
            then
               echo -e \" \\n #### Base quality recalibration \\n \"
               java -jar $GATK -T BaseRecalibrator -nct {0[cpu]} --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R $REF -I $TMPDIR/$NAME.realigned.dm.bam -knownSites $DBSNP --downsampling_type NONE -o $TMPDIR/$NAME.recal_data.grp # -L $EXOME
               java -jar $GATK -T PrintReads -R $REF -I $TMPDIR/$NAME.realigned.dm.bam -BQSR $TMPDIR/$NAME.recal_data.grp -o $OUTF/$NAME.realigned.dm.recalibrated.bam # -L $EXOME
            
            else
               echo $NAME.realigned.dm.bam not found >&2
               exit 1
            fi
            
            ### Cleanup
            set +e
                   ## rm $OUTF/$NAME.realigned.bam
                   ## rm $OUTF/$NAME.realigned.bai
                   ## rm $OUTF/$NAME.realigned.dm.bam
                   ## rm $OUTF/$NAME.realigned.dm.bam.bai
            
                    rm $TMPDIR/$NAME.bam
                    rm $TMPDIR/$NAME.realigned.bam
                    rm $TMPDIR/$NAME.realigned.bai
                    rm $TMPDIR/$NAME.realigned.dm.bam
                    rm $TMPDIR/$NAME.realigned.dm.bam.bai
            set -e
            
    """.format(in_vars))
    #condition_list.append([ in_vars["outfolder"]+'/'+in_vars["sdir"]+".realigned.dm.recalibrated.bam", None])
    condition_list.append(None)
    
    textlist.append("""### GATK: Call SNPs and Indels with the GATK Haplotype Caller
            if [ -s $OUTF/$NAME.realigned.dm.recalibrated.bam ];
            then
               echo -e \"\\n #### GATK: Call SNPs and Indels with the GATK Unified Genotyper \\n\"
               #java  -jar $GATK -T UnifiedGenotyper -nt {0[cpu]} -R $REF -I $OUTF/$NAME.realigned.dm.recalibrated.bam -o $OUTF/GATK.both.raw.vcf -glm BOTH --downsampling_type NONE
               java -jar $GATK -T HaplotypeCaller -nct 1 -R $REF --dbsnp $DBSNP -I $OUTF/$NAME.realigned.dm.recalibrated.bam -o $OUTF/GATK.both.raw.vcf -L $EXOME
               #java -jar $GATK -T HaplotypeCaller -nct 1 -R $REF --dbsnp $DBSNP --emitRefConfidence GVCF  -variant_index_type LINEAR -variant_index_parameter 128000 -I $OUTF/$NAME.realigned.dm.recalibrated.bam -o $OUTF/GATK.both.raw.g.vcf -L $EXOME
               #java -jar $GATK -T GenotypeGVCFs -R $REF --dbsnp $DBSNP --variant $OUTF/GATK.both.raw.g.vcf --out $OUTF/GATK.GATK.both.raw.vcf
            
               echo -e \"\\n #### GATK: Split SNPs and Indels \\n\"
               java  -jar $GATK -T SelectVariants -R $REF --variant $OUTF/GATK.both.raw.vcf -o $OUTF/GATK.snps.raw.vcf -selectType SNP
               java  -jar $GATK -T SelectVariants -R $REF --variant $OUTF/GATK.both.raw.vcf -o $OUTF/GATK.indel.raw.vcf -selectType INDEL
            
            else
               echo $NAME.realigned.dm.recalibrated.bam not found >&2
               exit 1
            fi
            
            if [ ! -s $OUTF/GATK.snps.raw.vcf ];
            then
               echo GATK.snps.raw.vcf not found
               exit  1
            fi
            
    """.format(in_vars))
    condition_list.append(None)
    
    textlist.append("""### Filter and compare SNP calls from 2 different pipelines
            # Filtering
            mkdir -p $OUTF/SNP_Intersection
            
            java -jar $GATK -T VariantFiltration -R $REF -o $OUTF/SNP_Intersection/GATK.snps.filtered.vcf --variant $OUTF/GATK.snps.raw.vcf --mask $OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 25.0 \" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > {0[max_cov]} || GQ < 15\" --genotypeFilterName CRGg
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK Filtering VariantFiltration >&2
                exit 1
            fi

    """.format(in_vars))
    condition_list.append(None)
    
    textlist.append(
        """ ## isolate PASSed variants
        grep -E \'^#|PASS\' $OUTF/SNP_Intersection/GATK.snps.filtered.vcf | grep -v CRGg > $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
            # Correct sample names in VFC files
            # sed -i -e \"s/FORMAT\\t$NAME/FORMAT\\t$NAME-GATK/\" $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
            # sed -i -e \"s/FORMAT\\t$NAME/FORMAT\\t$NAME-SHORE/\" $OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf
            
            #if [[ ! ( -s $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf && -s $OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf  ]];
            #then
            #   echo GATK.snps.filtered.cleaned.vcf or SHORE.snps.filtered.cleaned.vcf not found >&2
            #   exit 1
            #fi
            
    """.format(in_vars))
    condition_list.append(None)
    
    textlist.append("""## Intersecting
            #  java -jar $GATK -T CombineVariants -R $REF -genotypeMergeOptions PRIORITIZE -V:SHORE $OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf -V:GATK $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -priority GATK,SHORE -o $OUTF/SNP_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING
            
            cp $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf $OUTF/SNP_Intersection/merged.vcf
            
            # Evaluation
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP  -o $OUTF/SNP_Intersection/report.all.txt --eval $OUTF/SNP_Intersection/merged.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK VariantEval >&2
                exit 1
            fi
            
    """.format(in_vars))
    condition_list.append(None)
    
    textlist.append("""## Annotate Enrichment
            $BEDTOOLS/intersectBed -a $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -b $EXOME > $OUTF/SNP_Intersection/merged.all.vcf
            
            # borrow header from GATK vcf file
            grep '^#' $OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf > $OUTF/SNP_Intersection/merged.enriched.vcf
            cat $OUTF/SNP_Intersection/merged.all.vcf >> $OUTF/SNP_Intersection/merged.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP -o $OUTF/SNP_Intersection/report.enriched.txt --eval $OUTF/SNP_Intersection/merged.enriched.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK
            STATUS="${{?}}"
            echo $STATUS
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in Annotate VariantEval GATK  >&2
                exit 1
            fi
    """.format(in_vars))
    condition_list.append(None)
    
    
    text = '\n'.join(textlist)
    script.write(text)
    
    return (textlist,condition_list)
def indel_caller_pipelne(script,in_paths,in_vars):
    text = str()



    if in_vars["indel"] == "clindel":


        text =( """### SHORE: Prepare format map.list
            ### SHORE: Prepare format map.list
            mkdir -p $OUTF/shore
            $CLINDEL convert --sort -r $REF -n 6 -g 2 -e 100 Alignment2Maplist $OUTF/$NAME.realigned.dm.recalibrated.bam $OUTF/shore/map.list.gz
            
            if [ ! -s $OUTF/shore/map.list.gz ];
            then
               echo shore/map.list.gz not found >&2
               exit 1
            fi
            
            ### Clindel call indels
            $CLINDEL qIndel -n $NAME -f $SHOREREF -o $OUTF/shore/indels -i $OUTF/shore/map.list.gz -k $DBINDEL -c 3 -d 3 -C 400 -r 0.3 -a 0.20 -b 5
            if [ ! -s $OUTF/shore/indels ];
            then
               echo shore/indels is empty >&2
               exit 1
            fi
            
            ### Clean up
            set +e
            
            rm -r $OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
            rm $OUTF/shore/Variants/ConsensusAnalysis/reference.shore
            rm $OUTF/shore/map.list.gz
            
            set -e
            
            mkdir -p $OUTF/Indel_Intersection
            
            
            # select PASSed variants
            grep -E \"^#|PASS\" $OUTF/shore/indels/indels.vcf > $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf
            
            
            # check for success of the preceeding steps
            if [[ ! -s $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf ]];
            then
               echo CLINDEL.indel.filtered.cleaned.vcf not found >&2
               exit 1
            fi 
            
            # Evaluation
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"CLINDEL\"' -selectName CLINDEL -o $OUTF/Indel_Intersection/report.all.txt --eval $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -l INFO
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK Evaluation >&2
                exit 1
            fi
            # filter Indels for enriched regions
            $BEDTOOLS/intersectBed -a $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -b $EXOME > $OUTF/Indel_Intersection/merged.all.vcf
            
            # borrow header from GATK file
            grep '^#' $OUTF/GATK.indel.raw.vcf > $OUTF/Indel_Intersection/merged.enriched.vcf
            cat $OUTF/Indel_Intersection/merged.all.vcf >> $OUTF/Indel_Intersection/merged.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"CLINDEL\"' -selectName CLINDEL -o $OUTF/Indel_Intersection/report.enriched.txt --eval $OUTF/Indel_Intersection/merged.enriched.vcf -l INFO
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in filter Indels for enriched regions >&2
                exit 1
            fi
        """.format(in_vars))
 
    elif in_vars["indel"] == "gatk":


        text = ("""###GATK Indel preparation
            mkdir -p $OUTF/Indel_Intersection
            
            java -jar -Xmx4g $GATK -T VariantFiltration -R $REF -o $OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant $OUTF/GATK.indel.raw.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > {0[max_cov]} || GQ < 15\" --genotypeFilterName LOWQ
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK VariantFiltration  >&2
                exit 1
            fi
            # select PASSed variants
            grep -v \"CRG\" $OUTF/Indel_Intersection/GATK.indel.filtered.vcf > $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
            
            # check for success of the preceeding steps
            if [[ ! -s $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf ]];
            then
               echo GATK.indel.filtered.cleaned.vcf not found >&2
               exit 1
            fi
            
            # filter Indels for enriched regions
            $BEDTOOLS/intersectBed -a $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -b $EXOME > $OUTF/Indel_Intersection/merged.all.vcf
            
            grep '^#' $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf > $OUTF/Indel_Intersection/merged.enriched.vcf
            cat $OUTF/Indel_Intersection/merged.all.vcf >> $OUTF/Indel_Intersection/merged.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"GATK\"' -selectName GATK -o $OUTF/Indel_Intersection/report.enriched.txt --eval $OUTF/Indel_Intersection/merged.enriched.vcf -l INFO
            STATUS="${{?}}"
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in filter Indels for enriched regions >&2
                exit 1
            fi
        """.format(in_vars))



    elif in_vars["indel"] == "both":
        text =("""##BOTH indel routine
            ### SHORE: Prepare format map.list
            mkdir -p $OUTF/shore
            echo $CLINDEL convert --sort -r $REF -n 6 -g 2 -e 100 Alignment2Maplist $OUTF/$NAME.realigned.dm.recalibrated.bam $OUTF/shore/map.list.gz
            $CLINDEL convert --sort -r $REF -n 6 -g 2 -e 100 Alignment2Maplist $OUTF/$NAME.realigned.dm.recalibrated.bam $OUTF/shore/map.list.gz
            
            if [ ! -s $OUTF/shore/map.list.gz ];
            then
               echo shore/map.list.gz not found >&2
               exit 1
            fi
            
            ### Clindel call indels
            $CLINDEL qIndel -n $NAME -f $SHOREREF -o $OUTF/shore/indels -i $OUTF/shore/map.list.gz -k $DBINDEL -c 3 -d 3 -C 400 -r 0.3 -a 0.20 -b 5
            
            if [ ! -s $OUTF/shore/indels ];
            then
               echo shore/indels is empty >&2
               exit 1
            fi
            
            ### Clean up
            set +e
            
            rm -r $OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
            rm $OUTF/shore/Variants/ConsensusAnalysis/reference.shore
            rm $OUTF/shore/map.list.gz
            
            set -e
            
            
            ##GATK indel section
            mkdir -p $OUTF/Indel_Intersection
            
            # filter GATK variants (Clindel variants are produced pre-filtered)
            java -jar -Xmx4g $GATK -T VariantFiltration -R $REF -o $OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant $OUTF/GATK.indel.raw.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > {0[max_cov]} || GQ < 15\" --genotypeFilterName LOWQ
            STATUS="${{?}}"
            echo Clindel variants produced pre-filtered $STATUS
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK filter GATK variants Clindel variants are produced pre-filtered >&2
                exit 1
            fi
            # select PASSed variants
            grep -v \"CRG\" $OUTF/Indel_Intersection/GATK.indel.filtered.vcf > $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
            
            grep -E \"^#|PASS\" $OUTF/shore/indels/indels.vcf > $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf
            
            
            
            # Correct sample names in VFC files --- not necessary for tool switching
            sed -i -e "s/FORMAT\t$NAME/FORMAT\t$NAME-GATK/" $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
            sed -i -e "s/FORMAT\t$NAME/FORMAT\t$NAME-CLINDEL/" $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf
            sed -i -e "s/##FORMAT=<ID=GQ,Number=1,Type=String,Description=/##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=/" $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf
            
            # check for success of the preceeding steps
            if [[ ! ( -s $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf && -s $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf ) ]];
            then
               echo GATK.indel.filtered.cleaned.vcf or CLINDEL.indel.filtered.cleaned.vcf not found >&2
               exit 1
            fi
             
            # Intersecting
            java -Xmx5g -jar $GATK -T CombineVariants -R $REF -genotypeMergeOptions PRIORITIZE -V:GATK $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -V:CLINDEL $OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -priority GATK,CLINDEL -o $OUTF/Indel_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING
            STATUS="${{?}}"
            echo GATK intersecting :  $STATUS
            if [ "$STATUS" -gt 0 ];
            then
                echo GATK Intersecting >&2
                exit 1
            fi
            # Evaluation
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"CLINDEL\"' -selectName CLINDEL -select 'set==\"GATK\"'  -selectName GATK_CLINDEL -o $OUTF/Indel_Intersection/report.all.txt --eval $OUTF/Indel_Intersection/merged.vcf -l INFO
            STATUS="${{?}}"
            echo GATK Evaluation : $STATUS
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK Evaluation >&2
                exit 1
            fi
            # filter Indels for enriched regions
            $BEDTOOLS/intersectBed -a $OUTF/Indel_Intersection/merged.vcf -b $EXOME > $OUTF/Indel_Intersection/merged.all.vcf
            
            grep '^#' $OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf > $OUTF/Indel_Intersection/merged.enriched.vcf
            cat $OUTF/Indel_Intersection/merged.all.vcf >> $OUTF/Indel_Intersection/merged.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"CLINDEL\"' -selectName CLINDEL -select 'set==\"GATK\"' -selectName GATK -o $OUTF/Indel_Intersection/report.enriched.txt --eval $OUTF/Indel_Intersection/merged.enriched.vcf -l INFO
            STATUS="${{?}}"
            echo GATK  filter indels $STATUS
            if [ "$STATUS" -gt 0 ];
            then
                echo Error in GATK filter Indels for enriched regions >&2
                exit 1
            fi
        """.format(in_vars))
    else:
        print("The indel caller is not recognized: should be gatk - clindel - both")
        raise
    script.write(text)
    return text

def fusevariants_pipelne(script,in_paths,in_vars):
    text="""#FUSEVARIANTS    
        # fuse indel and snp calls into one file
        java -Xmx5g -jar $GATK -T CombineVariants -R $REF --variant $OUTF/SNP_Intersection/merged.enriched.vcf --variant $OUTF/Indel_Intersection/merged.enriched.vcf -o $OUTF/all_variants.vcf -U LENIENT_VCF_PROCESSING
        STATUS="${{?}}"
        if [ "$STATUS" -gt 0 ];
        then
            echo Error in fusevariants >&2
            exit 1
        fi
    """.format(in_vars)
    script.write(text)
    return text


def cleanup(script,in_paths,in_vars):
    text="""CLEANUP
        rm -r $OUTF/intermediate_files/*
        """
    script.write(text)
    return text

# End of prediction_support_functions.py

