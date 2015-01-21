#!/usr/bin/env python
import pprint
import argparse
import os.path
import re
import readline
import subprocess
curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
predictpath = '/'.join(curpath[:-2])
predictpath+='/Predict/'
curpath = '/'.join(curpath[:-1])
curpath += '/'

python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'


## Collect the script location and use it for the  execution of scripts
def parse_config_file(msg):
    '''
    This function helps the user to specify where a config file is, and if it does not exist, creates it
    Depending if the msg is setup or family, the generated file will be different
    '''
    token =1
    while token:
        answer = raw_input("Do you already have a %s config file ? please print 'y' or 'n'\n"%msg)
        if answer == 'y':
            token = 0
        elif answer =='n':
            while token:
                readline.parse_and_bind("tab: complete")
                readline.set_completer_delims(' \t\n;')
                text = raw_input("please insert the full path to your new %s config file\n >"%msg)
                if len(text)>0:
                    if text[0] =='~':
                        text = os.path.expanduser("~") + "/"+ text[1:] 
                    else :
                        pass                                            
                if (msg == ' overall '):
                    cmd =  ("%s %s/setup.py --newconfig %s"%(python_path,predictpath,text))
                    cmd1 =("%s/setup.py"%(predictpath))
                    cmd2 = '--newconfig %s'%text
                else:
                    cmd =  ("%s %s/collect_family_info.py --outfile %s"%(python_path,curpath,text))
                    cmd1= (" %s/collect_family_info.py"%(curpath))
                    cmd2 = '--outfile %s'%text
                print cmd
                p = subprocess.Popen("%s"%cmd,shell=True)
                
                stdout,stderr = p.communicate()
                p.wait()
                #setup.main(True,text)
                token = 0
    token =1
    prompt_0 = "Please insert the %s config file \n >"%msg
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
                prompt = "Error: the chosen path does not exist. Please correct it" + prompt_0 + " \nPath: " + text
        else :
            prompt = "Please insert something: \n " + prompt_0
            
    print (">The chosen configuration file is %s" % (text))
    #Here we load the file and check for all the parameters:
    #init default variables:
    return text

def get_path(prompt,outfolder=True):
    '''
    This function gathers an existing path as outfolder variable
    '''
    prompt_0 = prompt
    token =1
    while token : 
        #outputs a prompt, parse it to a full path and  checks that it works
        readline.parse_and_bind("tab: complete")
        readline.set_completer_delims(' \t\n;')
        text  = raw_input(prompt+"\n>")
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

def yes_no(text,state):
    '''
    Limits the options to yes or no when reading and input
    '''
    possibilities = ["yes", "no"]
    results = [x for x in possibilities if x.startswith(text)]
    return results[state]

def yes_no_alternative(fuse_def,msg):
    '''
    used to catch the user preferences for binary variables, it returns True or False
    '''
    print msg
    token = 1
    out_str = None
    dictionary =["yes","no"]
    readline.parse_and_bind("tab: complete")
    readline.set_completer(yes_no)
    while token:
        out_str = raw_input("%s  Default is %s %s"%(msg,fuse_def, "\n>"))
        if len(out_str)<1:
            out_str = fuse_def
            token =0
        elif (out_str in dictionary):
            token = 0
    print (">Selected: %s"%(out_str))
    if out_str =='yes':
        return True
    else:
        return False

def family(text,state):
    '''
    Used to limit the user choices for family type
    '''
    possibilities = ["trio", "family"]
    results = [x for x in possibilities if x.startswith(text)]

    return results[state]

def family_type(def_,msg):
    '''
    Guides the user to choose the family type
    '''
    print msg
    token = 1
    out_str = None
    dictionary =["trio","family"]
    readline.parse_and_bind("tab: complete")
    readline.set_completer(family)
    while token:
        out_str = raw_input("%s  Default is %s %s"%(msg,def_, "\n>"))
        if len(out_str)<1:
            out_str = def_
            token =0
        elif (out_str in dictionary):
            token = 0
    print (">Selected: %s"%(out_str))
    return out_str

def inherits(text,state):
    '''
    Limits the options to dictionary when reading and input
    '''
    possibilities = ["dominant_inherited", "dominant_denovo","recessive","Xlinked","compound"]
    results = [x for x in possibilities if x.startswith(text)]
    return results[state]

def inheritance_parse():
    '''
    Guides the user to choose the inheritanc types
    '''
    inherit_list = list()
   
    msg =  "Please choose an inherit type among the available ones . hit tab to show them"
    token = 1
    out_token=1
    out_str = None
    dictionary =["dominant_inherited", "dominant_denovo","recessive","Xlinked","compound"]

    while out_token:
        token =1
        while token:
            readline.parse_and_bind("tab: complete")
            readline.set_completer(inherits)
            out_str = raw_input("%s \n>"%msg)
            if len(out_str)<1:
                out_str = ""
                token =0
                out_token = 0
            elif (out_str in dictionary):
                inherit_list.append(out_str)
                token = 0
        print (">Selected: %s"%(out_str))
        if out_token ==0 and len(inherit_list)==0:
            print "Please choose at least one "
            out_token =1
    return inherit_list

def parse_args():
    '''
    Parses the initial arguments, if they are some that are not specified, it will automatically ask for them
    -If no config file or family config file are present, it allows the user to interactively build them on the fly
    -Anyways it is adivisable to have them ready beforehand.
    
    '''
    parser = argparse.ArgumentParser(description = 'Create a file that will be run in the cluster to fully run the prioritization in eDiVa.')
    parser.add_argument('--config',  type=argparse.FileType('r'), dest='config',  required=False, help='The config file you created using setup.py, containing all paths of GATK, SAMTOOLS, etc.')
    parser.add_argument('--family',  type=argparse.FileType('r'), dest='family',  required=False, help='A file containing all paths to vcf and bam files, including the affection status. Can be created with collect_family_info.py')
    parser.add_argument('--outfolder', type=str, dest='outfolder', required=False, help='The folder where the output data should be written to.')
    parser.add_argument('--inheritance', choices=['dominant_denovo', 'dominant_inherited', 'recessive', 'Xlinked', 'compound'], dest='inheritance', required=False, action = 'append', help="""choose a inheritance model [required]
    This option can be called multiple times.
    dominant_inherited: used for families
    dominant_denovo: apply to novel variants seen in the affected individuals
    
    recessive: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
    Xlinked: used for X linked recessive variants in trios only
    compound: detect compound heterozygous recessive variants
    """)
    parser.add_argument('--familytype', type=str, dest='familytype', required=False, help="choose if the data you provide is a trio or a larger family")
    parser.add_argument('--multisample', required=False, action='store_true', help="if your input variants have been called using GATK multi sample calling, you should toggle this option.")
    parser.add_argument('--qsubname', type=str, dest='qsub_name', required=False, help='name of the script you want to start later on')
    parser.add_argument('--jobname', type=str, dest='job_name', required=False, default='prioritize', help='name of the job that will be run on your cluster [default: prioritize]')
    parser.add_argument('--force', action = 'store_true', help='Enable eDiVa to overwrite old output')
    parser.add_argument('--qoptions', type=str, dest='qoptions', required=False, help='--qoptions \t qsub options with this format: -X val,-Y val,-Z, -XX value\n')
    args = parser.parse_args()
    #now we collect the variables
    if args.config ==None  :
        #config
        args.config = open(parse_config_file(" overall "),'r')
    if args.family == None:
        #family config
        args.family = open(parse_config_file(" family "),'r')
    if args.outfolder == None:
        #outfolder
        args.outfolder = get_path("Please insert the outfolder directory")
    if args.familytype == None:
        #familytype
        args.familytype = family_type('family',"Please choose which family type: 'trio' - 'family'\n")
    if args.multisample == None:
        #multisample
        args.multisample = yes_no_alternative('no',"Please write if the -multisample option was used for the sample analisys: 'yes' - 'no'\n")
    if args.qsub_name == None:
        #qsubname
        args.qsub_name =  raw_input("Please insert the desired qsub name\n>")
    if args.inheritance == None:
        #inheritance
        args.inheritance = inheritance_parse()
    if args.job_name == None:
        #jobname
        args.job_name = raw_input("Please insert the desired job name\n>")
    if args.force == None:
        #force
        args.force = yes_no_alternative('no',"Please write if you want to enable the -force option: 'yes' - 'no'\n")       
        
    
    # does the output folder exist?
    if not os.path.exists(args.outfolder):
        print "Creating output folder %s" % args.outfolder
        os.makedirs(args.outfolder)
    
    # check if pipeline was run in that folder already
    if os.path.isfile("%s/combined.variants.vcf" % args.outfolder) and not args.force:
        print("The output files of prior run still exist in the specified folder. Please use the --force option or choose another folder.")
    
    return args

