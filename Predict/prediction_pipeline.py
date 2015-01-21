#!/usr/bin/python -u
import readline
import os
import sys
import prediction_support_functions
if not os.environ.get('DRMAA_LIBRARY_PATH'):
    print "Adding the DRMAA library path to the environment"
    os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
#import ../pipeline_control/qsubclass
#import ../pipeline_control/pipeline_element
import pickle
import subprocess
import datetime
import imp
curpath= os.path.realpath(__file__)
curpath = curpath.split('/')
curpath = '/'.join(curpath[:-2])
curpath += '/'

qsubclass = imp.load_source('qsubclass', os.path.abspath(curpath + 'pipeline_control/qsubclass.py'))
pipeline_element = imp.load_source('pipeline_element', os.path.abspath(curpath +'pipeline_control/pipeline_element.py'))

python_path  = '/software/so/el6.3/PythonPackages-2.7.6/bin/python'
pipe_script  = os.path.abspath(curpath +'/pipeline_control/pipeline_element.py')
qsub_script = os.path.abspath(curpath +'pipeline_control/qsubclass.py')
## main execution

## Init variables that will be used later
infolder        = ""
outfolder       = ""
qsub_name       = "predict.sh"
max_cov         = 400 
nameStart       = 1
nameLength      = 6
firstreadextension = ""
secondreadextension= ""
indelcaller     = "clindel"
fusevariants    = 'no'
cpu             = 4
mem             = 20
help_c          = 0
cursor ="\n>"
qoptions_def = False
sample_list = None

(commandLine,inargs) = prediction_support_functions.parse_commandline()

if commandLine:
    parsed_config = prediction_support_functions.process_config_file(inargs["config"])
    infolder = inargs["infolder"]+'/'
    outfolder = inargs["outfolder"]+'/'
    qsub_name = inargs["qsubname"]
    nameStart = inargs["namestart"]
    nameLength = inargs["namelength"]
    firstreadextension = inargs["firstreadextension"]
    secondreadextension = inargs["secondreadextension"]
    if inargs.get("max_coverage",False):
        max_cov=inargs["max_coverage"]
    if inargs.get("indelcaller",False) :
        indelcaller=inargs["indelcaller"] 
    if inargs.get("cpu",False) :
        cpu = inargs["cpu"] 
    if inargs.get("mem",False):
        mem = inargs["mem"]
    if inargs.get("fusevariants",False):
        fusevariants = inargs["fusevariants"]
    if inargs.get("qoptions",False):
        qoptions_def = qsubclass.parse_command_options(inargs["qoptions"],outfolder,qsub_name)
        qoptions_input = True
    else:
        qoptions_input = False
    if inargs.get("sample_list",False):
        sample_list = inargs["sample_list"]
else:
    ### Gather and parse input file
    parsed_config   = prediction_support_functions.parse_config_file()
    
    #load infolder
    infolder  = prediction_support_functions.get_path("Please insert the path where your input files are located:%s"%(cursor),0)
    infolder +='/'
    #outfolder
    outfolder = prediction_support_functions.get_path("Please insert the path where you want your results to be stored:%s"%(cursor),1)
    outfolder +='/'
    
    #name_parameters
    (nameStart,nameLength) = prediction_support_functions.get_namelength(infolder)
    
    #extension parameters
    (firstreadextension,secondreadextension) = prediction_support_functions.get_extensions(infolder)
        
    #cpu params
    print "Now the cpu parameters: default is cpu = 4 and 20G of ram"
    print "If it's ok for you hit return"
    cpu = prediction_support_functions.parse_input("How many cpus? Enter for default = %d %s"%(cpu,cursor),1,32,cpu)
    print(">CPU cores : %d"%(cpu))
    mem = prediction_support_functions.parse_input("How many GB RAM? Enter for default = %d %s"%(mem,cursor),1,100,mem)
    print(">RAM : %d GB" %(mem))
    #test for indel caller
    indelcaller = prediction_support_functions.indel_parser(indelcaller)
    #max coverage
    max_cov = prediction_support_functions.parse_input("Max Coverage? Enter for default = %d %s"%(max_cov,cursor),1,1000,max_cov)
    print(">Max Coverage : %d x" %(max_cov))
    # fusevariants
    fusevariants = prediction_support_functions.fuse_parser(fusevariants)
    #finally load the qsubname
    print("Please insert a name for the qsub script to run in the cluster.  \n"
          "The default value is %s \n Hit return to leave as it is"
          %(qsub_name))
    tmp = raw_input(">")
    if len(tmp):
        if (tmp[0].isdigit()):
            qsub_name="s"+tmp
            print("qsub should start with a letter. Changed to : %s"%(tmp))
        else:
            qsub_name=tmp
    print("The chosen script name is %s"%(qsub_name))
    

#--------------------------------------------------------
# Now start the pipeline text file for real
#--------------------------------------------------------


#Run a scan in the input folder to find for all the samples
# A Sample is each pair of XXX.firstreadextension and XXX.secondreadextension
# Failing to meet that criteria will result in nothing
# For each sample it builds a XXX directory inside outfolder: outfolder/XXXX/
# Inside each outfolder/XXX/qsubname.sh

#pseudocode:
#   Find samples
#       For each sample
#       Build the output subfolder
#       Fill the pipeline with the correct parameters
#   End.

samples = prediction_support_functions.find_samples(infolder,nameStart,nameLength,firstreadextension,secondreadextension,sample_list)


#qclass = qsubclass.disclaimer()
qclass = True

if qclass:
    qlist =list()
    logfile = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + "qsub_log.log"
    if qoptions_def == False:
        qoptions_def,logfile = qsubclass.getOptions()
    else:
        print "Automatic Logfile:%s"%logfile
    print qoptions_def

#for i in range(3,4):
for full_sample in samples:
    qoptions = qoptions_def
 #   full_sample = samples[i]
    sample = full_sample[nameStart-1:nameLength]
    #mkdir
    prediction_support_functions.make_dirs(outfolder,sample)
    #buildpipe
    sample_info = {"sdir":sample,"sname":full_sample, "qsubname" : qsub_name, "first_e":firstreadextension,"second_e":secondreadextension,
                   "outfolder":outfolder, "cpu":str(cpu),"mem":str(mem), "indel":indelcaller, "max_cov":max_cov,"infolder":infolder}
    script_name =outfolder+sample+'/'+qsub_name  
    
    with open(outfolder+sample+'/'+qsub_name,'w')  as script: # Trying to create a new file or open one
        pipe = list()
        if qclass :
            if qoptions_input==False:
                qoptions['-e'] = list()
                qoptions['-e'].append(outfolder+'/'+sample+'/')
                qoptions['-o'] = list()
                qoptions['-o'].append(outfolder+'/'+sample+'/')
                qoptions['-N'] = [qsub_name]
                if qoptions.get('-l',False):
                    qoptions['-l'].append("virtual_free=%dG"%mem)
                else:
                    qoptions['-l'] = list()
                    qoptions['-l'].append("virtual_free=%dG,h_rt=51600"%mem)
                qoptions['-pe'] = list()
                qoptions['-pe'].append("smp %d"%cpu)
            command_pipe = python_path +' '+ pipe_script+ ' ' + script_name+'.pipe ' + logfile
            #print command_pipe
            qlist.append(qsubclass.qsubCall(command_pipe,qoptions,list(),logfile))
            
        ## Here now is where the pipeline is written in the "script" file        
        header,env_var = prediction_support_functions.write_header(script,parsed_config,sample_info)
        p_element = pipeline_element.pipeline_element(header,"Header writing")
        p_element.set_error("Error in header executions Please refer to SGE job error file")
        pipe.append(p_element)
        
        text = prediction_support_functions.seq_pipeline(script,parsed_config,sample_info)
        p_element = pipeline_element.pipeline_element(env_var + text,"Alignment etc ")
        p_element.set_error("Error in main body Please refer to SGE job error file")
        pipe.append(p_element)
        
        text = prediction_support_functions.indel_caller_pipelne(script,parsed_config,sample_info)
        p_element = pipeline_element.pipeline_element(env_var+ text,"Indel caller")
        p_element.set_error("Error in indel caller Please refer to SGE job error file")
        pipe.append(p_element)
        
        
        # Fusepipe
        if fusevariants=="yes":
            text = prediction_support_functions.fusevariants_pipelne(script,parsed_config,sample_info)
            p_element = pipeline_element.pipeline_element(env_var+text,"Fusevariant")
            p_element.set_error("Error in fusevariant executions Please refer to SGE job error file")
            pipe.append(p_element)
        pipeline_element.save_pipeline(pipe,script_name+'.pipe')
 
if qclass:
    try:
        print "Now I'm saving the queue list in %s file"%(logfile+'.qlist')
        with open(logfile +'.qlist','wb') as outfile:
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
    nohup_script = ("""import imp;import os;os.environ['DRMAA_LIBRARY_PATH'] = '/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0';"""+
    "qsubclass = imp.load_source('qsubclass','%s');qsubclass.run_qlist('%s.qlist') "%(qsub_script,logfile))
   # print nohup_script
    #os.system('nohup '+ python_path+' -c '+ '"'+nohup_script+'" ')
    qsubclass.run_qlist(logfile+'.qlist')
    print "Finished, have a nice day :)"

    