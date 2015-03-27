#!/usr/bin/env python
# Try to build a class of qsub call objects
# They must have attributes and defaults to call the qsub
# Then to call qsub with specific changes
# Reruns
# Status check
# Printing method
# Print to log_file
""" module documentation here """

import subprocess
import os
import sys
import datetime
import pickle
import time
import readline
if not os.environ.get('DRMAA_LIBRARY_PATH'):
    os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
    import drmaa
else:
    import drmaa
import datetime    
#try :
#    import drmaa
#except:
#    subprocess.call('export DRMAA_LIBRARY_PATH=/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0',shell=True)

    
qsub_options = ['-q','-m','-M','-l','-v','-e','-o','-N','-b','-shell','-pe']
qsub_flags   = ['-hard','-cwd']

def disclaimer():
    ''' beta test disclaimer - now not used'''
    qclass = False
    possibilities =['yes','no']
    token = 1
    msg = """
    Hi, I'm a beta function to automatize execution of the produced scripts.
    I will execut the scripts, re-execute them in case of error and produce
    a log file with a synthetic description of what happened...
    well, I'll be more complete in a recent future thanks to your feedback!
    I will execute automatically the scripts for you but you should not
    close this terminal session to ensure a full execution of the scripts on
    the SGE.
    Whould you like to try me? Please write 'yes' or 'no'
    """
    while token:
        instr = raw_input(msg)
        if instr in possibilities:
            token = 0
        else:
            msg ="Please write 'yes' or 'no'" 
    if instr =='yes':
        return True
    else:
        return False
    

def getOptions(outfolder = "",qsub_name =""):
    """ This function helps the user in gathering the options for SGE calls """
    msg ="""
    Hi!, Now I need two things before starting the execuiton:
        1) A logfile name which will be stored in this folder where I will put the
            progress
        2) A set of options for the SGE qsub call, to resemble your normal course
            of actions:
            
        Let's start with the logfile.
        Please insert a name for the logfile:
    """
    logfile = ""
    while len(logfile)==0:
        logfile = raw_input(msg)
    
    msg = """
    Thanks.
    Now please insert the execution options in the following format:
    >  -m value,-g value,-q value
    Please help me and don't add any spaces before or after or I'll fail...
    
    An example could be:
    -m abe,-M name.surname@crg.es,-hard !,-q so-el6,-v CUSTOM_EMAIL=yes
    
    Ah btw, you can leave this empty and will run with default SGE options
    """
    token = 1
    opts = dict()
    while token:
        try:
            instr = raw_input(msg)
            if len(instr)>0:
                opts = parse_command_options(instr,outfolder,qsub_name)
                token = 0
            else:
                opts = dict()
                token = 0
        except:
            print"There is an error, try to change the string format or leave it empty"
    
    return opts,logfile

def parse_call(command,options):
    """Produces the qsub call line
        It needs command = full path to shell script
        options : dictionary with {-option,value}
        The accepted options are in qsub_option global variable"""
    ##Later on, here could
    tmp = list()
    commandline = ""
    for x in qsub_flags:
        if(options.get(x,False)):
            for i in (options[x]):
                tmp.append(' '+x+' ')
    for x in qsub_options:
        if(options.get(x,False)):
            for i in (options[x]):
                tmp.append(' '+x+' '+''.join(i))

    commandline+= ''.join(tmp)
    commandline+=" -b y -shell y -cwd"#%(command)
    commandline = '-hard ' + commandline
    return commandline

def parse_command_options(instring,outfolder,qsub_name,old_dict=dict()):
    """Function to produce the options dictionary from a formatted string.
    The format is "-x value, -y value, -z value" and ideally could be sent from command line in new modules-functions"""
#   input is a string "-x value, -y value, -z value"
    
    out = dict()
    
    if len(outfolder)>0 and len(qsub_name)>1:
        out['-e'] = list().append(outfolder+'/')#+qsub_name
        out['-o'] = list().append(outfolder+'/')#+qsub_name
    
    fields = instring.split(',')
    for field in fields:
        field = field.strip()
        try:
            (key,val) = field.split(' ')
            if key in qsub_options:
                if out.get(key,False):
                    out[key].append(''.join(val))
                else:
                    out[key] = list()
                    out[key].append(''.join(val))
        except:
            if field in qsub_flags:
                if out.get(field,False):
                    out[field].append('')
                else:
                    out[field] = list()
                    out[field].append('')
    
    old_dict.update(out)
       
    
    return old_dict

def run_qlist(filename,qopts=None):
    """It loads a file where a list of qsubCall instances have been stored
        It checks the status and run each element if needed
        It takes care of rerun, cleanup, status update
        It finally stores the updated list in the very same file for a future execution"""
    qlist = list()
    try:
        with open(filename,'rb') as infile:
            qlist = pickle.load(infile)
    except :
        print("The chosen file does not exist or does not contain the qclass list \n : %s",filename)
    try:
        reruns = qlist[0].reruns
    except:
        reruns = 0
        
        
    #build sessions for drmaa
    with drmaa.Session() as session:    
        #session.initialize()
        jt = [None]*len(qlist)
        jobid = [None]*len(qlist)
        
        for i in range(reruns):
            
            #Run once all of them
            for ql in range(len(qlist)):

                if qopts != None:
                    if qlist[ql].options != None:
                        #qlist[ql].qoptions = parse_command_options(qopts,"","")
                        qlist[ql].options = parse_command_options(qopts,"","",qlist[ql].options)
                    else:
                        
                        qlist[ql].options = parse_command_options(qopts,"","")

                session,jt[ql],jobid[ql] = qlist[ql].just_launch(session)
                
            #Save file with jobid information
            with open(filename,'wb') as outfile:
                pickle.dump(qlist, outfile)
                print "saved the updated qlist with jobid"
            #Capture the outcome
            for ql in range(len(qlist)):
                #print ql
                qlist[ql].hanlde_and_close(session,jt[ql],jobid[ql])
            
        try:    
            with open(filename,'wb') as outfile:
                pickle.dump(qlist, outfile)
        except :
            print('error saving pipeline list to file')

def build_qlist(file_list, options,cleanup_list,logfile):
    """routine to build and store a list of qsubCall instances"""
    return None

class qsubCall :
    # Fields:
    # command  : shell filename
    # options  : qsub options as a dictionary ['m':12,'cpu':3... will be limited and checked before therun]
    # status   : 0...max_rerun or -1 if ended well
    # logfile  : where to print the log details
    # cleanup_list
    """The qsubCall is a class built of objects needed to run a shell script in the sun grid engine computation system.
        It takes care of initializing a call, running, re-running if there is some problem, and cleaning up unwanted files
        at the end of the execution.
        It can be used jointly with the ediva pipeline code or also outside it thanks to the provided functions build_qlist and run_qlist
    """
    
    def __init__(self,command,options,cleanup_list,logfile):
        """Initialization:
            It needs:
            command = shell file to be executed with the full path
            options = qsub options in a dictionary format
            cleanup_list= list with complete path to those files to be removed after execution
            logfile  = a logfile where to store the execution progress: Recommended to use the same logfile for those files executed together
            reruns  = number of attempts before calling an error : default = 2
            status  = execution status: 0 = to be executed, -1 = executed successfully, > 0 = number of failed execution attempts
            
        """
        self.command = str(command)
        
        self.options = dict(options)
        self.cleanup_list = list(cleanup_list)
        self.logfile = str(logfile)
        self.reruns  = 1
        self.status  = 0
        self.jobid   = False
        
    #Print out the details of the string
    def __str__(self):
        """Prints a description of the class element"""
        return "qsub_command %s with options %s and logfile %s "%(self.command,self.options,self.logfile)
    
    
    def cleanup(self):
        """Deletes the cleanup_list files"""
        print("Cleanup")
        for i in self.cleanup_list:
            print i
        return None
    
    
    def run(self):
        """ Not really used at the moment
            Runs the scripts #times,
            It takes care of storing the progress in a logfile
            It takes care of cleaning up after the execution
            It updates the execution status of the element"""
        optionline = parse_call(self.command,self.options)
        with  open(self.logfile,'a') as logfile:
            #open logfile
            for i in range(self.reruns):           
                if self.status >=0:
                    self.status+=1
                    try:
                        logfile.write(str(datetime.datetime.now()))
                        logfile.write("  Running %s attempt number %d\n"%(optionline,i+1))
                        s = drmaa.Session()
                        s.initialize()
                        jt = s.createJobTemplate()
                        jt.remoteCommand = self.command
                        jt.joinFiles=False
                        jt.nativeSpecification  =optionline
                        jobid = s.runJob(jt)
                        s.deleteJobTemplate(jt)
                        s.exit()
                        if retval.hasExited and retval.exitStatus ==0 :
                            self.status = -1
                            token = 0
                        else: 
                            print "ERROR"
                            print retval.hasExited
                            print retval.exitStatus
                            raise SystemError
                    except Exception,e:
                        print "ERROR"
                        print str(e)
                        logfile.write(str(datetime.datetime.now()))
                        logfile.write(" %s" %str(e))
                    finally:
                        self.cleanup()
                elif sellf.status >= self.reruns:
                    print "couldnt complete the task in the max number of attempt"
                else:
                    print("Completed")
                    print self.status
                    logfile.write(str(datetime.datetime.now()))
                    logfile.write("  Completed\n")
                    break
        return self.status
    #!/usr/bin/env python
    def just_launch(self,s):
        ''' Launches the qsub call '''
        optionline = parse_call(self.command,self.options)
        jt = None
        jobid = 0
        session = None
        with  open(self.logfile,'a') as logfile:
            self.status=0
            if  self.status >=0 and self.status<self.reruns:
                self.status+=1
                logfile.write(str(datetime.datetime.now()))
                logfile.write("  Running %s attempt number %d\n"%(optionline,self.status))
                jt = s.createJobTemplate()
                jt.remoteCommand = self.command
                jt.joinFiles=False
                jt.nativeSpecification=optionline
                self.jobid = s.runJob(jt)
                
        return(s,jt,self.jobid)
    
    def hanlde_and_close(self,s,jt,jobid):
        with  open(self.logfile,'a') as logfile:
            retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            #retval = s.wait(jobid, 1)
            try   :
                #from here on we need: [self,jobid,jt,s]                
                print('Job: {0} submitted with status {1} and app status {2}'.format(retval.jobId, retval.hasExited,retval.exitStatus))
                s.deleteJobTemplate(jt)
                if retval.hasExited and retval.exitStatus ==0 :
                    self.status = -1
                else: 
                    raise SystemError
            except Exception,e:
                print "ERROR in job : {0}".format(retval.jobId)
                print "SGE call exited with status {0}".format(retval.hasExited)
                print "The application exited with {0}".format(retval.exitStatus)
                print str(e)
                logfile.write(str(datetime.datetime.now()))
                logfile.write("ERROR in job : {0}\n".format(retval.jobId))
                
            finally:
                pass
                #self.cleanup()
        return None
    


if __name__=='__main__':
    import subprocess
    import os
    import sys
    import datetime
    import pickle
    import time
    import readline
    if not os.environ.get('DRMAA_LIBRARY_PATH'):
        os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
        import drmaa
    else:
        os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
        print "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
        import drmaa
        
    #Idea is to load a file in argv[1] and run it:
    filename = os.path.abspath(sys.argv[1])
    try:
        qopts = sys.argv[2]
        if len(qopts)>0:
            run_qlist(filename,qopts)
    except:
        run_qlist(filename)
        
    