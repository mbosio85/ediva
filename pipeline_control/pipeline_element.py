#!/usr/bin/env python
import subprocess
from datetime import datetime
import pickle
import sys
import os

class pipeline_element :
    # Fields:
    # command  : command to be run by the shell
    # status   : Has it still to be run? How many times did it crash? Error status
    # err_msg  : Error message to output
    # alias    : easier name to remember when reading the logfile.    
    
    """
    The pipeline_element is a class containing a single element within a pipeline
    The idea behind it is to have a better control over a pipeline execution and
    the current status within a python wrapper.
    The functions it should provide are a constatn status checkup, quick resume and
    efficient log generation [external functions provided for that]
    """
    
    def __init__(self):
        '''empty initialization function'''
        self.command = ""
        self.status  = 0 # 0 = not run yet 1 = run and stopped -1 = run with success
        self.err_msg = ""
        self.alias = self.command
        self.condition = None
        
    def __init__(self,command,alias = None):
        '''Initialization with a command and optional Alias'''
        self.command = command
        self.status = 0
        self.condition = None
        self.err_msg = ""
        if alias != None:
            self.alias = alias
        else:
            self.alias = self.command
    def set_error(self,msg):
        self.err_msg = msg
    def set_alias(self,msg):
        self.alias = msg
    def set_status(self,status):
        ''' set the current status'''
        self.status = status
    def set_condition(self,condition):
        ''' set the output control condition'''
        self.condition=condition
    def get_info(self):
        ''' extract a info list for the function'''
        return self.command,self.status,self.err_msg,self.alias
    def __str__(self):
        '''Print management'''
        msg = "Pipeline Element: "+self.alias + " with status : "+ str(self.status)
        return msg
    def check_condition(self ):
        if self.condition!=None: #not None
            try:
                for i in range(0,len(self.condition),2):
                    filename = self.condition[i]
                    mod_time = self.condition[i+1]
                    #print filename
                    #print mod_time
                    if mod_time != None:
                        cur_time = os.path.getmtime(filename)
                        #print cur_time
                        if mod_time != cur_time:
                            print '%s has been modified since last execution, reset to 0'%filename
                            self.status=0
                    else:
                        self.status =0
                    
            except:
                print 'error in checking the condition in %s'%self.alias
                self.status = 0
                raise
    def update_condition(self):
        try:
            for i in range(0,len(self.condition),2):
                filename = self.condition[i]
                #print filename
                self.condition[i+1] = os.path.getmtime(filename)
                #print self.condition[i+1]
        except:
            print 'error in updating the condition in %s'%self.alias
            raise    
        
    
    
    def run(self):
        '''
        Run the pipeline command with a subprocess call
        Checks the correct execution
        If not, it updates the error message with a timestamp '''
        with open("err_p.log",'w+') as er , open("out_p.log",'w+')as ou:    
            if self.status>=0:
                try:
                    #print(self.command)
                    
                    #a = subprocess.check_output(self.command,shell=True,stderr=er)
                    #a = subprocess.check_output(self.command,shell=True)
                    a = subprocess.call(self.command,shell =True,stdout=ou,stderr=er)
                    #If success:
                    
                    
                    if a ==0:
                        self.status = -1
                    else:
                        print self.command
                        print 'error : status %d'%a
                        self.status+=1
                        raise
                except:
                    self.status = 1
                    self.err_msg = str(datetime.now()) + " "  + self.alias + " ended with : " + self.err_msg
                    print "There is an error with %s"%self.err_msg
                    sys.stderr.write("Here will go the error file\n\n")
                    er.seek(0)
                    for line in er:
                        sys.stderr.write(line)
                        raise Exception
                    pass
                finally:
                    #print ou
                    er.seek(0)
                    for line in er:
                        sys.stderr.write(line)
                    print "\t> "+ self.alias
                    #print self.status
                    pass
        return None


def save_pipeline(pipe,filename):
    ''' saves a list of pipe_element in a binary file'''
    try:
        #print "\n\nNow I'm saving the queue list in %s file"%(filename)
        with open(filename,'wb') as outfile:
            pickle.dump(pipe, outfile)
            #print outfile
            return 0
    except :
        print('error saving pipeline list to file')
        return 1    
def load_pipeline(filename):
    ''' load a binary file and returns a list of pipe_elements'''
    try:
        with open(filename,'rb') as infile:
            pipe = pickle.load(infile)
    except :
        print("The chosen file does not exist or does not contain the pipeline list \n : %s",filename)
        pipe =None
    return pipe

def run_pipeline(pipe,logf):
    
    with open(logf,'a') as logfile:
        for i in range(len(pipe)):
            p_el = pipe[i]
            #p_el.status=0
            logfile.write("Run " + p_el.alias +'\n')
            print("\n\n> Run " + p_el.alias )
            print ('status = %d'% p_el.status)
            p_el.check_condition()
            print ('\t> Step status = %d  \n\t\t[-1 = Successfully executed | 0 = to be run | 1+ = number of consecutive failed attempt]'% p_el.status)
            if p_el.status >=0:
                p_el.run()
                if p_el.status<0:
                    logfile.write("Success for "+p_el.alias+'\n')
                    if p_el.condition != None:
                        p_el.update_condition()
                    print("\t> Success for "+p_el.alias+'\n')
                    save_pipeline(pipe,filename)
                else:
                    print p_el.status
                    logfile.write(p_el.err_msg +'\n')
                    i = len(pipe)+1                    
                    sys.stderr.write(p_el.err_msg +'\n')
                    raise
            else:
                logfile.write("Skip this step "+ p_el.alias +" because already executed successfully \n")
                print("\t> Skip this step "+ p_el.alias +" because already executed successfully \n")
    return None
    
    
    
if __name__=='__main__':
    import sys
    import pickle
    from datetime import datetime
    import pipeline_element
    
    try:
        filename = sys.argv[1]
        logfile  = sys.argv[2]        
    except:
        filename = None
    print "Executing pipeline:"
    #print filename
    #print logfile
    
    if filename:
        pipe = load_pipeline(filename)      
        run_pipeline(pipe,logfile)
        save_pipeline(pipe,filename)
    else:
        print "Error: incorrect syntax '%s'" % filename
        pass