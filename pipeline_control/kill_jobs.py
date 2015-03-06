#!/usr/bin/env python

if __name__=='__main__':
    import subprocess
    import os
    import sys
    import datetime
    import pickle
    if not os.environ.get('DRMAA_LIBRARY_PATH'):
        os.environ['DRMAA_LIBRARY_PATH'] = "/usr/share/univage/lib/lx-amd64/libdrmaa.so.1.0"
        import drmaa
    else:
        import drmaa
    import qsubclass

    #Idea is to load a file in argv[1] and run it:
    filename = sys.argv[1]
    try:
        with open(filename,'rb') as infile:
            qlist = pickle.load(infile)
            for ql in qlist:
                if int(ql.jobid)>0:
                    jobid = int(ql.jobid)
                    #print "qdel %d "%jobid
                    subprocess.call("qdel %d"%jobid,shell=True)
                    print "Killing job: %d"%jobid
            print "Done!"
    except:
        print "Error killing the pipeline-"