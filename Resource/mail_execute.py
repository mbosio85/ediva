

if __name__=='__main__':

    import subprocess    
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-csv"  , dest='csv',required=True, type=str, help="csv file with mail ")
    parser.add_argument("-text"  , dest='t',required=True, type=str, help="Text to add")   
    #parser.add_argument("-infile"  , dest='infile',required=True, type=str, help="infile")   
    args = parser.parse_args()
    
   #/home/rrahman/soft/python-mailer/rank.html
    if os.path.isfile(args.csv):
        mailCmd = 'python /home/rrahman/soft/python-mailer/pymailer.py -s '+args.t+ ' '+  args.csv +' Ranking'
        print mailCmd
	subprocess.call(mailCmd,shell=True)
