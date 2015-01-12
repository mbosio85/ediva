#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description = 'Create a file that will be run ')
parser.add_argument('--infile',  type=str, dest='infile',  required=True, help='This is the input file')
parser.add_argument('--outfile',  type=str, dest='outfile',  required=True, help='This is the output file on which we will write')
parser.add_argument('--value', type=int, dest='value', required=False, default = 5, help='name of the script you want to start later on')
args = parser.parse_args()

         
#start reading the file and parsing the output to te output files

with open(args.infile,'r') as infile_read, open(args.outfile,'w') as out_write, open(args.outfile[0:-4] + ".lowQual.vcf",'w') as bad_write:
    for line in infile_read:
        flag = 0  # IF 0 -> line goes to bad file | If 1 -> goes to outfile
        if line[0] == '#':
            #header line, write it to both outputs
            bad_write.write(line)
            out_write.write(line)
        else:
            tabline = line.split('\t')
            #Now tabline [0...8] are useless
            #We analyze tabline[9+] looking for quality measures
            for i in range(9,len(tabline)):
                ## data format is= GT:AD:DP:GQ:PL
                datasplit = tabline[i].split(':')
                ## check only when zygosity is not 0/0 or ./.
                if (datasplit[0] != "0/0") and (datasplit[0] != "./."):
                    AD =  datasplit[1].split(',')
                    flag=1
                    break
            # Then here we decide where to put the vcf file line     
            if flag:
                out_write.write(line)
            else:
                bad_write.write(line)
                print('.')
         #!/usr/bin/env python

