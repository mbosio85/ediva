import pprint
import argparse
import readline
import os

pp = pprint.PrettyPrinter(indent = 5)

parser = argparse.ArgumentParser(description = 'Create a file containing all necessary information to run the eDiVa prioritization.')

parser.add_argument('--outfile',  type=str, dest='config',  required=True, help='The output file will contain all sample names, the affection status and the file locations of the vcf and bam file.')
parser.add_argument('--sampleinfo', type=str, dest='sampleinfo', required=False, action='append', help="""Instead of going through the wizzard, you can also provide all information in one line.
Please follow the pattern --sampleinfo ID:affected:vcf:bam
affected is either 0 (unaffected) or 1 (affected)
vcf is the full vcf file path
bam is the full bam file path [optional]

This option can be called multiple times.
If GATK multi-sample calls are supplied, just enter the same vcf file name for each sample and toggle the --multisample option creating the pipeline.

""")

args = parser.parse_args()
token = True
line_list = list()

# write header
with open(args.config,'w+') as config_file:
    config_file.write('\t'.join(['ID', 'status', 'vcf', 'bam']))
    config_file.write('\n')
    
    if args.sampleinfo:
        for line in args.sampleinfo:
            splitline = line.split(':')
            name = splitline[0]
            affect = splitline[1]
            try:
                vcf = splitline[2]
                vcf = os.path.expanduser(vcf)
            except Exception,e:
                print("There was an error reading the folder.")
                print(e)
            try:
                bam = splitline[3]
                try:
                    bam = os.path.expanduser(bam)
                except Exception,e:
                    print("There was an error reading the folder.")
                    print(e)
            except:
                bam = str()
            
            args.config.write('\t'.join([name, affect, vcf, bam]))
            args.config.write('\n')
    
    
    else:       
        ok_john_go_on = True
        while ok_john_go_on:#If list : then review
            # for each element ask questions with default value on enter
            if len(line_list)>0:
                print"-------------------------------"
                for l in range(len(line_list)):
                    print ("\nREVIEWING SAMPLE:\n")
                    params = line_list[l].split('\t')
                    print params
                    readline.parse_and_bind("tab: complete")
                    readline.set_completer_delims(' \t\n;')
                    name   = raw_input("\tReviewing sample ID: %s  Press enter to maintain it\n\t: "%params[0])or params[0]
                    token = 1
                    while token:
                        affect =  (raw_input("\tPlease enter if the sample is affected or not. [0 - unaffected, 1 - affected] \n\t" +  
                                             "Press enter for %s \n\t:"%params[1]) or params[1]  ) 
                        if affect in  ['0','1']:
                            token = 0
                    readline.parse_and_bind("tab: complete")
                    readline.set_completer_delims(' \t\n;')
                    print("\tPlease enter the location of the respective vcf file. \n")
                    vcf = raw_input("\tPress enter to keep %s\n\t:"%params[2]) or params[2]
                    
                    readline.parse_and_bind("tab: complete")
                    readline.set_completer_delims(' \t\n;')
                    print("\tPlease enter the location of the respective bam file. \n")
                    bam    = raw_input("\tPress enter to keep %s\n\t:"%params[3]) or params[3]
                    line_list[l]  = '\t'.join([name,affect,vcf,bam])
                    
            
            print"-------------------------------"
            while True:
                print ("\nNEW SAMPLE:\n")
                print("If you supply GATK multi-sample calls, then just mention the same file name, if you're asked for a vcf file.")
                print("Remember. If you want to stop entering data, just leave a field empty and hit enter.\n")
                readline.parse_and_bind("tab: complete")
                readline.set_completer_delims(' \t\n;')
                name   = raw_input("\tPlease enter the sample ID. The format should be same as the vcf file, you will provide later. \n\t: ")
                if name == '':
                    print("\tField empty. Quitting...")
                    break
                token = 1
                while token:
                    affect =  raw_input("\tPlease enter if the sample is affected or not. [0 - unaffected, 1 - affected] \n\t: ")
                    if affect == '':
                        print("\tField empty. Quitting...")
                        token = 0
                        break
                    elif affect in  ['0','1']:
                        token = 0
                if affect == '':
                        print("\tField empty. Quitting...")
                        break  
                readline.parse_and_bind("tab: complete")
                readline.set_completer_delims(' \t\n;')
                vcf    = raw_input("\tPlease enter the location of the respective vcf file. \n\t: ")
                if vcf == '':
                    print("Field empty. Quitting...")
                    break
                
                readline.parse_and_bind("tab: complete")
                readline.set_completer_delims(' \t\n;')
                bam    = raw_input("\tPlease enter the location of the respective bam file. \n\t: ")
                if bam == '':
                    print("\tField empty. Quitting...")
                    break
                line_list.append('\t'.join([name, affect, vcf, bam]))
                #config_file.write('\t'.join([name, affect, vcf, bam]))
                #config_file.write('\n')
            
            while True:
                print"-------------------------------"
                print ("\nDATA CONFIRMATION:\n")
                tmp = raw_input("\tDo you want to review your data?\n\t:'y' or 'n', Hit enter for 'n'\n\t:")or 'n'
                if tmp == 'y':
                    break
                elif tmp =='n':
                    ok_john_go_on=False
                    break
                
        for line in line_list:
            config_file.write(line)
            config_file.write('\n')
print"-------------------------------"            
print"\nRESUME of the CONFIG FILE \n"
with open(args.config,'r+') as config_file:
    for line in config_file:
        print line
print"-------------------------------"          