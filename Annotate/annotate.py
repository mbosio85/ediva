#!/usr/bin/env python

from datetime import datetime
import annotate_support_functions
import subprocess
import difflib, sys
import threading, Queue
import time
import os.path
start = time.time()
######################################
#
#	Task => Annotate variants 
#	infile => VCF file with complete sample wise genotype information
#	outfile => text file (csv) with complete annotation with sample wise genotype information
#	Extra outfile => text file (csv) without genic annotation with sample wise information for variants that are not bi-allelic (e.g tri-allelic) 
#
#######################################

#if needed print usage:

##############################################################################################
## SETTINGS
##############################################################################################

## variables
help_            = 0
infile          = ""
geneDef         = "refGene" ## gene Definition
sep             = "," ## separator for annotation outfile currently comma (,) is default
type_var        = "all" ## type of variants to annotate from input vcf file
gtMode          = "complete" ## type of variants to annotate from input vcf file
onlygenic       = False ## variable for only genic annotation 
forceDel        = False ## varibale for force deleting the output annotation file (if exists)
qlookup         = "NA" ## varibale for enabling quick lookup mode of the program
templocation    = "INPATH" ## scratch place for creating the temp files while annotating

variants        = dict() ## hash to hold input variants from VCF
not_biallelic_variants = dict() ## hash to hold input variants from VCF where the site is not bi-allelic and annovar annotation is absent
thrds           = list() ## list to hold threads

ediva           = dict() ## hash to hold input INDELs with ediva annotation  SHARED
Annovar         = dict() ## hash to hold input INDELs with ANNOVAR annotation SHARED
samples         = dict()  ## hash to hold sample information SHARED
edivaStr        = dict()  ## hash to hold simple tandem repeat data from ediva public omics database SHARED
headers         = list()

fileSuffix      = str(datetime.now().time())
fileSuffix      = fileSuffix.replace(':','_')
fileSuffix      = fileSuffix.replace('.','_')

## ANNOVAR settings
ANNOVAR         = "/users/GD/tools/eDiVaCommandLine/lib/Annovar"
TABIX           = "PATH=$PATH:/users/GD/tools/tabix/"
##############################################################################################
## INPUT PARSE - ANNOVAR CONFIGURATION CHECK - OUTPUT CREATION
##############################################################################################

parser_ = {"geneDef": geneDef, "type": type_var,"infile":infile,
           "onlygenic":onlygenic,"qlookup":qlookup,
           "forcedel":forceDel,"gtmode":gtMode,
           "templocation":templocation,"help":help_}# Arguments of the input parser

parser_ = annotate_support_functions.input_parse(parser_)
help_           = parser_["help"]
infile          = os.path.abspath(parser_["infile"])
geneDef         = parser_["geneDef"]
type_var        = parser_["type"]
gtMode          = parser_["gtmode"] 
onlygenic       = parser_["onlygenic"]
forceDel        = parser_["forcedel"]
qlookup         = parser_["qlookup"]
templocation    = parser_["templocation"]

try:
    (vcftoAnn,genicAnn) = annotate_support_functions.annovar_check(ANNOVAR)
except IOError:
    sys.exit(1)

try:
    (outFile,SortedoutFile,outFileIns,templocation) = annotate_support_functions.out_file_generate(infile,qlookup,templocation,forceDel,fileSuffix)

except IOError:
    sys.exit(1)

##############################################################################################
## MAIN starts
##############################################################################################


## start processing input VCF file
print "MESSAGE :: Processing input VCF file - $input "
## start browsing over the VCF file
try:
    (samples,variants,not_biallelic_variants,headers) = annotate_support_functions.vcf_processing(infile,qlookup,gtMode,type_var)
except IOError:
    sys.exit(1)
## Initialization completed
print "MESSAGE :: Finished processing input VCF file - %s "%(infile);

## prepare missing data handler for db annotation
(missandb,missandb_coordinate,missanndbindel) = annotate_support_functions.preparemissdb(sep) #<<<<<<<<<<< missdb etc values!


##?This should be set to multithread once we checked it works properly [at least for annotation inside the pipeline]
if (qlookup == "NA"):
    out_queue = Queue.Queue()
    if (onlygenic):
        ## start a sigle thread for annovar genic annotation
        print "MESSAGE :: Annotation starting "
        Annovar = annotate_support_functions.AnnovarAnnotation(infile,templocation,fileSuffix,geneDef,ANNOVAR,Annovar,TABIX) ## spawn a thread for Annovar annotation
    else:
        ## start threading and annotating
        print "MESSAGE :: Annotation starting "
        #ediva = annotate_support_functions.edivaAnnotation(variants,not_biallelic_variants,sep,missandb,missandb_coordinate,missanndbindel)## spawn a thread for ediva annotation
        #Annovar = annotate_support_functions.AnnovarAnnotation(infile,templocation,fileSuffix,geneDef,ANNOVAR,Annovar) ## spawn a thread for Annovar annotation
        #edivaStr = annotate_support_functions.edivaPublicOmics() ## spawn a thread for ediva public omics
        thread_ediva = threading.Thread(out_queue.put(annotate_support_functions.edivaAnnotation(variants,not_biallelic_variants,sep,missandb,missandb_coordinate,missanndbindel)))
        thread_annovar= threading.Thread(out_queue.put(annotate_support_functions.AnnovarAnnotation(infile,templocation,fileSuffix,geneDef,ANNOVAR,Annovar,TABIX))) ## spawn a thread for Annovar annotation
        #thread_pub = threading.Thread(out_queue.put(annotate_support_functions.edivaPublicOmics())) ## spawn a thread for ediva public omics
        thread_ediva.start()
        thread_annovar.start()
        #thread_pub.start()
        # now get the queue values and join threads
        thread_ediva.join()
        thread_annovar.join()
        #thread_pub.join()
        ediva    = out_queue.get()
        Annovar  = out_queue.get()
        #edivaStr = out_queue.get()
        print 'threading done'
else:
    out_queue = Queue.Queue()
    thread_ediva = threading.Thread(out_queue.put(annotate_support_functions.edivaAnnotation(variants,not_biallelic_variants,sep,missandb,missandb_coordinate,missanndbindel)))
    thread_ediva.start()
    #if os.path.exists(qlookup):
    #    thread_pub = threading.Thread(out_queue.put(annotate_support_functions.edivaPublicOmics())) ## spawn a thread for ediva public omic
    #    thread_pub.start()
    # now get the queue values and join threads
    thread_ediva.join()
    #if os.path.exists(qlookup):
        #thread_pub.join()
    ediva    = out_queue.get()
    #if os.path.exists(qlookup):
    #    edivaStr = out_queue.get()
    print 'threading done'


## join spawned threads
## write annotation to file or ender output
if qlookup == "NA":
    ## write annotaton in file 
    print "MESSAGE :: Writing annotation to output file %s" % (outFile)
    ## open file handler
    with open(outFile,'w+') as ANN:
        ## write header to output file
        headerOutputFile = annotate_support_functions.getHeader(onlygenic,geneDef,headers)
        ANN.write(headerOutputFile+'\n')
        ## write data lines to main output file
        for key, value in variants.items(): 
            (chr_col,position,ref,alt,aftoprint,qual,filter_) = value.split(';')
            annovarValueToMatch = ';'.join((chr_col,position,ref,alt))
            edivaannotationtoprint = ediva.get(key,"NA")
            #print edivaannotationtoprint
            annovarannotationtoprint = Annovar.get(annovarValueToMatch,"NA,"*3+"NA")
            samplewiseinfortoprint = samples.get(key,"NA")
            #edivapublicanntoprint = edivaStr.get(';'.join((chr_col,position)),"NA,NA")
            # write annotation to file

            write_str=(chr_col+sep+position+sep+ref+sep+alt+sep+
                       qual+sep+filter_+sep+
                       aftoprint+sep+                   
                       annovarannotationtoprint+sep+edivaannotationtoprint+sep+samplewiseinfortoprint)
                      #edivapublicanntoprint+sep+samplewiseinfortoprint)

            write_str.replace('\n','')

            ANN.write(write_str+'\n')
    print "MESSAGE :: Writing annotation to output file %s" % (outFileIns)
    with open(outFileIns,'w+') as ANNINS:
        ## write header for inconsistent file
        headerOutputFile = annotate_support_functions.getHeaderIns(headers)
        ANNINS.write(headerOutputFile+'\n')
        ## write data lines to main output file
        for key, value in not_biallelic_variants.items(): 
            edivaannotationtoprint,annovarannotationtoprint,samplewiseinfortoprint = ("NA","NA","NA")
            edivapublicanntoprint = "NA,NA"
            (chr_col,position,ref,alt,aftoprint,qual,filter_) = value.split(';')
            edivaannotationtoprint = ediva.get(key,"NA")
            samplewiseinfortoprint = samples.get(key,"NA")
            #edivapublicanntoprint = edivaStr.get(';'.join((chr_col,position)),"NA,NA")
            ## write annotation to fileprint
            
            write_str=(chr_col+sep+position+sep+ref+sep+alt+sep+
                       qual+sep+filter_+sep+
                       aftoprint+sep+
                       annovarannotationtoprint+sep+edivaannotationtoprint+sep+samplewiseinfortoprint)
                       #edivapublicanntoprint+sep+samplewiseinfortoprint)
            write_str.replace('\n','')
            ANNINS.write(write_str+'\n')

            
    ## sort the file
    srtCmm = "sort -k1,1 -n -k2,2 --field-separator=, %s > %s " %(outFile,SortedoutFile)
    subprocess.call(srtCmm,shell=True)
    ## writing completed
    print "MESSAGE :: Writing annotation completed "
    print "MESSAGE :: Your annotated file is %s " %(outFile)
    print "MESSAGE :: Your sorted annotated file is %s "%(SortedoutFile)
    print "MESSAGE :: Reported non bi-allelic sites are in %s " %(outFileIns)
    ## Finalize everything
    print "MESSAGE :: Finalizing annotation process "
    annotate_support_functions.finalize(templocation,fileSuffix)
    print "MESSAGE :: Finalization completed "
    print "MESSAGE :: Templocation %s cleared"%(templocation)
else:
    if os.path.isfile(qlookup):
        ## render annotation to output file
        with open(outFile,'a') as ANN:
            headerOutputFile = annotate_support_functions.getHeaderQlookup(headers)
            ANN.write(headerOutputFile+'\n')
            ## get header
            for key, value in variants.items():
                print value
                (chr_col,position,ref,alt,dummy) = value.split(';')
                annovarValueToMatch = ';'.join((chr_col,position,ref,alt))
                edivaannotationtoprint = ediva.get(key,"NA")
                #print edivaannotationtoprint
                edivapublicanntoprint = edivaStr.get(';'.join((chr_col,position)),"NA,NA")
                write_str=(chr_col+sep+position+sep+ref+sep+alt+sep+edivaannotationtoprint+sep+
                      edivapublicanntoprint)
                write_str.replace('\n','')
                ANN.write(write_str+'\n')
        ## sort the file
        srtCmm = "sort -k1,1 -n -k2,2 --field-separator=, %s > %s " %(outFile,SortedoutFile)
        subprocess.call(srtCmm,shell=True)
        ## writing completed
        print "MESSAGE :: Writing annotation completed"
        print "MESSAGE :: Your annotated file is %s " %outFile
        print "MESSAGE :: Your sorted annotated file is %s " %SortedoutFile
    else:
        for key, value in variants.items():
            (chr_col,position,ref,alt,aftoprint) = value.split(';')
            annovarValueToMatch = ';'.join((chr_col,position,ref,alt))
            edivaannotationtoprint = ediva.get(key,"NA")
            #print edivaannotationtoprint
            edivapublicanntoprint = annotate_support_functions.edivaPublicOmics_search(chr_col,position)
            #  edivaStr.get(';'.join((chr_col,position)),"NA,NA")
            edivavals =  edivaannotationtoprint.split(',')
            edivapublicvals = edivapublicanntoprint.split(',')
            ## render to command line output
            print "chromosome: %s "                     % chr_col
            print "position: %s "                       %position
            print "Reference: %s "                      %ref
            print "Alteration: %s "                     %alt
            print "dbSNP identifier: %s"                % edivavals[0]
            print "EVS european frequency: %s "         % edivavals[1]
            print "EVS african frequency: %s "          % edivavals[2]
            print "EVS total frequency: %s "            % edivavals[3]
            print "1000genomes european frequency: %s " % edivavals[4]
            print "1000genomes african frequency: %s "  % edivavals[5]
            print "1000genomes american frequency: %s " % edivavals[6]
            print "1000genomes asian frequency: %s "    % edivavals[7]
            print "1000genomes total frequency: %s "    % edivavals[8]
            print "Segment duplication: %s "            % edivavals[9]
            print "Placental mammal phyloP: %s "        % edivavals[10]
            print "Primates phyloP: %s "                % edivavals[11]
            print "Vertebrates phyloP: %s "             % edivavals[12]
            print "Placental mammal phastcons: %s"      % edivavals[13]
            print "Primates phastcons: %s "             % edivavals[14]
            print "Vertebrates phastcons: %s "          % edivavals[15]
            print "Gerp score1: %s "			        % edivavals[16]
            print "Gerp score2: %s "                    % edivavals[17]
            print "Sift: %s "			                % edivavals[18]
            print "polyphen2: %s "                      % edivavals[19]
            print "Mutationassessor: %s "		        % edivavals[20]
            print "Condel: %s "                         % edivavals[21]
            print "Cadd score1: %s "		            % edivavals[22]
            print "Cadd score2: %s "                    % edivavals[23]
            print "Simple tandem repeat region: %s "    % edivapublicvals[0]
            print "Simple tandem repeat length: %s "    % edivapublicvals[1]

end = time.time()
py_time = end-start
