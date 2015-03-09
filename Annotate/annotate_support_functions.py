
#!/usr/bin/env python
import os
import sys
import subprocess
#import mysql.connector
import MySQLdb
import re
import ntpath
import struct
import hashlib
import argparse
from Bio import bgzf



######################################
#
#	Task => Annotate variants 
#	infile => VCF file with complete sample wise genotype information
#	outfile => text file (csv) with complete annotation with sample wise genotype information
#	Extra outfile => text file (csv) without genic annotation with sample wise information for variants that are not bi-allelic (e.g tri-allelic) 
#
#######################################
i_str = 'VCF file containing the variants to annoate'
q_str = ('Option : Command line variant argument. Variant format => chr\:position\:reference\:alteration \n ' + 
		'\t\t\t Option : File containing variants (one variant per each line). '+
		'Variant format => chr\:position\:reference\:alteration \n\n')
t_str =  ('Temporary scratch location for temp files (all temp files with the current session will be removed '+
		' at the end of execution) \n\t\t\t default: input VCF location ')
g_str = 'Gene deifnition you want to select for genic annotation (ensGene,refGene,knownGene,all) \n\t\t\t default: refGene'
v_str = 'Type of variants to annotate from input VCF file (SNP,INDEL,all) \n\t\t\t default: all'
s_str = ('complete: reports all genotypes from the input VCF file\n\t\t\t compact: reports only '+
		 'heterozygous and homozygous alteration genotypes from the input VCF file \n\t\t\t '+
		 'none: exclude sample wise genotype information in annotation \n\t\t\t default: compact')
o_str = 'If set, then only genic annotation will be performed'
f_str = 'If set, then it will over-write existing output annotation file with the same name '
## subroutine for usage of the tool @@
def usage():
    print ("\n usage:\n" +
    "--input,-i \t\t" 			+ i_str + " \n\n" +
    "--quicklookup,-q \t "		+ q_str +  "\n\n" +
    "--tempDir,-t \t\t "		+ t_str + " \n\n" +
    "--geneDef,-g \t\t " 		+ g_str + " \n\n" +
    "--variantType,-v \t "		+v_str	+ " \n\n" +
    "--sampleGenotypeMode,-s  "	+ s_str	+ " \n\n" +
    "--onlyGenicAnnotation,-o "	+o_str	+ " \n\n" +
    "--forceNewFileCreate,-f  "	+f_str	+ " \n\n" +
    "--help,\t\t show help \n\n" )
    raise IOError
##############################################################################################
## COMMAND LINE OPTIONS done- to be tested
##############################################################################################
def input_parse(defaults):
    parser_ = defaults
        ## grab command line options
    parser = argparse.ArgumentParser(description = 'Setup ')
    parser.add_argument('-i','--input'  ,            type=str, dest ="infile",       required=False, 					 help= i_str)
    parser.add_argument('-t','--tempDir',            type=str, dest ="templocation", required=False, 					 help= t_str)
    parser.add_argument('-v','--variantType'  ,      type=str, dest ="type",         required=False, 					 help= v_str)
    parser.add_argument('-o','--onlyGenicAnnotation',		   dest ="onlygenic",    required=False, action='store_true',help= o_str)
    parser.add_argument('-f','--forceNewFileCreate', 	  	   dest ="forcedel",     required=False, action='store_true',help= f_str)
    parser.add_argument('-q','--quicklookup'  ,      type=str, dest ="qlookup",      required=False, 					 help= q_str)
    parser.add_argument('-s','--sampleGenotypeMode', type=str, dest ="gtmode",       required=False, 					 help= s_str)
    parser.add_argument('-g','--geneDef',            type=str, dest ="geneDef",      required=False, 					 help= g_str)  
    args = parser.parse_args()


    if args.infile:
        parser_["infile"] = args.infile
    if args.geneDef:
        parser_["geneDef"]=args.geneDef
    if args.type:
        parser_["type"] = args.type
    if args.gtmode:
        parser_["gtmode"] =args.gtmode
    if args.onlygenic:
        parser_["onlygenic"]=args.onlygenic
    if args.forcedel :
        parser_["forcedel"] =  args.forcedel        
    if args.qlookup:
        parser_["qlookup"] = args.qlookup
    if args.templocation:
        parser_["templocation"] = args.templocation
        

    ### check mandatory command line parameters and take necessary actions
    if  len(parser_["infile"] )==0 and parser_["qlookup"] == "NA":
        usage()
        sys.exit()

    ### check on both input file and quick lookup mode
    if len(parser_["infile"])>0 and parser_["qlookup"] != "NA":
        print "\n ERROR :: Either run the program with a VCF file or by using quick lookup command line arguments \n"
        sys.exit()
    #
    ### final check of gene definition
    if parser_["geneDef"] != "refGene" and parser_["geneDef"] != "ensGene" and parser_["geneDef"]!= "knownGene" and parser_["geneDef"] != "all":
        print "\n WARNING :: Not a valid gene definition. Please select a correct gene definition and if you are not sure then use the default settings !\n"
        usage()
        sys.exit()
    ### final check of variant type
    if  parser_["type"] != "SNP" and  parser_["type"]!= "INDEL" and  parser_["type"]!= "CNV" and  parser_["type"]!= "all":
        print "\nWARNING :: Not a valid variant type. Please select a correct variant type and if you are not sure then use the default settings !\n"
        usage()
        sys.exit()
    ### final check of sample genotype mode type
    if parser_["gtmode"]!= "complete" and parser_["gtmode"] != "compact" and parser_["gtmode"] != "none":
        print "\nWARNING :: Not a valid sample genptype mode value. Please select a correct value and if you are not sure then use the default settings !\n";
        usage()
        sys.exit()
    ### final check of temp location
    if parser_["templocation"] != "INPATH" and not(os.path.isdir(parser_["templocation"] )):
        print "\nWARNING :: Not a valid temporary location or the program does not have access to this location. If you are not sure then use the default settings !\n"
        usage()
        sys.exit()
    ### final check on quick lookup mode parameter
    if parser_["qlookup"] != "NA":
        if os.path.exists(parser_["qlookup"]):
            pass
            ## nothing to do ! its a file !! varinat format inside the file is checked later on !!
        else:
            if not(':' in parser_["qlookup"]):
                print "\nERROR :: Not a valid format. Correct format is chr:position:reference:alteration \n"
                sys.exit()
            
    
    return parser_

##############################################################################################
## CONFIGURATION VARIFY "done"
##############################################################################################
def annovar_check(ANNOVAR):
    ### check of annovar settings
    vcftoAnn = ANNOVAR+"/convert2annovar.pl"
    genicAnn = ANNOVAR+"/ediva_summarize_annovar.pl"
    if not(os.path.isdir(ANNOVAR)):    #{
        print "\nERROR :: Annovar library location is not found. Please set the Annovar library path correctly inside the program \n"
        raise IOError
    if not(os.path.exists(vcftoAnn)):
        print "\nERROR :: Program for converting VCF to Annovar format is not found in the Annovar library location \n"
        raise IOError

    if not(os.path.exists(genicAnn)):
        print "\nERROR :: Program for converting VCF to Annovar format is not found in the Annovar library location \n";
        raise IOError

    return (vcftoAnn,genicAnn)

##############################################################################################
## OUTPUT FILE(s) "done"
##############################################################################################
def out_file_generate(infile,qlookup,templocation,forceDel,tempfile):
    outFile         = ''
    sortedOutFile   = ''
    outFileIns      = ''
    
    if qlookup == "NA":
        filename = os.path.basename(infile)[:-4]
        pathname = os.path.dirname(infile)
        if len(pathname)>1:
            pathname = pathname +'/'
            if templocation == "INPATH":
                templocation = pathname
        elif templocation == "INPATH":
                templocation = "."
        
	
        outFile  = pathname+ filename+".annotated"
        sortedOutFile = pathname+ filename +".sorted.annotated"
        outFileIns = pathname+filename + ".inconsistent.annotated"       
        ## check on output file existence
        if os.path.exists(outFile) or os.path.exists(sortedOutFile):
            ## check for new file creation flag
            if forceDel:
                ## delete the files if exist
                try:
                    os.remove(outFile)
                except:
                    pass
                try:
                    os.remove(sortedOutFile)
                except:
                    pass
                
                print "MESSAGE :: Target output file(s) already exists. Removing them now \n"
            else:
                print "\n  WARNING :: Target output file(s) already exists. Either rename them, remove them or set the --forceNewFileCreate variable \n";
                raise IOError
    else:
        if (os.path.exists(qlookup)):
            filename = ntpath.basename(qlookup)
	    filename = '.'.join(filename.split('.')[:-1])
            pathname = ntpath.dirname(qlookup)
            #filename = ntpath.basename(qlookup)[:-4]
            #pathname = ntpath.dirname(qlookup)
            if len(pathname)>1:
                pathname = pathname +'/'
            outFile  = pathname+ filename+".annotated"
            sortedOutFile = pathname + filename +".sorted.annotated"
            outFileIns = pathname+ filename + ".inconsistent.annotated"
	    print outFile
            if os.path.exists(outFile) or os.path.exists(sortedOutFile):
            ## check for new file creation flag
                if forceDel:
                    ## delete the files if exist
                    try:
                        os.remove(outFile)
                    except:
                        pass
                    try:
                        os.remove(sortedOutFile)
                    except:
                        pass
                    
                        print "MESSAGE :: Target output file(s) already exists. Removing them now \n"
                else:
                    print "\nWARNING :: Target output file(s) already exists. Either rename them, remove them or set the --forceNewFileCreate variable \n";
                    raise IOError
    print("MESSAGE :: Your annotated file is: %s"%(outFile))
    print("MESSAGE :: Your sorted annotated file is: %s" %(sortedOutFile))
    print("MESSAGE :: Reported non bi-allelic sites are in: %s" %(outFileIns))
    return(outFile,sortedOutFile,outFileIns,templocation) ## files to write annotation

##############################################################################################
## SUBROUTINES "done"
##############################################################################################



## subroutine for unknown command line arguments @@
def unknownArguments():
    print "ERROR :: Unknown command line arguments were found \n Please see the following usage";
    usage()
    sys.exit()
    raise IOError


## subroutine for finalizing annotation process @@
def finalize(templocation,fileSuffix):
    ## clear the tmp directory for this session
    clearCmm = "rm -r " + templocation + "/*"+fileSuffix +"*"
    #print(clearCmm)
    subprocess.call(clearCmm,shell = True)
    return None

## sub for preparing missing db annotation @@
def preparemissdb(sep):
    missandb            = 'NA'
    missandb_coordinate = '0'
    missanndbindel       = 'NA'
    
    missandb = missandb + (sep+"0")*9
    missandb = missandb + (sep+"NA")*14
    missandb_coordinate = missandb_coordinate + (sep+"NA")*8
    missanndbindel = missanndbindel + (sep+"0")*8

    return(missandb,missandb_coordinate,missanndbindel)
#@@ should be done: check it a couple of times with examples
## sub for replacing commas inside double qoutes for annovar genic annotation lines
def replaceCommainQoute(in_str):
    split_str = in_str.split('"');
    tmplist =[]
    token = 1
    #if in_str.startswith('"'):
    #    token=0
    for s in split_str[token:-1:2]:
        tmplist.append(s.replace(',',';'))
    split_str[token:-1:2] = tmplist
    return('"'.join(split_str))
#@@Done
## subroutine for ediva public omics data fetch
def edivaPublicOmics():
    edivaStr = dict()
    ## DB parameters
    username    = "edivacrg"
    database    = "eDiVa_public_omics"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "FD5KrT3q"
     
    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base
    
    cur = db.cursor()
    sql = "select chr,pos,lengthofrepeat,region from eDiVa_public_omics.Table_simpleRepeat;"
        #sql = "select chr,pos,lengthofrepeat,copyNum,region from ediva_public_omics.Table_simpleRepeat;"

    cur.execute(sql)
    print('\t ediva public omics start')
    # Fill the edivaStr dictionary with the gathered values
    # Key is obtained by res[0];res[1]] (chr ; pos)
    # Value is res[4],res[2]  (region,lengthofrepeat)
    
    edivaStr= [(str(row[0])+';'+str(row[1]) , str(row[3])+','+str(row[2])   ) for row in cur]
    edivaStr = dict(edivaStr)
    cur.close()
    db.close()
    return edivaStr


def edivaPublicOmics_search(chr,pos):
    out = 'NA,NA'
    username    = "edivacrg"
    database    = "eDiVa_public_omics"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "FD5KrT3q"

    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base

    cur = db.cursor()
    sql = "select lengthofrepeat,region from eDiVa_public_omics.Table_simpleRepeat where chr =%s and pos = %s limit 1"%(chr,pos)
        #sql = "select chr,pos,lengthofrepeat,copyNum,region from ediva_public_omics.Table_simpleRepeat;"
    cur.execute(sql)

    for row in cur:
        out = str(row[1])+','+str(row[0])
    cur.close()
    db.close()
    #print out

    return out
#@@ Done
def process_db_entry(res,res2,ediva,k,sep,missanndb_coordinate,missanndbindel):
    if res is None:
        res = list()
    if res2 is None:
        res2 = list()
    if len(res) > 1 and len(res2) > 1: ## both returned database rows
        # load ediva hash from database
        for reselement in res:
            if ediva.get(k,0):
                ediva[k] = ediva[k] + (sep+str(reselement))
            else:
               ediva[k] = str(reselement)     
        for reselement2 in res2 :
            if ediva.get(k,False):
                ediva[k] = ediva[k] + (sep+str(reselement2))
            else:
               ediva[k] = str(reselement2)
        ## add NAs for damage potential scores and cadd scores for indels
        ediva[k] +=  ((sep+ "NA")*6)
        #ediva[k] = ediva[k] + ((sep+ "NA")*6)
    elif len(res) > 1 and len(res2) < 1: ## only indel table returned database row
        # load ediva hash from database
        for reselement in res :
            if ediva.get(k,False):
                ediva[k] += (sep+str(reselement))
            else:
               ediva[k] = str(reselement)  
        # take care of missing positional values
        ediva[k] += (sep + missanndb_coordinate)
        ## add NAs for damage potential scores and cadd scores for indels
        ediva[k] += ((sep+ "NA")*6)
    elif len(res) < 1 and len(res2) > 1: ## only snp table returned database row
        # take care of missing positional values			
        if ediva.get(k,False):
            ## this should never happen !!
            ediva[k] += (sep + missanndbindel)
        else:
            ediva[k] =  missanndbindel
        for reselement in res2:    
            if ediva.get(k,False):
                ediva[k] += (sep+str(reselement))
            else:
                ediva[k] =  str(reselement)
        ## add NAs for damage potential scores for indels
        ediva[k] += ((sep+ "NA")*6)
    else:
        ## no entry in the database for both the queries
        if (ediva.get(k,False)):
            ## this should never happen !!
            ediva[k] += (sep+missanndbindel)
        else:
            ediva[k] =   (missanndbindel)
        # take care of missing positional values			
        ediva[k] +=(sep+missanndb_coordinate)
        ## add NAs for damage potential scores for indels
        ediva[k] +=((sep+ "NA")*6)
    return ediva
def ediva_reselement(ediva,k,res,sep):
    for reselement in res:
        if ediva.get(k,False):
            ediva[k] += sep+str(reselement)
        else:
            ediva[k]=str(reselement)
    return ediva

#@@ Check all input and output variables: evaluate if it's better to pass dictionaries
#@@ rather than variables.
#@@ There is much room for improvement in the factorization of code 
## subroutine for Annovar annotation
## subroutine for ediva annotation
def edivaAnnotation(variants,not_biallelic_variants,sep,missanndb,missanndb_coordinate,missanndbindel):
    ediva = dict()
    ## DB parameters
    username    = "edivacrg"
    database    = "eDiVa_annotation"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "FD5KrT3q"
    
    ## open DB connection
#    my $dbh = DBI->connect('dbi:mysql:'.$database.';host='.$dbhost.'',$username,$pass) or die "Connection Error!!\n";
    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base
    # extract db result
    print("\t Running ediva annotaton with %d variants : "%(len(variants.items())))
    for k, v in (variants.items()) :
        
        chr_col,pos,ref,alt = k.rstrip('\n').split(';')
        ## decide for variant type
        lenref = len(ref)
        lenalt = len(alt)
                
        if ((lenref + lenalt) > 2):## INDEL
            sql     = ""
            sql2    = ""
            stmt    = ""
            stmt2   = ""
            res     = list()
            res2    = list()
            ## prepare query 
            sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0') from Table_Chr%s_indel where indelid = '%s' limit 1;"""%(chr_col,k)
                        
            sql2 = """select ifnull(SDIndel,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA') 
			from Table_Chr%s where position = %s limit 1;"""%(chr_col,pos)
            sql = sql.replace('\n','')
            sql2 = sql2.replace('\n','')
            sql = sql.replace('\t','')
            sql2 = sql2.replace('\t','')
            ## prepare statement and query#######################
            cur = db.cursor()
            cur2 = db.cursor()
            cur.execute(sql)
            res = cur.fetchone()
            cur.close()
            cur2.execute(sql2)          
            #process query result            
            res2 = cur2.fetchone()
            cur2.close()
            ## prepare statement and query#######################
            ediva  = process_db_entry(res,res2,ediva,k,sep,missanndb_coordinate,missanndbindel)
        else:##SNP
            sql = ""
            stmt= ""
            res = list()
            #$sql = "select annotateSNPGermline('chr_col',$pos,'$ref','$alt');";					
            sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
			ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
			from Table_Chr%s where position = %s and Reference = '%s' and Alt = '%s' limit 1;"""%(chr_col,pos,ref,alt)
            sql = sql.replace('\n','')
            sql = sql.replace('\t','')
            cur = db.cursor()
            cur.execute(sql)


            #process query result
            res = cur.fetchone()
            #ediva  = process_db_entry_snp(res)
            if res is None:
                res = list()
            ## prepare statement and query#######################
            #process_db_entry(res,res2,ediva,k)
            cur.close()
            ## prepare statement and query#######################
            if (len(res) > 1):
                # load ediva hash from database
                ediva =  ediva_reselement(ediva,k,res,sep)
            else:
		#Less stringent search driven only by chr and position
                sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
				ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
				ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
				ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
				from Table_Chr%s where position = %s limit 1;"""%(chr_col,pos)
                sql = sql.replace('\n','')
                sql = sql.replace('\t','')
                ## prepare statement and query#######################
                cur = db.cursor()
                cur.execute(sql)
                
                #process query result
                res = cur.fetchone()
                if res is None:
                    res = list()
                ## prepare statement and query#######################
                #process_db_entry(res,res2,ediva,k)
                cur.close()
                ## prepare statement and query#######################
                
                if (len(res) > 1):
                    # load ediva hash from database
                    ediva =  ediva_reselement(ediva,k,res,sep)
                else:
                    ## handle missing database annotation entry
                    ediva[k]=missanndb    
    ## end of if-else for variant decision  
    # extract db result
    print("\t Running ediva annotation with %d non biallelic variants:"%(len(not_biallelic_variants.items())))
    for k, v in (not_biallelic_variants.items() ) :
        chr_col,pos,ref,alt = k.rstrip('\n').split(';')
        ## decide for variant type
        lenref = len(ref)
        lenalt = len(alt)
        if (lenref + lenalt) > 2:## INDEL
            sql     = ""
            sql2    = ""
            stmt    = ""
            stmt2   = ""
            res     = list()
            res2    = list()
            ## prepare query 
            #$sql = "select annotateINDEL('$k','chr_col',$pos);";
            sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0') from Table_Chr%s_indel where indelid = '%s' limit 1;"""%(chr_col,k)
            
            sql2 = """select ifnull(SDIndel,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA') 
			from Table_Chr%s where position = %s limit 1;"""%(chr_col,pos)
            sql = sql.replace('\n','')
            sql2 = sql.replace('\n','')
            sql = sql.replace('\t','')
            sql2 = sql.replace('\t','')
            
            ## prepare statement and query
            cur = db.cursor()
            cur2 = db.cursor()
            cur.execute(sql)
            res = cur.fetchone()
            if res is None:
                res = list()
            cur2.execute(sql2)
            
            #process query result
            
            res2 = cur2.fetchone()
            if res2 == None:
                res2 = list()
            ## prepare statement and query#######################
            #Now looks for the key in both databases: When present, the value is updated with the db info
            #When not present, the key is filled with the default "not found" elements like: misseldb etc.
            
            if (len(res) > 1 and len(res2) > 1): ## both returned database rows
                    # load ediva hash from database
                    ediva =  ediva_reselement(ediva,k,res,sep)
                    ediva =  ediva_reselement(ediva,k,res2,sep)
                    ## add NAs for damage potential scores for indels
                    ediva[k] = ediva[k] + ((sep+ "NA")*6)
            elif(len(res) > 1 and len(res2) < 1): ## only indel table returned database row
                ediva =  ediva_reselement(ediva,k,res,sep)
                # take care of missing positional values			
                ediva[k] += sep+missanndb_coordinate
                ## add NAs for damage potential scores for indels
                ediva[k] += ((sep+ "NA")*6)
                
            elif(len(res) < 1 and len(res2) > 1): ## only snp table returned database row
                # take care of missing positional values			
                if ediva.get(k,False):
                ## this should never happen !!
                    ediva[k] +=sep+missanndbindel
                else:
                    ediva[k] = missanndbindel
                
                ediva =  ediva_reselement(ediva,k,res2,sep)
                ## add NAs for damage potential scores for indels
                ediva[k] +=((sep+ "NA")*6)
            else:
                if ediva.get(k,False):
                ## this should never happen !!
                    ediva[k] +=sep+missanndbindel
                else:
                    ediva[k] = missanndbindel
                # take care of missing positional values			
                ediva[k] +=sep+missanndb_coordinate
                ## add NAs for damage potential scores for indels
                ediva[k] += ((sep+ "NA")*6)
            cur.close()
            cur2.close()
        else: ##SNP
			sql = ""
			stmt =""
			res = list()
		
			#$sql = "select annotateSNPGermline('chr_col',$pos,'$ref','$alt');";					
			sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
			ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
			from Table_Chr%s where position = %s and Reference = '%s' and Alt = '%s' limit 1;"""%(chr_col,pos,ref,alt)
			sql = sql.replace('\n','')
			sql = sql.replace('\t','')
		
			## prepare statement and query
			cur = db.cursor()
			cur.execute(sql)
			#process query result
			res = cur.fetchone()
			if res == None:
				res = list()
		
			if (len(res) > 1):
				# load ediva hash from database
				ediva =  ediva_reselement(ediva,k,res,sep)
			else:
				sql = """select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
				ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
				ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
				ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
				from Table_Chr%s where position = %s limit 1;"""%(chr_col,pos)
				sql = sql.replace('\n','')
			
				## prepare statement and query
				cur.execute(sql);
				#process query result
				res = cur.fetchone()
				if res == None:
					res = list()
				if (len(res) > 1):
				# load ediva hash from database
					ediva =  ediva_reselement(ediva,k,res,sep)
				else:
				## handle missing database annotation entry
					ediva[k] =missanndb
			cur.close()					
        ## end of while on not_biallelic_variants                                 
    ## close DB connection
    db.close()

    return ediva
#@@ Check all input and output variables: evaluate if it's better to pass dictionaries
#@@ rather than variables.
#@@ There is much room for improvement in the factorization of code 
## subroutine for Annovar annotation

def fill_annovar(FILE,Annovar, geneDef, comparison):
    with open(FILE) as FILE_pointer:
        for line in FILE_pointer:
            if not(line.startswith("Func")):
                newAnnovarLine = replaceCommainQoute(line)
                dt = newAnnovarLine.rstrip('\n').split(',')
                dt[30].replace('"','')
                if (dt[30].find(';')>=0):
                        annalts = dt[3].split(';')
                        dt[30] = annalts[0]
                valueTOmatch = dt[26]+";"+dt[27]+";"+dt[29]+";"+dt[30]
                valueTOmatch=valueTOmatch.replace('"','')
                ## fix missing values
                if len(dt[0])==0:
                    dt[0] = 'NA'
                if len(dt[1])==0:
                    dt[1] = 'NA'
                if len(dt[2])==0:
                    dt[2] = 'NA'
                if len(dt[3])==0:
                    dt[3] = 'NA'
                annToPass = dt[0]+","+dt[1]+","+dt[2]+","+dt[3]
                annToPass  =annToPass.replace('"','')
                if geneDef == comparison:
                    Annovar[valueTOmatch] = annToPass   # Fill the Annovar Dictionary with values and keys
                else:
                    Annovar[valueTOmatch] = Annovar[valueTOmatch]+ "," + annToPass
    return Annovar


def AnnovarAnnotation(infile,templocation,fileSuffix,geneDef,ANNOVAR,Annovar):
    perl ="/usr/bin/perl "+" "
    ## prepare Annovar input
    annInCmm = perl + ANNOVAR+"/convert2annovar.pl --includeinfo -format vcf4 "+infile+" > "+templocation+"/annInfile"+fileSuffix+"   2> "+infile+".annovar.log"
    print "\t 627 MESSAGE :: Running Annovar command \n > %s" %annInCmm
    subprocess.call(annInCmm,shell=True)
    annFile = templocation+"/annInfile"+fileSuffix+""


    ## run Annovar annotation
    if geneDef == "ensGene":
        command = (ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile "+templocation+"/Ensembl"+
                        fileSuffix+ " "+ templocation+"/annInfile"+fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)
        print "\t 636 MESSAGE :: Running Annovar command \> %s"%command
    elif geneDef == "refGene":
        command =(ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19  --step 1 --outfile "+templocation+"/Refseq"+fileSuffix + " "+
                        templocation+"annInfile"+fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)        
        print "\t 641 MESSAGE :: Running Annovar command \> %s"%command
    elif geneDef == "knownGene":
        command =(ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile "+templocation+"/Known"+fileSuffix + " "+
                        templocation+"/annInfile"+fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)        
        print "\t 646 MESSAGE :: Running Annovar command \> %s"%command
    elif geneDef == "all":
        print "\t 648 MESSAGE :: No sepicific gene definition selected, hence Annovar is going to run on all definitions !\n";
        ## refgene
        command = (ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19   --step 1 --outfile "+templocation+"/Refseq"+fileSuffix + " "+
                        templocation+"/annInfile"+fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)        
        print "\t 653 MESSAGE :: Running Annovar command \> %s"%command
        ## ensgene
        command = (ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19   --genetype ensgene --step 1 --outfile "+templocation+"/Ensembl"+fileSuffix + " "+
                        templocation +"/annInfile"+ fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)        
        print "\t 658 MESSAGE :: Running Annovar command \> %s"%command
        ## knowngene
        command =(ANNOVAR+"/ediva_summarize_annovar.pl --buildver hg19   --genetype knowngene --step 1 --outfile "+templocation+
                        "/Known"+ fileSuffix+" "+templocation+"/annInfile"+ fileSuffix+" "+ANNOVAR+"/hg19/ 2> "+infile+".annovar.log")
        subprocess.call(command,shell=True)        
        print "\t 663 MESSAGE :: Running Annovar command \> %s"%command

    else:
        print "\t MESSAGE :: Not a valid gene definition ! Exiting .."
        finalize(templocation,fileSuffix)

    ## read annotation from annovar
    annFianlAnnE = templocation+"/Ensembl"+fileSuffix+".genome_summary.csv"
    annFianlAnnR = templocation+"/Refseq"+fileSuffix+".genome_summary.csv"
    annFianlAnnK = templocation+"/Known"+fileSuffix+".genome_summary.csv"
    if os.path.isfile(annFianlAnnR):
        Annovar = fill_annovar(annFianlAnnR,Annovar, geneDef, geneDef )
    if os.path.isfile(annFianlAnnE):
        Annovar = fill_annovar(annFianlAnnE,Annovar, geneDef, 'ensGene')
    if os.path.isfile(annFianlAnnK):
        Annovar = fill_annovar(annFianlAnnK,Annovar, geneDef, "knownGene")
    return Annovar

#@@?Almost done, there is to check the substitution value at the end of the routine
## subroutnine por providing header to the main annotation output file
def getHeader(onlygenic,geneDef):
    stringTOreturn=""
    if (onlygenic):
        ## only genic annotation header
        ## check for gene definiton and construct header according to that
        if geneDef == 'ensGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl)" +
                              ",AminoAcidChange(Ensembl),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        elif geneDef == 'refGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),"+
            "ExonicFunction(Refseq),AminoAcidChange(Refseq),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        elif geneDef == 'knownGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Known),Gene(Known),ExonicFunction(Known),"+
            "AminoAcidChange(Known),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        else:
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),"+
            "AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),AminoAcidChange(Ensembl),Function(Known),"+
            "Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
            
    else:
        ## normal annotation header
        ## check for gene definiton and construct header according to that
        if geneDef == 'ensGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),"+
                              "AminoAcidChange(Ensembl),dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,"+
                              "Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,"+
                              "PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,"+
                              "PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,"+
                              "SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        elif geneDef == 'refGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),"+
                              "dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency," +
                              "Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,"+
                              "PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,"+
                              "PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,"+
                              "SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        elif geneDef == 'knownGene':
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),"+
                              "dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,Afr1000GenomesFrequency,"+
                              "Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,"+
                              "PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,"+
                              "Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,"+
                              "samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
        else:
            stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),"+
                              "Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),"+
                              "AminoAcidChange(Known),dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,Afr1000GenomesFrequency,"+
                              "Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,"+
                              "VertebratesPhyloP,PlacentalMammalPhastCons,PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,"+
                              "MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")
    ## replace newlines with nothing at header line
    stringTOreturn.replace('\n','') 
    stringTOreturn.replace(' ','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    return stringTOreturn


#@@ done, same as before: match the /usr/bin/perl expression at the end
## subroutnine por providing header to the inconsistent annotation output file
def getHeaderIns():
    stringTOreturn = ("Chr,Position,Reference,Alteration,AlleleFrequency,GenicAnnotation,dbsnpIdentifier,EurEVSFrequecy,AfrEVSFrequecy,"+
                      "TotalEVSFrequecy,Eur1000GenomesFrequency,Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,"+
                      "Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,"+
                      "PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,"+
                      "SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)")

    stringTOreturn.replace('\n','') 
    stringTOreturn.replace(' ','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    return stringTOreturn


#@@ As before : match /usr/bin/perl output
## subroutnine por providing header to the quick look up mode annotation output file
def getHeaderQlookup():
    stringTOreturn = ("Chr,Position,Reference,Alteration,dbsnpIdentifier,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,"+
                      "Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,"+
                      "PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,"+
                      "SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength")

    ## replace newlines with nothing at header line
    stringTOreturn.replace('\n','') 
    #stringTOreturn.replace('\s+','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    return stringTOreturn

##############################################################################################
## VCF PROCESSING
##############################################################################################
def process_line(myline,i,adindex,gtindex):
    gts = list()
    if ':' in myline[i]:
        gts = myline[i].split(':')
        if adindex != "NF" and ',' in gts[adindex]:
            dps = gts[adindex].split(',')
            dpref = float(dps[0])
            dpalt = float(dps[1])
            if (dpref+dpalt > 0):
                samAf = dpalt/(dpref+dpalt)
                samAf = "{0:.4f}".format(samAf)
                samAf = samAf[:-1]
                samAf =samAf.rstrip('0').rstrip('.')
                #samAf =samAf[:-2].rstrip('0').rstrip('.')
                #samAf = samAf[:5]#(samAf,0,5)
            else:
                samAf = "0.000"
            dpref = "{:.0f}".format(dpref)
            dpalt = "{:.0f}".format(dpalt)
        else:
            dpref = "." #$gts[1];
            dpalt = "."
            samAf = "."
        ## for missing genotype or homozygous reference genotype set the AF to 0
        if gtindex != "NF":
            genotype = gts[gtindex]
        else:
            genotype = "."
    else:
        genotype = myline[i]
        dpref = "."
        dpalt = "."
        samAf = "."
    return genotype,dpref,dpalt,samAf

def samples_fill(samples,gtMode,s_key,s_val,genotype):
    if gtMode == "complete":
        if (samples.get(s_key,False)):
            samples[s_key] += ";"+s_val
        else:
            samples[s_key] = s_val
    else: ## compact
        ## kick out genotypes of '0/0','./.','0|0' and '.|.'
        if genotype != '0/0' and genotype != './.' and genotype != '.|.' and genotype != '0|0':
            if (samples.get(s_key,False)):
                samples[s_key] += ";"+s_val
            else:
                samples[s_key] = s_val
    return samples

def vcf_processing(infile,qlookup,gtMode,type_in):
    allowed_chr = list()
    for i in range(23):
	allowed_chr.append(str(i+1))
    allowed_chr+=(['X','Y','x','y'])
    skipped_chr = list()
    
    variants = dict()
    not_biallelic_variants = dict()
    samples = dict()
    var_counter = 0
    head_count = 0
    comment    =0
    headers = []
    if qlookup == 'NA':#  [later indent the  whole block if needed [depends if the lookup mode is
        # used or not
	tmp = open(infile,'r')
	magic_number = tmp.read(2)
	print 'MESSAGE:: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
	tmp.close()
        with open(infile) if magic_number!='\x1f\x8b'else bgzf.open(infile) as INFILE:
            for line in INFILE:
                if line.startswith('##'):
                    comment+=1
                    pass
                elif line.startswith('#CHROM'): ## grab sample names with other columns from the header line
                    headers = (line.replace('\n','')).split('\t')
                    ## check for malformed VCF file and take action
                    if len(headers) < 8:
                            print "ERROR:: Not a valid VCF format \n 856"
                            raise IOError
                    ## check for genotype format and sample information in the VCF and genotype mode parameter value and take action
                    if (len(headers) < 9 and not(gtMode ==  "none")):
                            print "ERROR:: No Genptype format column and no sample genotype information column present in the VCF. Please run the tool with the --sampleGenotypeMode parameter set to \"none\" \n";
                            raise IOError
                elif 'END=' in line:
		    pass #this are gvcf lines including not relevant annotations
		else: ## data lines
                        var_counter +=1
                        myline = line.rstrip('\n').split('\t')
                        ## check for malformed VCF file and take action
                        if (len(myline) < 8):
                                print "ERROR:: Not a valid VCF format \n 866"
                                raise IOError
                        chr_col     = myline[0]
                        position    = myline[1]
                        ref         = myline[3]
                        alt         = myline[4]
			alt_tmp = alt.split(',')
			while '<NON_REF>' in alt_tmp: alt_tmp.remove('<NON_REF>')
			alt = ','.join(alt_tmp)
                        infos       = myline[7].split(';')
                        AF          = "."
                        gtindex     = "NF"
                        adindex     = "NF"
                        ## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
                        if chr_col.startswith('chr') or chr_col.startswith('Chr'):
                                chr_col = chr_col[3:]
                        ## take care of chr 23 or 24 and convert to X or Y
                        if chr_col == "23":
                            chr_col='X'
                        if chr_col == "24":
                            chr_col='Y'
                        if chr_col == "25":
                            chr_col='MT'
                        keywords = ['M','T','m','t']
                        keywords_alt = ['A','T','G','C','-',',']
                        keychars = ['n','N']
                        if not chr_col in allowed_chr:
                            skipped_chr.append(chr_col)
                        #if any(k in chr_col for k in keywords):
                        #        print( "WARNING:: !! does not support chromosome %s currently. This variant line will be skipped in the final annotation output " % (chr_col))
                        #elif not(any( k in alt for k in keywords_alt)):
                        elif not(all(k in keywords_alt for k in alt)):
                                print( "WARNING:: Unknown alternate allele detected at %s and %s. This variant line will be skipped in the final annotation output " %( chr_col,position))
                        else:
                        ## grab the AF from the INFO field in the VCF
                            for info in infos:
                                if info.startswith('AF='):
                                    AF = info[3:]
                                    break
                            ## confirm AF extraction from the info tags
                            if AF == ".":
                                print "WARNING:: AF tag not found in the INFO column at chromosome %s and position %s. AF will be set to \".\" for this variant"%(chr_col,position)
				print "WARNING:: Info field: %s\n"%infos
                            ## always test for complete genotype format field consistency in the VCF; if abnormal report for that variant
                            if gtMode != "none":
                                if len(myline) > 8 and ':' in myline[8] :
                                    gtcheck = myline[8].split(':')
                                    for gti in range(0,len(gtcheck)):
                                            if gtcheck[gti] == "GT":
                                                    gtindex = gti
                                            if gtcheck[gti] == "AD":
                                                    adindex = gti
                                    ## check for the GT and AD fields
                                    if gtindex == "NF":
                                        print ("WARNING:: Weird genotype format %s found in $input at chromosome %s and position %s. No GT field present in %s \n"
                                               %(myline[8],chr_col,position,myline[8] ))
                                    if adindex == "NF":
                                        print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. No AD field present in %s \n"
                                        %(myline[8], infile,chr_col,position,myline[8]))
                                else:
                                    print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. Expected genotype format field separator \":\" \n"
                                           %(myline[8],infile,chr_col,position))
                        ## process based on alteration
			    if ',' in alt : ## section for all sites other than bi-allelic
				alts = alt.split(',')
				afs = AF.split(',')
				if len(afs)<len(alts):
				    afs=['.']*len(alts)
				    
				## process each alteration allele
				for j in range(0,len(alts)):
					al = alts[j]
					alfr = afs[j]
					 
					## decide for variant type
					lenref = len(ref)
					lenalt = len(al)
					if (lenref + lenalt) > 2:## INDEL
					    token_ref = "NA"
					    token_obs = "NA"
					    ## make indelID
					    ## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
					    if not(any(k in ref for k in keychars)) and not(any(k in al for k in keychars)): #($ref !~ m/[Nn]/ and $al !~ m/[Nn]/)
					    #try:
						hash_ref = hashlib.md5(str(ref).encode())
						hash_al = hashlib.md5(str(al).encode())
						token_ref = str(struct.unpack('<L', hash_ref.digest()[:4])[0])
						token_obs = str(struct.unpack('<L', hash_al.digest()[:4])[0]) ## hope it's what's supposed to do unpack('L', md5($ref));
						#except:
						#    print myline
						#    #print hash_ref
						#    print ref
						#    print alt
						#    print hash_ref.digest()
						#    hash_ref = hashlib.md5(str(ref).encode())
						#    raise
					    if type_in == 'INDEL' or type_in == 'all':
						   ## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
						   # e.g A,C in the alternate column in VCF will report only A in the main annotation file
						   ## we are doing this because we want to keep the annotation main file consistent
						    if j == 0 and token_ref != "NA" and token_obs != "NA":   #HERE WE FILL THE VARIANTS DICTIONARY
							variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+al+';'+alfr
						    else:
							not_biallelic_variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+al+';'+alfr
						     ## if sample wise information is present in the VCF then process ; otherwise skip
						     ## also check for none value in genotype mode parameter
						    if (len(myline) > 8) and (gtMode != "none"):
							for i in range(9,len(myline)):
							    (genotype,dpref,dpalt,samAf) = process_line(myline,i,adindex,gtindex)
							    ## check for sample genotype mode
							    s_key =chr_col+';'+position+';'+token_ref+';'+token_obs
							    s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
							    samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
					else: ## SNP
					    ## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
					    ## e.g A,C in the alternate column in VCF will report only A in the main annotation file
					    ## we are doing this because we want to keep the annotation main file consistent
					    print 'snp : %s ' % type_in
					    if type_in == 'SNP' or type_in == 'all':
						## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
						if j == 0 and  not(any(k in ref for k in keychars)) and not(any(k in al for k in keychars) ):
						    variants[ chr_col+';'+position+';'+ref+';'+al] = chr_col+';'+position+';'+ref+';'+al+';'+alfr
						else:
						    not_biallelic_variants[ chr_col+';'+position+';'+ref+';'+al] = chr_col+';'+position+';'+ref+';'+al+';'+alfr
					    ## if sample wise information is present in the VCF then process ; otherwise skip
					    ## also check for none value in genotype mode parameter
						if len(myline) > 8 and gtMode != "none":
						    for i in range(9,len(myline)):
							(genotype,dpref,dpalt,samAf) = process_line(myline,i,adindex,gtindex)
							## check for sample genotype mode
							s_key =chr_col+';'+position+';'+ref+';'+al
							s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
							samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
			    else: ## section for bi-allelic sites
				## decide for variant type
				lenref = len(ref)
				lenalt = len(alt)
				if (lenref + lenalt) > 2: ## INDEL
				    token_ref = "NA"
				    token_obs = "NA"
				    ## make indelID
				    ## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
				    if not(any(k in ref for k in keychars)) and not(any(k in alt for k in keychars)):
					# Assumes the default UTF-8
					hash_ref = hashlib.md5(str(ref).encode())
					hash_alt = hashlib.md5(str(alt).encode())
					token_ref = str(struct.unpack('<L', hash_ref.digest()[:4])[0]) ## hope it's what's supposed to do unpack('L', md5($ref));
					token_obs = str(struct.unpack('<L', hash_alt.digest()[:4])[0]) ## hope it's what's supposed to do unpack('L', md5($ref));
					#token_ref = ref
					#token_obs = alt
				    if type_in == 'INDEL' or type_in == 'all':
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
					if token_ref != "NA" and token_obs != "NA":
					    variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+alt+';'+AF
					else:
					    not_biallelic_variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+alt+';'+AF
					## if sample wise information is present in the VCF then process ; otherwise skip
					## also check for none value in genotype mode parameter
					if len(myline) > 8 and gtMode != "none":
					    for i in range(9,len(myline)):
						(genotype,dpref,dpalt,samAf) = process_line(myline,i,adindex,gtindex)
						## check for sample genotype mode
						s_key =chr_col+';'+position+';'+token_ref+';'+token_obs
						s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
						samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
				else:
				    ## SNP
				    if type_in == 'SNP' or type_in == 'all':
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
					if not(any(k in ref for k in keychars)) and not(any(k in alt for k in keychars)):
					    variants[ chr_col+';'+position+';'+ref+';'+alt] = chr_col+';'+position+';'+ref+';'+alt+';'+AF
					else:
					    not_biallelic_variants[ chr_col+';'+position+';'+ref+';'+alt] = chr_col+';'+position+';'+ref+';'+alt+';'+AF
					## if sample wise information is present in the VCF then process ; otherwise skip
					## also check for none value in genotype mode parameter
					if len(myline) > 8 and gtMode != "none":
					    for i in range(9,len(myline)):
						(genotype,dpref,dpalt,samAf) = process_line(myline,i,adindex,gtindex)
						## check for sample genotype mode
						s_key =chr_col+';'+position+';'+ref+';'+alt
						s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
						samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
	    #Reporting skipped lines
            if len(skipped_chr)>0:
                tmp =  dict((i,skipped_chr.count(i)) for i in set(skipped_chr))
                print "WARNING: there have been skipped lines due to unknown chr:"
                for i in list(set(skipped_chr)):
                   print "\tChr: %s - %d lines"%(i,tmp.get(i))
    else :
        ## Quick lookup mode
        ## check for input type in -q parameter
        ## decide for file or query string
        if os.path.isfile(qlookup):
            ## its a file
            # print "MESSAGE :: Processing input file - %s " %qlookup
            with open(qlookup) as FL:
                for line in FL:
                    var = line.rstrip('\n').split(':')
                    ## check for quick lookup data field consistency
                    if len(var) != 4:
                        print "ERROR :: Not a valid format. Correct format is chr:position:reference:alternate "
                        raise IOError
                    else:
                        ## check for simple checking of the quick loookup data fields
                        if re.search('[a-zA-Z]+',var[1]):
                            print "\ERROR :: Not a valid position value "
                            raise IOError
                        elif re.search('[0-9Nn]+',var[2]):
                            print "ERROR :: Not a valid reference allele value "
                            raise IOError
                        elif re.search('[0-9Nn]+',var[3]):
                            print "ERROR :: Not a valid alternate allele value "
                            raise IOError
                        else:
                            ## assign variation from quick look up format
                            chr_col = var[0]
                            pos     = var[1]
                            ref     = var[2]
                            alt     = var[3]
                            ## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
                            if chr_col.startswith('chr') or chr_col.startswith('Chr'):
                                chr_col = chr_col[3:]
                            ## take care of chr 23 or 24 and convert to X or Y
                            if chr_col == "23":
                                chr_col='X'
                            if chr_col == "24":
                                chr_col='Y'
                            if chr_col == "25":
                                chr_col='MT'
                            keychars = ['n','N']
                            keywords = ['M','T','m','t']
                            keywords_alt = ['A','T','G','C','-']
#                            if any(k in chr_col for k in keywords):
#				print( "WARNING:: !! does not support chromosome %s currently. This variant line will be skipped in the final annotation output " % (chr_col))
			    if not chr_col in allowed_chr:
				skipped_chr.append(chr_col)
                                #elif not(any( k in alt for k in keywords_alt)):
                            elif not(all(k in keywords_alt for k in alt)):
				print( "WARNING:: Unknown alternate allele detected at %s and %s. This variant line will be skipped in the final annotation output " %( chr_col,pos))
				print 'line : 1151'
				print alt
				raise
                            else:
                            ## decide for variant type
				lenref = len(ref)
				lenalt = len(alt)
				if (lenref + lenalt) > 2 :## INDEL
				    token_ref = "NA"
				    token_obs = "NA"
				    ## make indelID
				    ## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
				    if not(any(k in ref for k in keychars)) and not(any(k in alt for k in keychars)): #($ref !~ m/[Nn]/ and $al !~ m/[Nn]/)
					hash_ref = hashlib.md5(str(ref).encode())
					hash_alt = hashlib.md5(str(alt).encode())
					token_ref = str(struct.unpack('<L', hash_ref.digest()[:4])[0])
					token_obs = str(struct.unpack('<L', hash_alt.digest()[:4])[0]) ## hope it's what's supposed to do unpack('L', md5($ref));
				    variants[chr_col+';'+pos+';'+token_ref+';'+token_obs] = chr_col+';'+pos+';'+ref+';'+alt+';.'
				else: ## SNP
				    variants[chr_col+';'+pos+';'+ref+';'+alt]= chr_col+';'+pos+';'+ref+';'+alt+';.'
		#Reporting skipped lines
		if len(skipped_chr)>0:
		    tmp =  dict((i,skipped_chr.count(i)) for i in set(skipped_chr))
		    print "WARNING: there have been skipped lines due to unknown chr:"
		    for i in list(set(skipped_chr)):
			print "\tChr: %s - %d lines"%(i,tmp.get(i))
	else:
	    ## query string
	    var = qlookup.rstrip('\n').split(':')
	    if len(var) != 4:
		print "ERROR :: Not a valid format. Correct format is chr:position:reference:alteration "
		raise IOError
	    else:
		## check for simple checking of the quick loookup data fields
		if  re.search('[a-zA-Z]+',var[1]):
		    print "\ERROR :: Not a valid position value "
		    raise IOError
		elif re.search('[0-9Nn]+',var[2]):
		    print "ERROR :: Not a valid reference allele value "
		    raise IOError
		elif re.search('[0-9Nn]+',var[3]):
		    print "ERROR :: Not a valid alternate allele value "
		    raise IOError
		else:
		    ## assign variation from quick look up format
		    chr_col = var[0]
		    pos     = var[1]
		    ref     = var[2]
		    alt     = var[3]
		    ## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
		    if chr_col.startswith('chr') or chr_col.startswith('Chr'):
		        chr_col = chr_col[3:]
			## take care of chr 23 or 24 and convert to X or Y
		    if chr_col == "23":
			chr_col='X'
		    if chr_col == "24":
			chr_col='Y'
		    if chr_col == "25":
			chr_col='MT'
		    keychars = ['n','N']
		    keywords = ['M','T','m','t']
		    keywords_alt = ['A','T','G','C','-']
		    if any(k in chr_col for k in keywords):
			print( "WARNING:: !! does not support chromosome %s currently. This variant line will be skipped in the final annotation output \n" % (chr_col))
		    #elif not(any( k in alt for k in keywords_alt)):
                    elif not(all(k in keywords_alt for k in alt)):
			print( "WARNING:: Unknown alternate allele detected at %s and $position. This variant line will be skipped in the final annotation output \n" %( chr_col))
			print 'line 1218'
			print alt
			raise
		    else:
			lenref = len(ref)
			lenalt = len(alt)
	
			if (lenref + lenalt) > 2 :## INDEL
			    token_ref = "NA"
			    token_obs = "NA"
			    ## make indelID
			    ## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
			    if not(any(k in ref for k in keychars)) and not(any(k in alt for k in keychars)): #($ref !~ m/[Nn]/ and $al !~ m/[Nn]/)
				hash_ref = hashlib.md5(str(ref).encode())
				hash_alt = hashlib.md5(str(alt).encode())
				token_ref = str(struct.unpack('<L', hash_ref.digest()[:4])[0])
				token_obs = str(struct.unpack('<L', hash_alt.digest()[:4])[0]) ## hope it's what's supposed to do unpack('L', md5($ref));
			    variants[chr_col+';'+pos+';'+token_ref+';'+token_obs] = chr_col+';'+pos+';'+ref+';'+alt+';.'
			else: ## SNP
			    variants[chr_col+';'+pos+';'+ref+';'+alt]= chr_col+';'+pos+';'+ref+';'+alt+';.'

    return samples,variants,not_biallelic_variants