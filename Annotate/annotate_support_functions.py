
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
#import y_serial_v060 as y_serial


######################################
#
#       Task => Annotate variants
#       infile => VCF file with complete sample wise genotype information
#       outfile => text file (csv) with complete annotation with sample wise genotype information
#       Extra outfile => text file (csv) without genic annotation with sample wise information for variants that are not bi-allelic (e.g tri-allelic)
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
                 'vcflike: sample information is stored all in one column divided by ":" \n\t\t\t '+
                 'none: exclude sample wise genotype information in annotation \n\t\t\t default: complete')
o_str = 'If set, then only genic annotation will be performed'
f_str = 'If set, then it will over-write existing output annotation file with the same name '
## subroutine for usage of the tool @@
def usage():
    print ("\n usage:\n" +
    "--input,-i \t\t"                   + i_str + " \n\n" +
    "--quicklookup,-q \t "              + q_str +  "\n\n" +
    "--tempDir,-t \t\t "                + t_str + " \n\n" +
    "--geneDef,-g \t\t "                + g_str + " \n\n" +
    "--variantType,-v \t "              +v_str  + " \n\n" +
    "--sampleGenotypeMode,-s  " + s_str + " \n\n" +
    "--onlyGenicAnnotation,-o " +o_str  + " \n\n" +
    "--forceNewFileCreate,-f  " +f_str  + " \n\n" +
    "--help,\t\t show help \n\n" )
    raise IOError
##############################################################################################
## COMMAND LINE OPTIONS done- to be tested
##############################################################################################
def input_parse(defaults):
    ''' Parses the input arguments for further processing'''
    parser_ = defaults
        ## grab command line options
    parser = argparse.ArgumentParser(description = 'Setup ')
    parser.add_argument('-i','--input'  ,            type=str, dest ="infile",       required=False,                                     help= i_str)
    parser.add_argument('-t','--tempDir',            type=str, dest ="templocation", required=False,                                     help= t_str)
    parser.add_argument('-v','--variantType'  ,      type=str, dest ="type",         required=False,                                     help= v_str)
    parser.add_argument('-o','--onlyGenicAnnotation',              dest ="onlygenic",    required=False, action='store_true',help= o_str)
    parser.add_argument('-f','--forceNewFileCreate',               dest ="forcedel",     required=False, action='store_true',help= f_str)
    parser.add_argument('-q','--quicklookup'  ,      type=str, dest ="qlookup",      required=False,                                     help= q_str)
    parser.add_argument('-s','--sampleGenotypeMode', type=str, dest ="gtmode",       required=False,                                     help= s_str)
    parser.add_argument('-g','--geneDef',            type=str, dest ="geneDef",      required=False,                                     help= g_str)
    parser.add_argument('--csvfile', dest='csvfile', required=False, help='csv file with username and user email address. [optional]')
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
    if args.csvfile:
        parser_["csvfile"] = args.csvfile

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
    if parser_["gtmode"]!= "complete" and parser_["gtmode"] != "compact" and parser_["gtmode"] != "vcflike" and parser_["gtmode"] != "none":
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
    ''' Checks for Annovar settings'''
    ### check of annovar settings
    vcftoAnn = ANNOVAR+"/convert2annovar.pl"
    genicAnn = ANNOVAR+"/ediva_summarize_annovar.pl"
    maftoAnn = ANNOVAR+"/maf2annovar.pl"
    if not(os.path.isdir(ANNOVAR)):    #{
        print "\nERROR :: Annovar library location is not found. Please set the Annovar library path correctly inside the program \n"
        raise IOError
    if not(os.path.exists(maftoAnn)):
        print "\nERROR :: Program for converting MAF to Annovar format is not found in the Annovar library location \n"
        raise IOError

    if not(os.path.exists(genicAnn)):
        print "\nERROR :: Program for converting VCF to Annovar format is not found in the Annovar library location \n";
        raise IOError
    if not(os.path.exists(genicAnn)):
        print "\nERROR :: Program for converting VCF to Annovar format is not found in the Annovar library location \n";
        raise IOError

    return (vcftoAnn,genicAnn)

##############################################################################################
## OUTPUT FILE(s) "done"
##############################################################################################
def out_file_generate(infile,qlookup,templocation,forceDel,tempfile,MAF):
    ''' Generates the output files for annotation'''
    outFile         = ''
    sortedOutFile   = ''
    outFileIns      = ''
    filename = infile
    if qlookup == "NA":
        filename = os.path.basename(infile)[:-4]
        pathname = os.path.dirname(infile)
        if len(pathname)>1:
            pathname = pathname +'/'
            if templocation == "INPATH":
                templocation = pathname
        elif templocation == "INPATH":
            templocation = "."

        if MAF==0:

            outFile  = pathname+ filename+".annotated.csv"
            sortedOutFile = pathname+ filename +".sorted.annotated.csv"
            outFileIns = pathname+filename + ".inconsistent.annotated.csv"
        else:
            outFile  = pathname+ filename+".annotated.maf"
            sortedOutFile = pathname+ filename +".sorted.annotated.maf"
            outFileIns = pathname+filename + ".inconsistent.annotated.maf"
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
            #filename = '.'.join(filename.split('.')[:-1])
            pathname = ntpath.dirname(qlookup)
            #filename = ntpath.basename(qlookup)[:-4]
            #pathname = ntpath.dirname(qlookup)
            if len(pathname)>1:
                pathname = pathname +'/'
            if MAF==0:
                outFile  = pathname+ filename+".annotated.csv"
                sortedOutFile = pathname+ filename +".sorted.annotated.csv"
                outFileIns = pathname+filename + ".inconsistent.annotated.csv"
            else:
                outFile  = pathname+ filename+".annotated.maf"
                sortedOutFile = pathname+ filename +".sorted.annotated.maf"
                outFileIns = pathname+filename + ".inconsistent.annotated.maf"
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
    if MAF ==0:
        print("MESSAGE :: Your sorted annotated file is: %s" %(sortedOutFile))
        #print("MESSAGE :: Reported non bi-allelic sites are in: %s" %(outFileIns))
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
    ''' clear the tmp directory for this session'''
    clearCmm = "rm -r " + templocation + "/*"+fileSuffix +"*"
    #print(clearCmm)
    subprocess.call(clearCmm,shell = True)
    return None

## sub for preparing missing db annotation @@
def preparemissdb(sep):
    '''sub for preparing missing db annotation  '''
    missandb            = 'NA'
    missandb_coordinate = '0'
    missanndbindel       = 'NA'

    missandb = missandb + (sep+"0")*9
    missandb = missandb + (sep+"NA")*16#It was 14!! remember this because we changed it to include segrepeat and EIGEN
    missandb_coordinate = missandb_coordinate + (sep+"NA")*10
    missanndbindel = missanndbindel + (sep+"0")*8
    #print missandb


    #raise
    return(missandb,missandb_coordinate,missanndbindel)
#@@ should be done: check it a couple of times with examples
## sub for replacing commas inside double qoutes for annovar genic annotation lines
def replaceCommainQoute(in_str):
    '''sub for replacing commas inside double qoutes for annovar genic annotation lines'''
    split_str = in_str.split('"');
    tmplist =[]
    token = 1

    for s in split_str[token:-1:2]:
        tmplist.append(s.replace(',',';'))
    split_str[token:-1:2] = tmplist
    return('"'.join(split_str))



#@@ dumps the current dict content to the databse file specified in db
#@@ in the table specified via table
#@@ it uses y_serial to do it in a sqlite file
def dump_to_db(dict_to_dump,db_str,tablename):
    db_str= os.path.abspath(db_str)
    dir_str= os.path.dirname(db_str)
    ## problematic line here
    subprocess.call(['touch',db_str])
    subprocess.call(['chmod',  '777',db_str ])

    db = y_serial.Main(db_str)
    todump = list()
    for k,v in dict_to_dump.items():
        todump.append([v,k.replace(';','_')])
        #todump.append(['a','b'])
    print tablename
    db.ingenerator(todump,tablename)
    del dict_to_dump
    dict_to_dump=dict()
    #subprocess.call(['chmod',  '774',dir_str ])
    return dict_to_dump


#@@Done
## subroutine for ediva public omics data fetch
def edivaPublicOmics():
    ''' subroutine for ediva public omics data fetch'''
    edivaStr = dict()
    ## DB parameters
    username    = "edivacrg"
    database    = "eDiVa_public_omics"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "FD5KrT3q"

    ## DB parameters
    username = 'edivapublic';
    database = 'eDiVa_annotation';
    dbhost = 'mysqlsrv-ediva.linux.crg.es';
    passw = 'x86d2k1B';

    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base

    cur = db.cursor()
    sql = "select chr,pos,lengthofrepeat,region from eDiVa_public_omics.Table_simpleRepeat;"


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


def edivaPublicOmics_search(chr_,pos):
    out = 'NA,NA'
    username    = "edivacrg"
    database    = "eDiVa_public_omics"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "FD5KrT3q"

    ## DB parameters
    username = 'edivapublic';
    database = 'eDiVa_annotation';
    dbhost = 'mysqlsrv-ediva.linux.crg.es';
    passw = 'x86d2k1B';

    db = MySQLdb.connect(host=dbhost, # your host
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base

    cur = db.cursor()
    sql = "select lengthofrepeat,region from eDiVa_public_omics.Table_simpleRepeat where chr =%s and pos = %s limit 1"%(chr_,pos)
    cur.execute(sql)

    for row in cur:
        out = str(row[1])+','+str(row[0])
    cur.close()
    db.close()
    #print out

    return out
#@@ Done
def process_db_entry(res,res2,ediva,k,sep,missanndb_coordinate,missanndbindel,exac_value= None,clinvar_value=None):
    ''' Processes the database entry, fetching values and output the annotated line'''
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
        #ediva[k] +=  ((sep+ "NA")*6)
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
        #ediva[k] += ((sep+ "NA")*6)
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
        #ediva[k] += ((sep+ "NA")*6)
    else:
        ## no entry in the database for both the queries
        if (ediva.get(k,False)):
            ## this should never happen !!
            ediva[k] += (sep+missanndbindel)
        else:
            ediva[k] =   (missanndbindel)
        # take care of missing positional values
        ediva[k] +=(sep+missanndb_coordinate)
    ##extract the repeat region value and id because they go after exac
    temp=ediva[k].split(sep)
    repeat_region_info = sep.join(temp[-2:])
    ediva[k]=sep.join(temp[:-2])

    ## add NAs for damage potential scores for indels and ABB score Now just for SNPs
    ediva[k] +=((sep+ "NA")*8+sep+'0')
    if exac_value != None:
        #there is a result:
        for ex_element in exac_value:

            ediva[k] += (sep+str(ex_element))

    else:
        ediva[k] += ((sep+ "0")*9)

    ## Clinvar part
    if clinvar_value != None:
        #there is a result:
        for ex_element in clinvar_value:

            ediva[k] += (sep+str(ex_element))

    else:
        ediva[k] += ((sep+ "NA")*2)

    ediva[k]+=sep+repeat_region_info

    return ediva

def ediva_reselement(ediva,k,res,sep,exac_value=None,clinvar_value=None):
    '''Populate the ediva annotation dictionary '''
    if ediva.get(k,False):
        pass
    else:
        for reselement in res:
            if ediva.get(k,False):
                ediva[k] += sep+str(reselement)
            else:
                ediva[k]=str(reselement)
        ##extract the repeat region value and id because they go after exac
        temp=ediva[k].split(sep)
        repeat_region_info = sep.join(temp[-2:])
        ediva[k]=sep.join(temp[:-2])
        if exac_value==None:
            ediva[k] += ((sep+ "0")*9)
        elif exac_value=='indel':
            pass
        else:
            for ex_element in exac_value:
                ediva[k] += (sep+str(ex_element))
        ## Clinvar part
        if clinvar_value != None:
            #there is a result:
            for ex_element in clinvar_value:

                ediva[k] += (sep+str(ex_element))

        else:
            ediva[k] += ((sep+ "NA")*2)

        ##aggregate repeat region info after the exac value and clinvar
        ediva[k]+=sep+repeat_region_info
    return ediva

#@@ Takes care to do the queries for the ediva database for SNP and INDEl
#@@ IT calls reselement and process_db_entry to fix the format where needed
def query_db(k,v,ediva,db,sep,missanndb,missanndb_coordinate,missanndbindel):
    chr_col,pos,ref,alt = k.rstrip('\n').split(';')
    step1_common_items="""
    SELECT ifnull(varinfo.dbsnpid,'NA'),
        ifnull(varinfo.EurEVSFreq,'0'),
        ifnull(varinfo.AfrEVSFreq,'0'),
        ifnull(varinfo.TotalEVSFreq,'0'),
        ifnull(varinfo.EurAFKG,'0'),
        ifnull(varinfo.AfrAFKG,'0'),
        ifnull(varinfo.AsaAFKG,'0'),
        ifnull(varinfo.AmrAFKG,'0'),
        ifnull(varinfo.AFKG,'0')        """

    snp_1_specific=",\n ifnull(varinfo.SDSnp,'0'),"
    indel_1_specific="SELECT ifnull(varinfo.SDIndel,'0'),"

    step2_common_items="""
        ifnull(varinfo.`placentalMammal.phyloP`,'NA'),
        ifnull(varinfo.`primates.phyloP`,'NA'),
        ifnull(varinfo.`vertebrates.phyloP`,'NA'),
        ifnull(varinfo.`placentalMammal.phastCons`,'NA'),
        ifnull(varinfo.`primates.phastCons`,'NA'),
        ifnull(varinfo.`vertebrates.phastCons`,'NA'),
        ifnull(varinfo.gerp1,'NA'),
        ifnull(varinfo.gerp2,'NA'),
    """
    snp_ending ="""
        ifnull(varinfo.sift,'NA'),
        ifnull(varinfo.polyphen2,'NA'),
        ifnull(varinfo.mutationassessor,'NA'),
        ifnull(varinfo.condel,'NA'),
        ifnull(varinfo.cadd1,'NA'),
        ifnull(varinfo.cadd2,'NA'),
        ifnull(eigen.Eigen_raw,'NA'),
        ifnull(eigen.Eigen_phred,'NA'),
        ifnull(abb.ABB_score,'0'),"""

    exac_query="""  SELECT ifnull(exac.AF,'0'),
                    ifnull(exac.AF_adj,'0'),
                    ifnull(exac.AF_AFR,'0'),
                    ifnull(exac.AF_AMR,'0'),
                    ifnull(exac.AF_EAS,'0'),
                    ifnull(exac.AF_FIN,'0'),
                    ifnull(exac.AF_NFE,'0'),
                    ifnull(exac.AF_OTH,'0'),
                    ifnull(exac.AF_SAS,'0')"""
    clivar_query=""" SELECT ifnull(clinvar.clinical_significance,'NA'),
                     ifnull(clinvar.access_number,'NA') """

    repeat_query="""
                    ifnull(simplerep.lengthofrepeat,'NA'),
                    ifnull(simplerep.region,'NA')
                    """

    chr_col,pos,ref,alt = k.rstrip('\n').split(';')
    ## decide for variant type
    lenref = len(ref)
    lenalt = len(alt)
    sql     = ""
    sql2    = ""
    stmt    = ""
    stmt2   = ""
    res     = list()
    res2    = list()
    if ((lenref + lenalt) > 2):## INDEL
        sql=step1_common_items + """
        FROM eDiVa_annotation.Table_Chr%s_indel as varinfo WHERE indelid = '%s' LIMIT 1;"""%(chr_col,k)


        sql2=indel_1_specific + step2_common_items +   repeat_query+ """
        FROM Table_Chr%s as varinfo
        LEFT JOIN eDiVa_public_omics.Table_simpleRepeat AS simplerep ON (
            varinfo.chromosome = simplerep.chr AND
            varinfo.position = simplerep.pos AND
            varinfo.position = '%s'

        )
        WHERE varinfo.position = '%s' LIMIT 1;"""%(chr_col,pos,pos)
        exac_sql=""
        exac_sql=exac_query + """
        FROM `exac_indel` as exac WHERE exac.indelid= '%s' LIMIT 1;"""%(k)

        clinvar_sql=clivar_query +  """
            FROM eDiVa_annotation.Table_clinvar as clinvar WHERE clinvar.chr = '%s'
            AND clinvar.pos = %s AND clinvar.ref = '%s'
            AND clinvar.alt = '%s'
            LIMIT 1;"""%(chr_col,pos,ref,alt)


        ## prepare statement and query#######################
        #print sql2
        #raise
        sql = sql.replace('\n',' ')
        sql2 = sql2.replace('\n',' ')
        sql = sql.replace('\t',' ')
        sql2 = sql2.replace('\t',' ')
        exac_sql = exac_sql.replace('\n',' ')
        exac_sql = exac_sql.replace('\t',' ')
        clinvar_sql = clinvar_sql.replace('\n',' ')
        clinvar_sql = clinvar_sql.replace('\t',' ')
        cur = db.cursor()
        cur2 = db.cursor()
        exac_cur=db.cursor()
        clinvar_cur=db.cursor()
        cur.execute(sql)
        res = cur.fetchone()
        cur.close()
        cur2.execute(sql2)
        #process query result
        res2 = cur2.fetchone()
        cur2.close()


        #exac query execution
        exac_cur.execute(exac_sql)
        exac_result= exac_cur.fetchone()
        exac_cur.close()

        #clinvar query execution
        clinvar_cur.execute(clinvar_sql)
        clinvar_result = clinvar_cur.fetchone()
        clinvar_cur.close()

        ediva  = process_db_entry(res,res2,ediva,k,sep,missanndb_coordinate,missanndbindel,exac_result,clinvar_result)
    else:
        stringent= """ WHERE varinfo.position = %s AND varinfo.Reference = '%s' AND
        varinfo.Alt = '%s'
            LIMIT 1;"""%(pos,ref,alt)
        broad = """WHERE varinfo.position = %s LIMIT 1;"""%(pos)
        sql = step1_common_items + snp_1_specific + step2_common_items + snp_ending+  repeat_query + """
        FROM eDiVa_annotation.Table_Chr%s  as varinfo
        LEFT JOIN eDiVa_annotation.Eigein_coding AS eigen ON (
            varinfo.chromosome = eigen.chromosome AND
            varinfo.position = eigen.position AND
            varinfo.Reference = eigen.Reference AND
            varinfo.Alt = eigen.Alt AND
            varinfo.position = '%s')
        LEFT JOIN eDiVa_public_omics.ABB_score AS abb ON(
            varinfo.chromosome = abb.chr AND
            varinfo.position = abb.pos
        )
        LEFT JOIN eDiVa_public_omics.Table_simpleRepeat AS simplerep ON (
            varinfo.chromosome = simplerep.chr AND
            varinfo.position = simplerep.pos
        )
        """%(chr_col,pos) + stringent
        #print sql

        sql = sql.replace('\n',' ')
        sql = sql.replace('\t',' ')
        cur = db.cursor()
        #raise

        cur.execute(sql)
        res = cur.fetchone()
        if res is None:
            res = list()
        ## prepare statement and query#######################
        cur.close()
        exac_sql = exac_sql=exac_query + """
            FROM eDiVa_annotation.exac_snp as exac WHERE exac.chromosome = '%s'
            AND exac.position = %s AND exac.Reference = '%s'
            AND exac.Alt = '%s'
            LIMIT 1;"""%(chr_col,pos,ref,alt)

        exac_sql = exac_sql.replace('\n',' ')
        exac_sql = exac_sql.replace('\t',' ')
        exac_cur=db.cursor()
        #print exac_sql
        exac_cur.execute(exac_sql)
        exac_result= exac_cur.fetchone()
        exac_cur.close()


        clinvar_sql=clivar_query +  """
            FROM eDiVa_annotation.Table_clinvar as clinvar WHERE clinvar.chr = '%s'
            AND clinvar.pos = %s AND clinvar.ref = '%s'
            AND clinvar.alt = '%s'
            LIMIT 1;"""%(chr_col,pos,ref,alt)
        clinvar_sql = clinvar_sql.replace('\n',' ')
        clinvar_sql = clinvar_sql.replace('\t',' ')
        clinvar_cur=db.cursor()
        #print exac_sql
        clinvar_cur.execute(clinvar_sql)
        clinvar_result= clinvar_cur.fetchone()
        clinvar_cur.close()

        if (len(res) > 1):
            # load ediva hash from database
            ediva =  ediva_reselement(ediva,k,res,sep,exac_result,clinvar_result)
        else:
            #Less stringent search driven only by chr and position
            sql = step1_common_items + snp_1_specific+ step2_common_items +snp_ending +  repeat_query + """
        FROM eDiVa_annotation.Table_Chr%s  as varinfo
        LEFT JOIN eDiVa_annotation.Eigein_coding AS eigen ON (
            varinfo.chromosome = eigen.chromosome AND
            varinfo.position = eigen.position AND
            varinfo.Reference = eigen.Reference AND
            varinfo.Alt = eigen.Alt AND
            varinfo.position = '%s' )
        LEFT JOIN eDiVa_public_omics.ABB_score AS abb ON(
            varinfo.chromosome = abb.chr AND
            varinfo.position = abb.pos
        )
        LEFT JOIN eDiVa_public_omics.Table_simpleRepeat AS simplerep ON (
            varinfo.chromosome = simplerep.chr AND
            varinfo.position = simplerep.pos
        )
        """%(chr_col,pos) + broad
            #print sql
        #if (len(res) > 1):
        #    # load ediva hash from database
        #    ediva =  ediva_reselement(ediva,k,res,sep,exac_result)
        #else:
        #    #Less stringent search driven only by chr and position
        #    sql = sqlbase+broad
            sql = sql.replace('\n',' ')
            sql = sql.replace('\t',' ')
            ## prepare statement and query#######################
            cur = db.cursor()
            cur.execute(sql)

            #process query result
            res = cur.fetchone()
            if res is None:
                res = list()
            ## prepare statement and query#######################
            cur.close()
            if (len(res) > 1):
                # load ediva hash from database
                ediva =  ediva_reselement(ediva,k,res,sep,exac_result,clinvar_result)
            else:
                ## handle missing database annotation entry
                ediva[k]=missanndb
                ediva[k]+= ((sep+ "0")*10)# ABB +EXAC
                ediva[k] += ((sep + "NA")*4) #Clinvar+Tandem repeat

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

    ## DB parameters
    username = 'edivapublic';
    database = 'eDiVa_annotation';
    dbhost = 'mysqlsrv-ediva.linux.crg.es';
    passw = 'x86d2k1B';

    ## open DB connection
    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw, # your password
    db=database) # name of the data base
    # extract db result
    print("\t Running ediva annotaton with %d variants : "%(len(variants.items())))
    for k, v in (variants.items()) :
        ediva=query_db(k,v,ediva,db,sep,missanndb,missanndb_coordinate,missanndbindel)

    ## end of if-else for variant decision
    # extract db result
    #print("\t Running ediva annotation with %d non biallelic variants:"%(len(not_biallelic_variants.items())))
    #for k, v in (not_biallelic_variants.items() ) :
    #    ediva=query_db(k,v,ediva,db,sep,missanndb,missanndb_coordinate,missanndbindel)
    #    ## end of while on not_biallelic_variants
    ### close DB connection
    db.close()
    return ediva


#@@ Check all input and output variables: evaluate if it's better to pass dictionaries
#@@ rather than variables.
## subroutine for Annovar annotation

def fill_annovar(FILE,Annovar, geneDef, comparison):
    '''Run Annovar for annotate variants and fill the data dictionary'''


    with open(FILE) as FILE_pointer:
        chr_idx = 26
        pos_idx = 27
        ref_idx = 29
        alt_idx = 30
        func_idx= 0
        gene_idx= 1
        exf_idx =2
        AAC_idx =3
        done  = []
        for line in FILE_pointer:
            if (line.startswith("Func")):
                #fields = line.strip().split(',')
                ##look for indexes positions
                #chr_idx = fields.index('Chr')
                #pos_idx = fields.index('Start')
                #ref_idx = fields.index('Ref')
                #alt_idx = fields.index('Obs')
                #func_idx= fields.index('Func')
                #gene_idx= fields.index('Gene')
                #exf_idx =fields.index('ExonicFunc')
                #AAC_idx =fields.index('AAChange')
                #print alt_idx
                #print '.....'
                pass
            else:
                newAnnovarLine = replaceCommainQoute(line)
                fields = newAnnovarLine.rstrip('\n').split(',')
                fields[alt_idx]=fields[alt_idx].replace('"','')
                if (fields[alt_idx].find(';')>=0):
                    annalts = fields[alt_idx].split(';')
                    fields[alt_idx] = annalts[0]
                chr_=fields[chr_idx].replace('chr','')
                chr_=chr_.replace('Chr','')
                valueTOmatch = chr_+";"+fields[pos_idx]+";"+fields[ref_idx]+";"+fields[alt_idx]
                valueTOmatch=valueTOmatch.replace('"','')
                ## fix missing values
                if len(fields[func_idx])==0:
                    fields[func_idx] = 'NA'
                if len(fields[gene_idx])==0:
                    fields[gene_idx] = 'NA'
                if len(fields[exf_idx])==0:
                    fields[exf_idx] = 'NA'
                if len(fields[AAC_idx])==0:
                    fields[AAC_idx] = 'NA'
                annToPass = fields[func_idx]+","+fields[gene_idx]+","+fields[exf_idx]+","+fields[AAC_idx]
                annToPass  =annToPass.replace('"','')
                if valueTOmatch in done:
                    pass
                    #print valueTOmatch
                else:
                    if geneDef == comparison:
                        Annovar[valueTOmatch] = annToPass   # Fill the Annovar Dictionary with values and keys
                    else:
                        Annovar[valueTOmatch] = Annovar[valueTOmatch]+ "," + annToPass
                    done.append(valueTOmatch)
        #       print "\n "
        #       print line
        #       print valueTOmatch
        #       print fields[0]
        #       print annToPass
        #       raise
    return Annovar


def AnnovarAnnotation(infile,templocation,fileSuffix,geneDef,ANNOVAR,Annovar,TABIX,MAF=0):
    perl ="/usr/bin/perl "+" "
    ## prepare Annovar input
    tabix_cmd = "PATH=$PATH:%s\n"%TABIX
    tmp = open(infile,'r')
    magic_number = tmp.read(2)
    zipped = magic_number=='\x1f\x8b'
    tmp.close()
    if zipped:
        infile_vcf=infile[:-3]+'.vcf'
        tabix_cmd+="bgzip -d %s -c > %s\n"%(infile,infile_vcf)
        remove_temp_cmd="rm %s\n"%(infile_vcf)
        print "MESSAGE :: NEED TO DECOMPRESS THE FILE , IT WILL TAKE A LITTLE LONGER"
    else:
        infile_vcf = infile
        remove_temp_cmd=''
    if MAF==0: #VCF
        annInCmm = perl + ANNOVAR+"/convert2annovar.pl --includeinfo -format vcf4 "+infile_vcf+" > "+templocation+"/annInfile"+fileSuffix+"   2> "+infile_vcf+".annovar.log\n"
    else:
        annInCmm = perl + ANNOVAR+"/maf2annovar.pl "+infile_vcf+" > "+templocation+"/annInfile"+fileSuffix+"   2> "+infile_vcf+".annovar.log\n"

    print "627 MESSAGE :: Running Annovar command \n > %s" %(tabix_cmd+annInCmm+remove_temp_cmd)
    subprocess.call(tabix_cmd+annInCmm+remove_temp_cmd,shell=True)
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
def header_defaults():
    sep=','
    head_common = ['Chr','Position','Reference','Alteration','QUAL','FILTER','AlleleFrequency' ]
    basic_head =['Chr','Position','Reference','Alteration']
    #head_common =basic_head.extend([ 'QUAL','FILTER','AlleleFrequency']  )
    common_fields = ['dbsnpIdentifier','EurEVSFrequency','AfrEVSFrequency','TotalEVSFrequency','Eur1000GenomesFrequency',
                 'Afr1000GenomesFrequency','Asia1000GenomesFrequency','Amr1000GenomesFrequency',
                 'Total1000GenomesFrequency','SegMentDup','PlacentalMammalPhyloP','PrimatesPhyloP',
                 'VertebratesPhyloP','PlacentalMammalPhastCons','PrimatesPhastCons',
                 'VertebratesPhastCons','Score1GERP++','Score2GERP++','SIFTScore',
                 'polyphen2','MutAss','Condel','Cadd1','Cadd2', 'Eigen_raw','Eigen_Phred','ABB_score',
                 'ExAC_AF','ExAC_adjusted_AF','ExAC_AFR','ExAC_AMR','ExAC_EAS',
                 'ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS','Clinvar Significance','Clinvar ID']
    repeat_fields=[     'SimpleTandemRepeatRegion','SimpleTandemRepeatLength']
    ensembl_annot=['Function(Ensembl)','Gene(Ensembl)','ExonicFunction(Ensembl)','AminoAcidChange(Ensembl)']
    refseq_annot =['Function(Refseq)','Gene(Refseq)','ExonicFunction(Refseq)','AminoAcidChange(Refseq)']
    known_annot=['Function(Known)','Gene(Known)','ExonicFunction(Known)','AminoAcidChange(Known)']
    annot=['Function','Gene','ExonicFunction','AminoAcidChange']
    return (sep,head_common,common_fields,repeat_fields,ensembl_annot,refseq_annot,known_annot,annot,basic_head)
## subroutnine por providing header to the main annotation output file
def getHeader(onlygenic,geneDef,headers,gtMode):
    stringTOreturn=""
    (sep,head_common,common_fields,repeat_fields,ensembl_annot,refseq_annot,known_annot,annot,basic_haed)=header_defaults()
    if (onlygenic):
            ## only genic annotation header
            ## check for gene definiton and construct header according to that
            # old removed :samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)
        if geneDef == 'ensGene':
            stringTOreturn = sep.join(head_common+ensembl_annot+ repeat_fields)
        elif geneDef == 'refGene':
            stringTOreturn = sep.join(head_common+refseq_annot+ repeat_fields)
        elif geneDef == 'knownGene':
            stringTOreturn = sep.join(head_common+known_annot+ repeat_fields)
        else:
            stringTOreturn = sep.join(head_common+refseq_annot+ensembl_annot+
                                       known_annot+ repeat_fields)
    else:
        if geneDef == 'ensGene':
            stringTOreturn = sep.join(head_common+ensembl_annot+common_fields+repeat_fields)+sep

        elif geneDef == 'refGene':
            stringTOreturn = sep.join(head_common+refseq_annot+common_fields+repeat_fields)

        elif geneDef == 'knownGene':
            stringTOreturn = sep.join(head_common+known_annot+common_fields+repeat_fields)
        else:
            stringTOreturn = sep.join(head_common+refseq_annot+ensembl_annot+
                                      known_annot+common_fields+repeat_fields)

    ## replace newlines with nothing at header line
    stringTOreturn=stringTOreturn.replace(sep,sep+'#')
    if len(headers)>9:
        #stringTOreturn+=sep+sep.join(headers[9:])
        if gtMode =='vcflike': sep_sample = ':'
        else: sep_sample=sep
        stringTOreturn+=sep+sep.join([sep_sample.join([x,'DP','REF','ALT','AF','GQ']) for x in headers[9:]])
    else:
        stringTOreturn+=sep+'samples'
    stringTOreturn.replace('\n','')
    stringTOreturn.replace(' ','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    stringTOreturn='#'+stringTOreturn
    return stringTOreturn



## subroutnine por providing header to the inconsistent annotation output file
def getHeaderIns(headers):
    stringTOreturn=""
    (sep,head_common,common_fields,repeat_fields,ensembl_annot,refseq_annot,known_annot,annot,basic_head)=header_defaults()

    stringTOreturn = sep.join(head_common+annot+common_fields+repeat_fields)
    stringTOreturn=stringTOreturn.replace(',',',#')
    if len(headers) > 9:
        stringTOreturn+=sep+sep.join(headers[9:len(headers)])
    else:
        stringTOreturn+=sep + 'samples'
    stringTOreturn.replace('\n','')
    stringTOreturn.replace(' ','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    stringTOreturn='#'+stringTOreturn
    return stringTOreturn


#@@ As before : match /usr/bin/perl output
## subroutnine por providing header to the quick look up mode annotation output file
def getHeaderQlookup(headers):
    stringTOreturn=""
    (sep,head_common,common_fields,repeat_fields,ensembl_annot,refseq_annot,known_annot,annot,basic_head)=header_defaults()
    #stringTOreturn = sep.join(head_common+annot+common_fields+repeat_fields)
    stringTOreturn = sep.join(basic_head+common_fields+repeat_fields)
    ## replace newlines with nothing at header line
    stringTOreturn.replace('\n','')
    #stringTOreturn.replace('\s+','')
    #=~ s/\n|\s+//g;  verify that this is correct!!
    return stringTOreturn

##############################################################################################
## VCF PROCESSING
##############################################################################################
def process_line(myline,i,adindex,gtindex,dpindex,gqindex):
    ''' process vcf input file line by line'''
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
        if gtindex != "NF":     genotype = gts[gtindex]
        else:                   genotype = "."
        if dpindex != "NF" and dpindex<len(gts):
            DP = gts[dpindex]
        else :
            DP='0'
        try:
            if gqindex != "NF" and len(gts)>gqindex:GQ = gts[gqindex]
            else :                      GQ='0'
        except:
            print gqindex
            print dpindex
            print gts
            raise

    else:
        genotype = myline[i]
        dpref = "."
        dpalt = "."
        samAf = "."
        DP='0'
        GQ='0'
    return (genotype,DP,dpref,dpalt,samAf,GQ)

def samples_fill(samples,gtMode,s_key,s_val,genotype):
    if gtMode == "complete":
        if (samples.get(s_key,False)):
            samples[s_key] += ","+s_val
        else:
            samples[s_key] = s_val
    elif gtMode == "vcflike":
        if (samples.get(s_key,False)):
            samples[s_key] += ","+s_val.replace(',',':')
        else:
            samples[s_key] = s_val.replace(',',':')
    else: ## compact
        ## kick out genotypes of '0/0','./.','0|0' and '.|.'
        if genotype != '0/0' and genotype != './.' and genotype != '.|.' and genotype != '0|0':
            if (samples.get(s_key,False)):
                samples[s_key] += ","+s_val
            else:
                samples[s_key] = s_val
    return samples


def check_for_genotype(chr_col,position,token_ref,token_obs,sample_info,gtMode,samples):
    ''' checks for sample genotype and adds it to the sample array , processed '''
    (genotype,DP,dpref,dpalt,samAf,GQ) = sample_info
    s_key =chr_col+';'+position+';'+token_ref+';'+token_obs
    s_val = ','.join([genotype,DP,dpref,dpalt,samAf,GQ])
    samples = samples_fill(samples,gtMode,s_key,s_val,genotype)

    return samples

def process_fileline(infile,myline,ref,al,alfr,variants,not_biallelic_variants,samples,type_in,gtMode,j):
    ## decide for variant type
    chr_col     = myline[0]
    if chr_col.startswith('chr') or chr_col.startswith('Chr'):chr_col=chr_col[3:]
    position    = myline[1]
    qual        = myline[5]
    filter_         = myline[6].replace(';','|')
    gtindex     = "NF"
    adindex     = "NF"
    gqindex     = "NF"
    dpindex     = "NF"
    lenref = len(ref)
    lenalt = len(al)

    keychars = list()#['n','N']
    if gtMode != "none":
        if len(myline) > 8 and ':' in myline[8] :
            gtcheck = myline[8].split(':')
            for gti in range(0,len(gtcheck)):
                if gtcheck[gti] == "GT":
                    gtindex = gti
                if gtcheck[gti] == "AD":
                    adindex = gti
                if gtcheck[gti] == "DP":
                    dpindex = gti
                if gtcheck[gti] == "GQ":
                    gqindex = gti
            # check for the GT and AD fields
            if gtindex == "NF":
                print ("WARNING:: Weird genotype format %s found in $input at chromosome %s and position %s. No GT field present in %s \n"
                       %(myline[8],chr_col,position,myline[8] ))
            if adindex == "NF":
                print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. No AD field present in %s \n"
                %(myline[8], infile,chr_col,position,myline[8]))
            if dpindex == "NF":
                print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. No DP field present in %s \n"
                %(myline[8], infile,chr_col,position,myline[8]))
            if gqindex == "NF":
                print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. No GQ field present in %s \n"
                %(myline[8], infile,chr_col,position,myline[8]))
        else:
            pass
            #print ("WARNING:: Weird genotype format %s found in %s at chromosome %s and position %s. Expected genotype format field separator \":\" \n"
            #%(myline[8],infile,chr_col,position))
## process based on alteration
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
                variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+al+';'+alfr+';'+qual+';'+filter_
            else:
                not_biallelic_variants[ chr_col+';'+position+';'+token_ref+';'+token_obs] = chr_col+';'+position+';'+ref+';'+al+';'+alfr+';'+qual+';'+filter_
             ## if sample wise information is present in the VCF then process ; otherwise skip
             ## also check for none value in genotype mode parameter
            if (len(myline) > 8) and (gtMode != "none"):
                s_key =chr_col+';'+position+';'+token_ref+';'+token_obs
                if (samples.get(s_key,False)):
                    samples[s_key]=""

                for i in range(9,len(myline)):
                    #(genotype,dpref,dpalt,samAf,GQ) = process_line(myline,i,adindex,gtindex)
                    sample_info = process_line(myline,i,adindex,gtindex,dpindex,gqindex)
                    ## check for sample genotype mode
                    samples = check_for_genotype(chr_col,position,token_ref,token_obs,sample_info,gtMode,samples)
                    #s_key =chr_col+';'+position+';'+token_ref+';'+token_obs
                    #s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
                    #samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
    else: ## SNP
        ## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
        ## e.g A,C in the alternate column in VCF will report only A in the main annotation file
        ## we are doing this because we want to keep the annotation main file consistent
        #print 'snp : %s ' % type_in
        if type_in == 'SNP' or type_in == 'all':
            ## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
            if j == 0 and  not(any(k in ref for k in keychars)) and not(any(k in al for k in keychars) ):
                variants[ chr_col+';'+position+';'+ref+';'+al] = chr_col+';'+position+';'+ref+';'+al+';'+alfr+';'+qual+';'+filter_
            else:
                not_biallelic_variants[ chr_col+';'+position+';'+ref+';'+al] = chr_col+';'+position+';'+ref+';'+al+';'+alfr+';'+qual+';'+filter_
        ## if sample wise information is present in the VCF then process ; otherwise skip
        ## also check for none value in genotype mode parameter
                s_key =chr_col+';'+position+';'+ref+';'+al
                if (samples.get(s_key,False)):
                    samples[s_key]=""
            if len(myline) > 8 and gtMode != "none":
                for i in range(9,len(myline)):
                    sample_info = process_line(myline,i,adindex,gtindex,dpindex,gqindex)
                    ## check for sample genotype mode
                    samples = check_for_genotype(chr_col,position,ref,al,sample_info,gtMode,samples)
                    #s_key =chr_col+';'+position+';'+ref+';'+al
                    #s_val = headers[i]+">"+genotype+">"+dpref+">"+dpalt+">"+samAf
                    #samples = samples_fill(samples,gtMode,s_key,s_val,genotype)
    return(variants,not_biallelic_variants,samples)



#@@ Processes the vcf file and extract information about biallelic variants, non biallelic variants and samples
def vcf_processing(infile,qlookup,gtMode,type_in):
    ''' vcf file processing '''
    allowed_chr = list()
    for i in range(23):
        allowed_chr.append(str(i+1))
    allowed_chr+=(['X','Y','x','y'])
    skipped_chr = list()
    db = 'temp.sqlite'#/tmp/test.sqlite'
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
        print 'MESSAGE :: BGZIP file  : %s'%(magic_number=='\x1f\x8b')
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
                    myline = line.rstrip('\n').rstrip('\t').split('\t')
                    ## check for malformed VCF file and take action
                    if (len(myline) < 8):
                        print "ERROR:: Not a valid VCF format \n 866"
                        raise IOError
                    chr_col     = myline[0]
                    if chr_col.startswith('chr') or chr_col.startswith('Chr'):chr_col=chr_col[3:]
                    position    = myline[1]
                    ref         = myline[3]
                    alt         = myline[4]
                    qual        = myline[5]
                    filter_     = myline[6]
                    alt_tmp = alt.split(',')
                    while '<NON_REF>' in alt_tmp: alt_tmp.remove('<NON_REF>')
                    alt = ','.join(alt_tmp)
                    infos       = myline[7].split(';')
                    AF          = "."
                    gtindex     = "NF"
                    adindex     = "NF"
                    gqindex     = "NF"
                    dpindex     = "NF"
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
                    keywords_alt = ['A','T','G','C','-',',','n','N']

                    if not chr_col in allowed_chr:
                        skipped_chr.append(chr_col)
                    #if any(k in chr_col for k in keywords):
                    #        print( "WARNING:: !! does not support chromosome %s currently. This variant line will be skipped in the final annotation output " % (chr_col))
                    #elif not(any( k in alt for k in keywords_alt)):
                    elif not(all(k in keywords_alt for k in alt)):
                        print( "WARNING:: Unknown alternate allele %s detected at %s and %s. This variant line will be skipped in the final annotation output " %(alt, chr_col,position))
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
#                        ## process based on alteration
                        if ',' in alt : ## section for all sites other than bi-allelic
                            alts = alt.split(',')
                            afs = AF.split(',')
                            if len(afs)<len(alts):
                                afs=['.']*len(alts)

                            ## process each alteration allele
                            for j in range(0,1):##,len(alts)):
                                al = alts[j]
                                alfr = afs[j]
                                (variants,not_biallelic_variants,samples) = process_fileline(infile,myline,ref,al,alfr,variants,not_biallelic_variants,samples,type_in,gtMode,j)
                        else: ## section for bi-allelic sites
                            (variants,not_biallelic_variants,samples) = process_fileline(infile,myline,ref,alt,AF,variants,not_biallelic_variants,samples,type_in,gtMode,0)
                #if len(variants) >= 100000:
                #    variants = dump_to_db(variants,db,'variants')
                #if len(not_biallelic_variants) >= 100000:
                #    not_biallelic_variants = dump_to_db(not_biallelic_variants,db,'notbvariants')
                #if len(samples) >= 100000:
                #    samples = dump_to_db(samples,db,'samples')


            #Reporting skipped lines
            if len(skipped_chr)>0:
                tmp =  dict((i,skipped_chr.count(i)) for i in set(skipped_chr))
                print "WARNING: there have been skipped lines due to unknown chr:"
                for i in list(set(skipped_chr)):
                    print "\tChr: %s - %d lines"%(i,tmp.get(i))
        #    if len(variants) > 0:
        #       variants = dump_to_db(variants,db,'variants')
        #    if len(not_biallelic_variants) > 0:
        #       not_biallelic_variants = dump_to_db(not_biallelic_variants,db,'notbvariants')
        #    if len(samples) > 0:
        #       samples = dump_to_db(samples,db,'samples')
    else :
        ## Quick lookup mode
        ## check for input type in -q parameter
        ## decide for file or query string
        if os.path.isfile(qlookup):
            ## its a file
            # print "MESSAGE :: Processing input file - %s " %qlookup
            ##check if it is a MAF file or a quick lookup file
            MAF = 0
            with  open(qlookup) as rd:
                for line in rd:
                    tmp = line.strip()
                    if not(line.startswith('#')):
                        if line.count(':') != 3 and line.count('\t')>1:
                            print 'Reading file as a MAF file'
                            MAF=1
                        break

            counter = 0
            with open(qlookup) as FL:
                for line in FL:
                    if MAF ==0:
                        var = line.rstrip('\n').rstrip('\t').split(':')
                    else:
                        if counter==0  :
                            if not(line.startswith('#')):
                                counter +=1
                            continue
                        fields = line.strip().split('\t')
                        if fields[11] == fields[10]:
                            var = [fields[4],fields[5],fields[10],fields[12]]
                        else:
                            var = [fields[4],fields[5],fields[10],fields[11]]
                        # Chr Pos Refernce 2nd Tumor allele
                    counter +=1
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
#                               print( "WARNING:: !! does not support chromosome %s currently. This variant line will be skipped in the final annotation output " % (chr_col))
                            if not chr_col in allowed_chr:
                                skipped_chr.append(chr_col)
                                #elif not(any( k in alt for k in keywords_alt)):
                            elif not(all(k in keywords_alt for k in alt)):
                                print( "WARNING:: Unknown alternate allele detected at %s and %s. This variant line will be skipped in the final annotation output " %( chr_col,pos))
                                print 'line : 1151'
                                print alt
                                #raise
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

    return samples,variants,not_biallelic_variants,headers
