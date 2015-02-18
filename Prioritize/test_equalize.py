#!/usr/bin/env python
import scipy as sp
import scipy.stats
import pprint
import argparse
import csv
import re
import os
import MySQLdb
try:
    import xlsxwriter
    import xlrd
    writeXLS = True
    print 'yes XLS module imported'
except:
    writeXLS = False
    print 'No XLS module imported'
    
### Note header
#[   'Chr',
#    'Position',
#    'Reference',
#    'Alteration',
#    'Function(Refseq)',
#    'Gene(Refseq)',
#    'ExonicFunction(Refseq)',
#    'AminoAcidChange(Refseq)',
#    'Function(Ensembl)',
#    'Gene(Ensembl)',
#    'ExonicFunction(Ensembl)',
#    'AminoAcidChange(Ensembl)',
#    'Function(Known)',
#    'Gene(Known)',
#    'ExonicFunction(Known)',
#    'AminoAcidChange(Known)',
#    'dbsnpIdentifier',
#    'dbSNPfrequency',
#    'EurEVSFrequency',
#    'AfrEVSFrequency',
#    'TotalEVSFrequency',
#    'Eur1000GenomesFrequency',
#    'Afr1000GenomesFrequency',
#    'Amr1000GenomesFrequency',
#    'Asia1000GenomesFrequency',
#    'Total1000GenomesFrequency',
#    'SugMentDup',
#    'PlacentalMammalPhyloP',
#    'PrimatesPhyloP',
#    'VertebratesPhyloP',
#    'PlacentalMammalPhastCons',
#    'PrimatesPhastCons',
#    'VertebratesPhastCons',
#    'Score1GERP++',
#    'Score2GERP++',
#    'SIFTScore',
#    'polyphen2',
#    'MutAss',
#    'Condel',
#    'samples(sampleid>zygosity>Cov>AF)']
###


def main ():
    
    sub_pp = pprint.PrettyPrinter(indent = 4)    
    parser = argparse.ArgumentParser(description = 'rank SNPs according to their mutation properties')    
    parser.add_argument('--infile', type=argparse.FileType('r'), dest='infile', required=True, help='comma separated list of SNPs annotated with mutation impact data. [required]')
    parser.add_argument('--outfile', type=argparse.FileType('w+'), dest='outfile', required=True, help='comma separated list of SNPs annotated with ranks. [required]')
    args = parser.parse_args()
    alldata = list(csv.reader(args.infile))
    
    header = alldata.pop(0)
    header_dict = dict(zip(header,range(len(header))))
    index_varfunction = header_dict.get('ExonicFunction(Refseq)','NA')#identifycolumns(header, 'ExonicFunction(Refseq)')
    index_genicfunction = header_dict.get('Function(Refseq)','NA')#identifycolumns(header, 'Function(Refseq)')
    index_allele_freq  = header_dict.get('AlleleFrequency','NA')
    
    
    #alldata_clean = [ line for line in alldata if not line[index_varfunction] == 'synonymous SNV' ]
    alldata_transpose = zip(*alldata)       
    values2investigate = ['Total1000GenomesFrequency', 'TotalEVSFrequency', 'SegMentDup', 'Condel', 'VertebratesPhyloP', 'VertebratesPhastCons', 'Cadd2', 'SIFTScore']
    
    ## This one is to produce a rank-filter by Allele Frequency with a sample tool
    while True: # Allele frequency filtering
        print('The allele Frequency filtering with tukey window for data input between 0-1')
        AF_data = alldata_transpose[index_allele_freq]
        #AF_data = sp.linspace(0,1-0.05,20).tolist()
        AF_data = [0,00.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00]
        AF_filter = [tukeywin(x,0.85) for x in AF_data]
        #print AF_data
        sub_pp.pprint( zip(AF_data ,AF_filter))
        break
    
    
    # binning, because otherwise subtle differences get too much weight
    print('binning')
    while True:
        #binned_values = binning(alldata_transpose, header, 'MAF')    
        _1000G        =sp.array([0 if x == 'NA' else x for x in alldata_transpose[header_dict.get(values2investigate[0])]],dtype=sp.float32)
        _EVS_Freq     =sp.array([0 if x == 'NA' else x for x in alldata_transpose[header_dict.get(values2investigate[1])]],dtype=sp.float32)
        MAF           = (0.5*_1000G+0.5*_EVS_Freq).tolist()
        binned_values = bin_list(MAF)
        ranked_maf    = (binned_values)
        #print binned_values
        
        #binned_values = binning(alldata_transpose, header, 'segdup')
        binned_values = bin_list(alldata_transpose[header_dict.get(values2investigate[2])])
        ranked_segdup = (binned_values)
        #print binned_values
        
        #binned_values = binning(alldata_transpose, header, 'condel')
        binned_values = bin_list(alldata_transpose[header_dict.get(values2investigate[3])],None,None,True)
        ranked_condel = (binned_values)
        #print binned_values
        
        #binned_values = binning(alldata_transpose, header, 'PhyloP')
        binned_values = bin_list(alldata_transpose[header_dict.get(values2investigate[4])],None,None,True)
        ranked_phylop = (binned_values)
    
        #print binned_values
        
        #binned_values = binning(alldata_transpose, header, 'PhastCons')
        binned_values = bin_list(alldata_transpose[header_dict.get(values2investigate[5])])
        ranked_phastcons = (binned_values)
    
        #print binned_values
        
        #binned_values = binning(alldata_transpose, header, 'Cadd2')
        binned_values = bin_list(alldata_transpose[header_dict.get(values2investigate[6])])
        ranked_cadd   = binned_values
    
        #print binned_values
        
        #sub_pp.pprint(binned_values)
        #exit(0)
        break
    print('Rank product ')
    # calculate rank product
    rank_product_list = list()    
    max_maf_rank       = max(ranked_maf)
    max_segdup_rank    = max(ranked_segdup)
    max_condel_rank    = max(ranked_condel)
    max_phastcons_rank = max(ranked_phastcons)
    max_cadd_rank      = max(ranked_cadd)
    div_factor = 100**4
    
    for i in range( len(binned_values) ):
        # skip synonymous variants
        if ( alldata[i][index_varfunction] == 'synonymous SNV' or alldata[i][index_varfunction] == 'NA' ) and not alldata[i][index_genicfunction] == 'splicing':
            # synonymous SNVs get the maximum rank and are downgraded by that
            rank_product = float( ( max_maf_rank * max_segdup_rank * max_condel_rank * max_phastcons_rank * max_cadd_rank ) ) / ( div_factor)
        else:            
            rank_product = float( ranked_maf[i] * ranked_segdup[i] * ranked_condel[i] * ranked_phastcons[i] * ranked_cadd[i] ) / ( div_factor ) # 4 tools deliver information, decrease the numeric value to more graspable values ### currently deleted * ranked_phylop[i]
        
        rank_product_list.append(rank_product)
    
    # all rank products get a rank for more facile overview
    rankrank = scipy.stats.rankdata(rank_product_list)
    rank_idx = sp.argsort(rank_product_list)

    
    outcsv = csv.writer(args.outfile)
    header.append('rank') #,'maf','segdup','condel','condelbin','phastcons','product'])
    outcsv.writerow(header)
    for i in rank_idx:
        tmp = alldata[i]
        tmp.append(int(sp.ceil(100*rank_product_list[i])))
        outcsv.writerow(tmp)
    ##############################################
    # MEDIAN FILTERING
    all_ranked_array = sp.matrix([ranked_maf,ranked_segdup,ranked_condel,ranked_cadd,ranked_phastcons])
    mean_rank = sp.mean(all_ranked_array,axis=0)
    #sub_pp.pprint(mean_rank)
    median_rank = [sp.median(all_ranked_array[:,i].T.tolist()) for i in range(len(ranked_maf))]

    
    
    ##############################################
    print('Database search for OMIM values for a subset of variants')
    db_search(alldata_transpose[6],alldata_transpose[0],alldata_transpose[1],alldata_transpose[2],alldata_transpose[3])
    
    #excel writing
    ### write an xls output by reading the existing output file if it exists
    ### so that we can then add a new worksheet
    
    print('Excel writing')
    if writeXLS == True:
        # open output file for re-reading
        excel_name = args.outfile.name + ".xlsx"
        tmp_name = 'tmp.xlsx'
        
        # open xls file for writing
        xls = xlsxwriter.Workbook(tmp_name)
        worksheet = xls.add_worksheet('ediva_filtered3')
        row = 0
        args.outfile.seek(0)
        # read line by line and transform to xls
        for line in args.outfile:
            #line.rstrip('\n')
            data = line.split(',')
            worksheet.write_row(row, 0, data)
            
            row += 1
           
            
        #check if already exist
        with open(excel_name,'r') as old_excel:
            workbook_rd = xlrd.open_workbook(excel_name)
            worksheets = workbook_rd.sheet_names()
            for worksheet_name in worksheets:
                try:
                    worksheet_rd = workbook_rd.sheet_by_name(worksheet_name)
                    worksheet_wr = xls.add_worksheet(worksheet_name)
                    
                    num_rows = worksheet_rd.nrows - 1
                    curr_row = -1
                    while curr_row < num_rows:
                            curr_row += 1
                            row_content = worksheet_rd.row(curr_row)
                            row_content_values = [x.value for x in row_content]
                            worksheet_wr.write_row(curr_row, 0, row_content_values)
                            #print row
                except:
                    print "There was a problem in passing the %s data. \nIt may be because the sheet was already there before"%(worksheet_name)
        xls.close()
        os.remove(excel_name)
        os.rename(tmp_name, excel_name)
    
    return #exit(0)

def rank (binned_list):
    # rank returns values from 1 to 100
    rank_list = list()
    for i in range( len(binned_list) ):
        rankvalue = binned_list[i] % 101 + 1 # 101, because otherwise 0 and 100 % 100 calculates to 0
        rank_list.append(rankvalue)
    
    #rank_list = scipy.stats.rankdata(binned_list) # this does produce the wrong order
    
    return(rank_list)

def bin_list(in_list,th_1=None,th_2=None,invert=False,n_bins=100):
    """ I want to build a function that takes in:
        - list : of floats and takes care of 'NA'
        - th_1 : small threshold: all values <= th_1 are assigned value 1
        - th_2 : high threshold : all values >= th_2 are assigned value 100
        - invert:decide if the input list must be inverted in terms of output values 1-100 -> 100-1
        - n_bins =100: number of total available scores
        
        - The output must be an associated integer value 1-100 for each input list_element
    """
    sub_pp = pprint.PrettyPrinter(indent = 4)       
    #First of all take care of NA values by the mean array_value
    n_bins_0 = n_bins
    offset = 1

    in_array = sp.array(in_list)
    in_array[in_array=='NA'] = sp.nan
    in_array = in_array.astype(sp.float32)
    in_mean  = scipy.stats.nanmean(in_array)
    nan_idx  = sp.where(sp.isnan(in_array))
    in_array[nan_idx] = in_mean    
    out_array = sp.zeros(len(in_list))

    
    #Generate a sorted_list of the values between the thresholds:
    temp_array = in_array
    if th_1 != None:
        n_bins -= 1
        offset += 1
        th_1_indexes = sp.where(in_array<= th_1)
        out_array[th_1_indexes] = 1
        temp_array[th_1_indexes] = sp.nan
    if th_2 != None:
        n_bins -= 1
        th_2_indexes = sp.where(in_array>=th_2)
        out_array[th_2_indexes] = n_bins
        temp_array[th_2_indexes] = sp.nan
    
    #Establish the number of element per bin
    #temp_array.astype(sp.float32)
    
    num_elements_within_threshold = len(temp_array[~sp.isnan(temp_array)])

    #important check if the elements are more than n_bins
    if num_elements_within_threshold<n_bins:
        n_bins = num_elements_within_threshold
        print "WARNING: the elements within thresholds are less than the number of bins"
        print "WARNING: Set the within_threshold bin number to %d"%n_bins
    num_el_per_bin = int(sp.floor(num_elements_within_threshold/float(n_bins)))
    #print '\n num_within:%d'%num_elements_within_threshold
    #print 'num bins  :%d'%n_bins
    #print 'n_per_bin    :%d'%num_el_per_bin
    #print 'offset  : %d'%offset 
    sort_indices = sp.argsort(temp_array)
    sort_indices = sort_indices[0:num_elements_within_threshold]
    sorted_val   = sp.sort(temp_array)
    #build the max_value_per_bin
    max_per_bin = sp.zeros(n_bins)
    for i in range(n_bins):
        index = int(sp.floor((i+1)*num_elements_within_threshold/float(n_bins) )) 
        max_per_bin[i] = sorted_val[index-1]

    for i in range(len(temp_array)):
        if ~sp.isnan(temp_array[i]):
            for bin_i in range(n_bins):
                if temp_array[i] <= max_per_bin[bin_i]:
                    out_array[i] = offset + bin_i
                    break
            
    
    
    
    #bin_values   = offset+sp.floor(sp.linspace(0,num_elements_within_threshold-1,num_elements_within_threshold)/float(num_el_per_bin))
    #out_array  = out_array[sort_indices.astype(int).tolist()] = bin_values
    

    #Manage the case of inverted values
    if invert:
        out_array = n_bins_0 + 1 - out_array
    out_array = out_array.astype(int)
    out_array = [min(n_bins_0,x)for x in out_array]
    #aa = [out_array.count(x) for x in range(1,n_bins_0+1)]
    ##sub_pp.pprint(out_array)
    #dd = dict(zip(range(1,n_bins_0+1),aa))
    #sub_pp.pprint(dd)
    #print(sum(aa))
    #print '----------------------------\n\n'
    return out_array#.astype(int).tolist()
  
  
def db_search(variant_list,chr_,pos,ref,alt):
    ''' Function to search in the edivaweb database the omim links for each of the ranked variants
    And then store the result in a new table.
    It's a test to see wether we can get the data from that database or not'''
    
    
    outStr = dict()
    ## DB parameters
    username    = "rrahman"
    database    = "eDiVa_scratch"
    dbhost      = "mysqlsrv-ediva.linux.crg.es"
    passw       = "mDWr41J1"
     
    db = MySQLdb.connect(host=dbhost, # your host, usually localhost
    user=username, # your username
    passwd=passw)#, # your password
    #db=database) # name of the data base
    
    cur = db.cursor()
    for i in range(10):#len(variant_list)):
        gene_name = variant_list[i]
        
        sql = ("SELECT gene_name , title_mim_number ,details_mim_number "+
               "FROM eDiVa_scratch.Table_gene2mim, eDiVa_scratch.Table_omim "+
               "where eDiVa_scratch.Table_gene2mim.mim_number = eDiVa_scratch.Table_omim.mim_number "+
               "and eDiVa_scratch.Table_gene2mim.gene_name ='%s';"%gene_name)
            #sql = "select chr,pos,lengthofrepeat,copyNum,region from ediva_public_omics.Table_simpleRepeat;"
    
        cur.execute(sql)

        count = 0
        omim_disease =""
        omim_2 = ""
        for row in cur:
            count +=1
            print row
   
            omim_2+=row[2]+ " | "
            omim_disease+=row[1]+ " | "
        print omim_disease
        print omim_2
        print '\n\n'
        if count >1:
            print("%s  : %d mim_terms"%(gene_name,count))
        
        
        sql_clinvar =("SELECT clinical_significance, disease_name, clinical_review, access_number "+
                      "FROM eDiVa_scratch.Table_clinvar "+
                      "WHERE chr='%s' and pos='%s' and ref='%s'and alt='%s'"
                      %(chr_[i],pos[i],ref[i],alt[i] )
        )
        cur.execute(sql_clinvar)

        count = 0
        for row in cur:
            count +=1
            #print row[0:2]
        if count >1:
            print("%s  : %d clinvar_count \n"%(gene_name,count))
    cur.close()
    db.close()
    
def tukeywin(x, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.
 
    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output
 
    Reference
    ---------
 
http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    '''
    # Special cases
    if alpha <= 0:
        return 1 #rectangular window
    elif alpha >= 1:
        return 1
 
    # Normal case
    #x = sp.linspace(0, 1, window_length)
    #w = sp.ones(x.shape)
    w = 1
    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    if first_condition:
        w   = 0.5 * (1 + sp.cos(2*sp.pi/alpha * (2*(x - alpha/2) )))
 
    # second condition already taken care of
 
    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>(1 - alpha/2)
    if third_condition:
        w    = 0.5 * (1 + sp.cos(2*sp.pi/alpha * (2*(x - 1 + alpha/2))) )
 
    return round(w    ,2)
    
main()
