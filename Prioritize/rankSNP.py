import scipy as sp
import scipy.stats
import pprint
import argparse
import csv
import re

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
    parser.add_argument('--outfile', type=argparse.FileType('w'), dest='outfile', required=True, help='comma separated list of SNPs annotated with ranks. [required]')
    
    args = parser.parse_args()
    
    alldata = list(csv.reader(args.infile))

    header = alldata.pop(0)
    index_varfunction = identifycolumns(header, 'ExonicFunction(Refseq)')
    index_genicfunction = identifycolumns(header, 'Function(Refseq)')
    
    #alldata_clean = [ line for line in alldata if not line[index_varfunction] == 'synonymous SNV' ]
    print header 
    alldata_transpose = zip(*alldata)
    
    
    values2investigate = ['Total1000GenomesFrequency', 'TotalEVSFrequency', 'SegMentDup', 'Condel', 'VertebratesPhyloP', 'VertebratesPhastCons', 'SIFTScore']
    
    # binning, because otherwise subtle differences get too much weight
    binned_values = binning(alldata_transpose, header, 'MAF')
    ranked_maf    = rank(binned_values)
    
    binned_values = binning(alldata_transpose, header, 'segdup')
    ranked_segdup = rank(binned_values)
    
    binned_values = binning(alldata_transpose, header, 'condel')
    ranked_condel = rank(binned_values)
    
    binned_values = binning(alldata_transpose, header, 'PhyloP')
    ranked_phylop = rank(binned_values)
    
    binned_values = binning(alldata_transpose, header, 'PhastCons')
    ranked_phastcons = rank(binned_values)
    
    binned_values = binning(alldata_transpose, header, 'Cadd2')
    ranked_cadd   = binned_values
    
    #sub_pp.pprint(binned_values)
    #exit(0)
    
    # calculate rank product
    rank_product_list = list()
    
    max_maf_rank       = max(ranked_maf)
    max_segdup_rank    = max(ranked_segdup)
    max_condel_rank    = max(ranked_condel)
    max_phastcons_rank = max(ranked_phastcons)
    max_cadd_rank      = max(ranked_cadd)
    
    
    for i in range( len(binned_values) ):
        
        
        # skip synonymous variants
        if ( alldata[i][index_varfunction] == 'synonymous SNV' or alldata[i][index_varfunction] == 'NA' ) and not alldata[i][index_genicfunction] == 'splicing':
            # synonymous SNVs get the maximum rank and are downgraded by that
            rank_product = float( ( max_maf_rank * max_segdup_rank * max_condel_rank * max_phastcons_rank * max_cadd_rank ) ) / ( 100**2 )
        else:
            
            rank_product = float( ranked_maf[i] * ranked_segdup[i] * ranked_condel[i] * ranked_phastcons[i] * ranked_cadd[i] ) / ( 100**2 ) # 4 tools deliver information, decrease the numeric value to more graspable values ### currently deleted * ranked_phylop[i]
        
        rank_product_list.append(rank_product)
    
    # all rank products get a rank for more facile overview
    rankrank = scipy.stats.rankdata(rank_product_list)
    
    for i in range( len(alldata) ):
        # bin the final ranks into groups, to decrease resolution and have a max of about 1000
        binned_rank = int(rankrank[i] / 100)
        #alldata[i].append( binned_rank )
        alldata[i].append( rank_product_list[i] )
        #alldata[i].append(int(rankrank[i])) #, ranked_maf[i], ranked_segdup[i], ranked_condel[i], debug_saver[i], ranked_phastcons[i], rank_product_list[i]])
    
    outcsv = csv.writer(args.outfile)
    header.append('rank') #,'maf','segdup','condel','condelbin','phastcons','product'])
    outcsv.writerow(header)
    for line in alldata:
        outcsv.writerow(line)
    
    exit(0)

###########################
###### subroutines ########
###########################

def binning (alldata, header, parameter):
    
    #sub_pp = pprint.PrettyPrinter(indent = 6) #DEBUG
    binned_values = list()
    
    # MAF
    if parameter == 'MAF':
        index_1000G  = identifycolumns(header, 'Total1000GenomesFrequency')
        index_EVS    = identifycolumns(header, 'TotalEVSFrequency')
        column_1000G = list(alldata[index_1000G])
        column_EVS   = list(alldata[index_EVS])

        
        # clean out closet (substitute NA by 0):
        for i in range(len(column_1000G)):
            if column_1000G[i] == 'NA':
                column_1000G[i] = str(0)
            if column_EVS[i] == 'NA':
                column_EVS[i]   = str(0)


        array_1000G = sp.array(column_1000G, dtype=sp.float32)
        mean_1000G  = sp.mean(array_1000G)
        
        array_EVS = sp.array(column_EVS, dtype=sp.float32)
        mean_EVS  = sp.mean(array_EVS)
        
        
        for i in range( len(column_1000G) ):
            MAF = max( round(float(column_1000G[i]), 2), round(float(column_EVS[i]), 2) )
            bin_value = int(MAF * 100) # (good 1-100 bad) small values should get small ranks
            binned_values.append(bin_value)
    
    # Segmental duplications
    elif parameter == 'segdup':
        index_segdup  = identifycolumns(header, 'SegMentDup')
        column_segdup = alldata[index_segdup]
                
        mean_segdup = getmean(column_segdup)
        
        for i in range( len(column_segdup) ):
            if column_segdup[i] == 'NA':
                segdup = round(mean_segdup, 2)
                bin_value = int(segdup * 100)  
                binned_values.append(bin_value) # (good 1-100 bad)
            else:
                segdup = round(float(column_segdup[i]), 2)
                bin_value = int(segdup * 100) # small values get small ranks
                binned_values.append(bin_value) # (good 1-100 bad)
    
    # Condel
    elif parameter == 'condel':
        index_condel  = identifycolumns(header, 'Condel')
        column_condel = alldata[index_condel]
        
        index_function  = identifycolumns(header, 'Function(Refseq)')
        column_function = alldata[index_function]
        
        mean_condel = getmean(column_condel)

        for i in range( len(column_condel) ):
            if column_condel[i] == 'NA' and not column_function[i] == 'splicing':
                condel = round(mean_condel, 2)
            
            elif column_condel[i] == 'NA' and column_function[i] == 'splicing':
                condel = 1
                
            else:
                condel = round(float(column_condel[i]), 2)
                
            
            bin_value = int(condel * 100) # binning of condel values (bad 1-100 good)
            bin_value_rev = abs(bin_value - 100 ) # reverse order
            binned_values.append(bin_value_rev)
    
    # PhyloP --- this value is not used in the rank product!!!
    elif parameter == 'PhyloP':
        index_phylop  = identifycolumns(header, 'VertebratesPhyloP')
        column_phylop = alldata[index_phylop]
        
        mean_phylop = getmean(column_phylop)
        
        for i in range( len(column_phylop) ):
            if column_phylop[i] == 'NA':
                #NAcount += 1 # DEBUG
                #print 'found NA in phylop' #DEBUG
                phylop = round(mean_phylop, 2)
                bin_value = int(phylop * 100) # [bad 1 - 100 good]
                bin_value_rev = abs(bin_value - 100 ) # reverse order
                binned_values.append(bin_value_rev)
            else:
                phylop = round(float(column_phylop[i]), 2)
                bin_value = int(phylop * 100)  # [bad 1 - 100 good]
                bin_value_rev = abs(bin_value - 100 ) # reverse order
                binned_values.append(bin_value_rev)
    
    # PhastCons
    # The phastCons scores, by contrast, represent probabilities of negative selection and range between 0 and 1.
    # --- PhastCons = 1, could be a disease variant.
    elif parameter == 'PhastCons':
        index_phastcons  = identifycolumns(header, 'VertebratesPhastCons')
        column_phastcons = alldata[index_phastcons]
        
        mean_phastcons = getmean(column_phastcons)
        
        for i in range( len(column_phastcons) ):
            if column_phastcons[i] == 'NA':
                phastcons = round(mean_phastcons, 2)
                bin_value = int(phastcons * 100)
                bin_value_rev = abs(bin_value - 100 )
                binned_values.append(bin_value_rev) # [good 1 - 100 bad]
            else:
                phastcons = round(float(column_phastcons[i]), 2)
                if phastcons == 0: phastcons = 0.001
                bin_value = int(phastcons * 100)  # [good 1 - 100 bad] a value of 1 should get a small rank
                bin_value_rev = abs(bin_value - 100 )
                binned_values.append(bin_value_rev)
    
    # Cadd
    # Cadd already represents ranked scores ranging from
    elif parameter == 'Cadd2':
        index_cadd = identifycolumns(header, 'Cadd2')
        column_cadd = alldata[index_cadd]
        
        mean_cadd = getmean(column_cadd)
        
        for i in range( len(column_cadd) ):
            if column_cadd[i] == 'NA':
                cadd = round(mean_cadd, 2)
                
                bin_value_rev = abs(cadd - 100 )
                binned_values.append(bin_value_rev) # [good 1 - 100 bad]
            else:
                
                bin_value_rev = abs(float(column_cadd[i]) - 100 )
                binned_values.append(bin_value_rev)
                
    return(binned_values)
            
def getmean (datalist):
    counter = 0
    allsum  = float(0)
    
    for i in range( len(datalist) ):
        if not datalist[i] == 'NA':
            counter += 1
            allsum += float(datalist[i])
    
    mean = allsum / counter
    
    return(mean)

def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        index = 'NA'
    return(index)

def rank (binned_list):
    # rank returns values from 1 to 100
    rank_list = list()
    for i in range( len(binned_list) ):
        rankvalue = binned_list[i] % 101 + 1 # 101, because otherwise 0 and 100 % 100 calculates to 0
        rank_list.append(rankvalue)
    
    #rank_list = scipy.stats.rankdata(binned_list) # this does produce the wrong order
    
    return(rank_list)


main()