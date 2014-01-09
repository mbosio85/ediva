import scipy as sp
#import scipy.stats
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
#    'EurEVSFrequecy',
#    'AfrEVSFrequecy',
#    'TotalEVSFrequecy',
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
    ### for DEBUG
    pp = pprint.PrettyPrinter(indent = 4)
    ###
    
    parser = argparse.ArgumentParser(description = 'rank SNPs according to their mutation properties')
    
    parser.add_argument('--infile', type=argparse.FileType('r'), dest='infile', required=True, help='comma separated list of SNPs annotated with mutation impact data. [required]')
    parser.add_argument('--outfile', type=argparse.FileType('w'), dest='outfile', required=True, help='comma separated list of SNPs annotated with ranks. [required]')
    
    args = parser.parse_args()
    
    alldata = list(csv.reader(args.infile))
    
    header = alldata.pop(0)
    
    #pp.pprint(alldata[0])
    #print type()
    
    alldata_transpose = zip(*alldata)
    #pp.pprint(alldata_transpose[0])
    
    values2investigate = ['Total1000GenomesFrequency', 'TotalEVSFrequecy', 'SegMentDup', 'Condel', 'VertebratesPhyloP', 'VertebratesPhastCons', 'SIFTScore']
    
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
    
    # print type(binned_values)
    # pp.pprint(binned_values)
    
    rank_product_list = list()
    for i in range( len(binned_values) ):
        rank_product = float( ( ranked_maf[i] * ranked_segdup[i] * ranked_condel[i] * ranked_phastcons[i]) ) / ( 100**2 ) # 4 tools deliver information, decrease the numeric value to more graspable values ### currently deleted * ranked_phylop[i]
        #print 'ranked_maf[i] %f, ranked_segdup[i] %f, ranked_condel[i] %f, ranked_phylop[i] %f' % (ranked_maf[i], ranked_segdup[i], ranked_condel[i], ranked_phylop[i])
        #pp.pprint(rank_product)
        rank_product_list.append(rank_product)
    
    for i in range( len(alldata) ):
        alldata[i].append(rank_product_list[i])
    
    outcsv = csv.writer(args.outfile)
    header.append('rank')
    outcsv.writerow(header)
    for line in alldata:
        outcsv.writerow(line)
    
    exit(0) #DEBUG

def binning (alldata, header, parameter):
    
    sub_pp = pprint.PrettyPrinter(indent = 6) #DEBUG
    binned_values = list()
    
    # MAF
    if parameter == 'MAF':
        index_1000G  = identifycolumns(header, 'Total1000GenomesFrequency')
        index_EVS    = identifycolumns(header, 'TotalEVSFrequecy')
        column_1000G = alldata[index_1000G]
        column_EVS   = alldata[index_EVS]
        
        array_1000G = sp.array(column_1000G, dtype=sp.float32)
        mean_1000G  = sp.mean(array_1000G)
        
        array_EVS = sp.array(column_EVS, dtype=sp.float32)
        mean_EVS  = sp.mean(array_EVS)
        
        
        for i in range( len(column_1000G) ):
            MAF = max( round(float(column_1000G[i]), 2), round(float(column_EVS[i]), 2) )
            bin_value = int(MAF * 100) # (good 1-100 bad)
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
                bin_value = int(segdup * 100)
                binned_values.append(bin_value) # (good 1-100 bad)
    
    # Condel
    elif parameter == 'condel':
        index_condel  = identifycolumns(header, 'Condel')
        column_condel = alldata[index_condel]
        
        mean_condel = getmean(column_condel)

        for i in range( len(column_condel) ):
            if column_condel[i] == 'NA':
                condel = round(mean_condel, 2)
                bin_value = int(condel * 100) # binning of condel values (bad 1-100 good)
                bin_value = abs(bin_value - 101 ) # reverse order
                binned_values.append(bin_value) 
            else:
                condel = round(float(column_condel[i]), 2)
                bin_value = int(condel * 100) # binning of condel values (bad 1-100 good)
                bin_value = abs(bin_value - 101 ) # reverse order
                binned_values.append(bin_value)
    
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
                bin_value = abs(bin_value - 101 ) # reverse order
                binned_values.append(bin_value)
            else:
                phylop = round(float(column_phylop[i]), 2)
                bin_value = int(phylop * 100)  # [bad 1 - 100 good]
                bin_value = abs(bin_value - 101 ) # reverse order
                binned_values.append(bin_value)
    
    # PhastCons
    elif parameter == 'PhastCons':
        index_phastcons  = identifycolumns(header, 'VertebratesPhastCons')
        column_phastcons = alldata[index_phastcons]
        
        mean_phastcons = getmean(column_phastcons)
        #print mean_phastcons # DEBUG
        
        #NAcount = 0 # DEBUG
        for i in range( len(column_phastcons) ):
            if column_phastcons[i] == 'NA':
                #NAcount += 1 # DEBUG
                #print 'found NA in phastcons' #DEBUG
                phastcons = round(mean_phastcons, 2)
                bin_value = int(phastcons * 100)
                binned_values.append(bin_value) # [good 1 - 100 bad]
            else:
                phastcons = round(float(column_phastcons[i]), 2)
                if phastcons == 0: phastcons = 0.001
                bin_value = int(phastcons * 100)  # [good 1 - 100 bad] a value of 1 should get a small rank
                bin_value = abs(bin_value - 101 ) # reverse order
                binned_values.append(bin_value)
        #print 'phastcons NA: %d, all: %d' % (NAcount, len(column_phastcons)) # DEBUG
    
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
        rankvalue = binned_list[i] % 100 + 1 
        rank_list.append(rankvalue)
    
    return(rank_list)


main()