import argparse
import csv
import pprint
import re
from scipy.stats import poisson

# Note header:
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
#    'SegMentDup',
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
#    'pos_samples(sampleid>zygosity>Cov>AF)',
#    'neg_samples(sampleid>zygosity>Cov>AF)',
#    'rank']


parser = argparse.ArgumentParser(description = 'filter SNPs according to their family distribution')
    
parser.add_argument('--infile', type=argparse.FileType('r'), dest='infile', required=True, help='comma separated list of SNPs annotated with mutation impact data. [required]')
parser.add_argument('--outfile', type=argparse.FileType('w'), dest='outfile', required=True, help='comma separated list of SNPs annotated with inheritance pattern. [required]')
parser.add_argument('--filteredoutfile', type=argparse.FileType('w'), dest='filteredfile', required=True, help='filtered comma separated list of SNPs annotated with inheritance pattern, only reporting the requested variants. [required]')
parser.add_argument('--family',  type=argparse.FileType('r'), dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')
parser.add_argument('--inheritance', choices=['dominant_denovo', 'dominant_inherited', 'recessive', 'Xlinked', 'compound'], dest='inheritance', required=True, help="""choose a inheritance model [required]
dominant_inherited: used for families
dominant_denovo: apply to novel variants seen in the affected individuals

recessive: detect recessive, homozygous variants (if trio is specified the script will require that all non-affected are heterozygous)
Xlinked: used for X linked recessive variants in trios only
compound: detect compound heterozygous recessive variants
""")
parser.add_argument('--familytype', choices=['trio', 'family'], dest='familytype', required=True, help="choose if the data you provide is a trio or a larger family")
parser.add_argument('--geneexclusion',  type=argparse.FileType('r'), dest='geneexclusion', required=False, help='[Analysis of DNA sequence variants detected by high-throughput sequencing; DOI: 10.1002/humu.22035]. [required]')

args = parser.parse_args()

def main (args):
    pp = pprint.PrettyPrinter( indent=4) # DEBUG
    
    # read the gene exclusion list
    # an empty set will not filter out anything, if gene exclusion list is not provided
    genes2exclude = set()
    if args.geneexclusion:
        for gene in args.geneexclusion:
            gene = gene.rstrip()
            genes2exclude.add(gene)
    
    # read family relationships
    family = dict()
    for line in args.famfile:
        if line.startswith('sample'):
            continue
        line = line.rstrip('\n')
        splitline = line.split('\t')
        family[splitline[0]] = splitline[1]
    
    ####
    # example
    # sample	affected
    # VH017	1
    # VH018	0
    # VH019	0
    ####
    
    # read all data
    alldata = list(csv.reader(args.infile))
    
    header = alldata.pop(0)
    
    
    out         = csv.writer(args.outfile)
    outfiltered = csv.writer(args.filteredfile)

    ############
    # both the output file headers should be consistent
    #############
    
    header.append('inheritance')
    header.append('filter')
    outfiltered.writerow(header)
    
    out.writerow(header)
    
    ############
    # Filter out SNPs from filtered file where : 
    # kick when variant function : synonymous,unknown => line[ExonicFunction(Refseq)]
    # keep when function should be exonic,splicing => line[Function(Refseq)]
    # kick when segmentdup should be 0 => line[SegMentDup]
    # only consider Refseq annotation for performing the first criteria
    #############

    index_sample      = identifycolumns(header, 'samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)') #samples(sampleid>zygosity>DPRef>DPAlt>AF)
    index_MAF1k       = identifycolumns(header, 'Total1000GenomesFrequency')
    index_MAFevs      = identifycolumns(header, 'TotalEVSFrequency')
    index_function    = identifycolumns(header, 'Function(Refseq)')
    index_varfunction = identifycolumns(header, 'ExonicFunction(Refseq)')
    index_segdup      = identifycolumns(header, 'SegMentDup')
    index_gene        = identifycolumns(header, 'Gene(Refseq)')
    
    # init for compound
    initializer = 0
    #if args.inheritance == 'compound':
    #    compound_gene_storage = []
    #    pp.pprint(line)
    #    pp.pprint([line[index_gene], index_gene])
    #    old_gene   = re.sub('\(.*?\)','',line[index_gene]) # PKHD1L1(NM_177531:exon75:c.12330+1G>A) transformed to PKHD1L1. Also works for ";" separated multiple annotations
    #    pp.pprint(old_gene)
    #    old_gene_set = set(old_gene.split(';')) # try to remove double occurences of the same gene names
    #    pp.pprint(old_gene_set)
    #    new_gene   = re.sub('\(.*?\)','',line[index_gene])
    #    new_gene_set = set(new_gene.split(';'))
    
    # start reading data
    for line in alldata:
        
         # init for compound
        if args.inheritance == 'compound' and initializer == 0:
            compound_gene_storage = []
            old_gene   = re.sub('\(.*?\)','',line[index_gene]) # PKHD1L1(NM_177531:exon75:c.12330+1G>A) transformed to PKHD1L1. Also works for ";" separated multiple annotations
            old_gene_set = set(old_gene.split(';')) # try to remove double occurences of the same gene names
            new_gene   = re.sub('\(.*?\)','',line[index_gene])
            new_gene_set = set(new_gene.split(';'))
            initializer = 1
        
        
        try:
            MAF1k      = line[index_MAF1k]
            MAFevs     = line[index_MAFevs]
        except:
            pp.pprint(line)
        try:
            # avoiding problems with NAs
            MAF    = max(float(MAF1k), float(MAFevs))
        except:
            MAF    = 0
        sampledata = line[index_sample]
        
        # filter out genes, that are on the gene exclusion list.
        genenames = set()
        if args.geneexclusion:
            genecolumn   = re.sub('\(.*?\)','',line[index_gene])
            genenames = set(genecolumn.split(';'))
        
        judgement = int()
        ###
        # look for de novo variants
        ###
        if args.inheritance == 'dominant_denovo':
            pp.pprint(line) # DEBUG
            judgement = denovo(sampledata, family)
            pp.pprint(judgement)
            # top SNP
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                out.writerow(line)

            elif judgement == 1 and MAF <= 0.01:
                
                line.append('denovo')
                line.append('pass')
                out.writerow(line)
                
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                        if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                                if (line[index_segdup] == '0'):
                                        outfiltered.writerow(line)
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.01:
                line.append('denovo')
                line.append('filtered')
                out.writerow(line)
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                out.writerow(line)
        
        ###
        # look for familial dominant variants. (being tolerant for missing values)
        ###
        elif args.inheritance == 'dominant_inherited':
            judgement = dominant(sampledata, family)
            # top SNP
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                out.writerow(line)

            elif judgement == 1 and MAF <= 0.05:
                line.append('dominant')
                line.append('pass')
                out.writerow(line)
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                        if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                                if (line[index_segdup] == '0'):
                                        outfiltered.writerow(line)
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.05:
                line.append('dominant')
                line.append('filtered')
                out.writerow(line)
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                out.writerow(line)
            pass
        
        ###
        # look for recessive variants (be aware of trio and family inheritance)
        ###
        elif args.inheritance == 'recessive':
            judgement = recessive(sampledata, family, args.familytype)
            
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                out.writerow(line)

            elif judgement == 1 and MAF <= 0.03:
                line.append('recessive')
                line.append('pass')
                out.writerow(line)
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                        if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                                if (line[index_segdup] == '0'):
                                    outfiltered.writerow(line)
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.03:
                line.append('recessive')
                line.append('filtered')
                out.writerow(line)
            
            # does not fit anything
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                out.writerow(line)
        
        ###
        # look for X linked recessive variants in trios
        ###
        elif args.inheritance == 'Xlinked':
            
            index_chromosome = identifycolumns(header, 'Chr')
            
            # skip all variants not located 
            if line[index_chromosome].lower() == 'x' or line[index_chromosome] == '23':
                pass
            else:
                continue
            
            judgement = xlinked(sampledata, family)
            
            # X linked should only work on trios
            if not args.familytype == 'trio':
                line.append('Trio_only')
                line.append('filtered')
                out.writerow(line)
            
            elif len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                out.writerow(line)
            
            elif judgement == 1 and MAF <= 0.01:
                line.append('Xlinked')
                line.append('pass')
                out.writerow(line)
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                        if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                                if (line[index_segdup] == '0'):
                                    outfiltered.writerow(line)
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.01:
                line.append('Xlinked')
                line.append('filtered')
                out.writerow(line)
            
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                out.writerow(line)
        
        ###
        # look for compound heterozygous variants
        ###
        elif args.inheritance == 'compound':
            """
            valid combinations:
            child	parent	parent
            1/1 invalid
            0/1	0/0	0/0 (denovo)
            0/1	0/1	0/0
            0/1	0/0	0/1
            0/1	0/1	0/1
            
            invalid combinations:
            2&4 (parent 50% chance sick)
            3&4 (parent 50% chance sick)
            """
            
            new_gene = re.sub('\(.*?\)','',line[index_gene])
            new_gene_set = set(new_gene.split(';'))
            
            # if not old_gene == new_gene:
            # check if the names are the same
            # sometimes gene looks like 'PKHD1L1;PKHD1L1'
            if len(old_gene_set - new_gene_set) > 0:
                
                #pp.pprint(['old: ', old_gene, 'new: ',new_gene, 'orig: ', line[index_gene]])
                
                comp_judgement = compoundizer(compound_gene_storage, family, index_sample)
            
                if len(compound_gene_storage) == 1:
                    compound_gene_storage[0].append('NOT_compound')
                    compound_gene_storage[0].append('filtered')
                    out.writerow(compound_gene_storage[0]) # there is only one line
                
                elif comp_judgement == 1:
                    for row in compound_gene_storage:
                        row.append('compound')
                        row.append('pass')
                        out.writerow(row)
                        outfiltered.writerow(row)
                else:
                    for row in compound_gene_storage:
                        row.append('compound')
                        row.append('filtered')
                        out.writerow(row)
                
                # reset values
                compound_gene_storage = []
                old_gene     = new_gene
                old_gene_set = new_gene_set
                pass
            
            judgement = compound(sampledata, family)
            
            # top SNP
            if len(genes2exclude & genenames) > 0:
                line.append('NOT_' + args.inheritance)
                line.append('exclusionlist')
                out.writerow(line)

            elif judgement == 1 and MAF <= 0.03:
                if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                        if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                                if (line[index_segdup] == '0'):
                                    compound_gene_storage.append(line)
                                    pass
            
            # fits inheritance, but is too frequent in the population
            elif judgement == 1 and MAF > 0.03:
                line.append('compound')
                line.append('filtered')
                out.writerow(line)
            
            # does not fit anything
            else:
                
                #pp.pprint(judgement)
                #pp.pprint(line)
                
                line.append('NOT_' + args.inheritance)
                line.append('filtered')
                out.writerow(line)
            
            pass
    
    else:
        # clean up for last gene
        if args.inheritance == 'compound':
            
            comp_judgement = compoundizer(compound_gene_storage, family, index_sample)
            
            if len(compound_gene_storage) == 1:
                compound_gene_storage[0].append('NOT_compound')
                compound_gene_storage[0].append('filtered')
                out.writerow(compound_gene_storage[0])
            
            elif comp_judgement == 1:
                for row in compound_gene_storage:
                    row.append('compound')
                    row.append('pass')
                    out.writerow(row)
                    outfiltered.writerow(row)
            else:
                for row in compound_gene_storage:
                    row.append('compound')
                    row.append('filtered')
                    out.writerow(row)
        

###################################################
# sub routines
###################################################

def compound(sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')

    check_samples = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
        features    = sam.split('>')
        # error catching because of wrong splitting, e.g. 40ACVi>0/1>99>0.333;0.167,40ACVm>0/1>99>0.333;0.167,40ACVp>0/2>99>0.333;0.167
        if len(features) == 1:
            continue
        
        name        = features[0]
        zygosity    = features[1]
        refcoverage = features[2] # could be numeric or .
        altcoverage = features[3] # could be numeric or .
        
        #sub_pp.pprint([refcoverage, altcoverage])
        
        if check_samples.has_key(name):
            check_samples[name] = 1

        # homo alt is not expected in compound
        if zygosity == '1/1' and family[name] == '1':
            #sub_pp.pprint("dropped out in 1/1 1")
            judgement = 0
            break
        
        # het is good (could be an inherited variant or de novo)
        if zygosity == '0/1' and family[name] == '1':
            #sub_pp.pprint("dropped out in 0/1 1")
            judgement = 1
            continue
        
        # het or hom ref for parents is good
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '0':
            #sub_pp.pprint("dropped out in 0/0 0")
            judgement = 1
            continue
        
        # parents must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            #sub_pp.pprint("dropped out in 1/1 0")
            judgement = 0
            break
        
        # offspring should have the variant
        elif zygosity == '0/0' and family[name] == '1':
            #sub_pp.pprint("dropped out in 0/0 1")
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            # use poisson distribution
            # poisson_miss = poisson.cdf(0.0, float(ND_coverage)/2)
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue
            
            try:
                int(refcoverage)
            except:
                sub_pp.pprint(sampledata)
                exit(0)
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            #    #sub_pp.pprint("dropped out in ./. 0 8 0")
            #    judgement = 1
            #    continue
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            # poisson for low coverage and percentage for high coverage
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 1
                continue
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            #    #sub_pp.pprint("dropped out in ./. 0 0 8")
            #    judgement = 0
            #    break
            
            # hom alt
            elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                break
            
            # coverage too low?
            else:
                #sub_pp.pprint("dropped out in ./. 0 else")
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':
            #sub_pp.pprint("dropped out in ./. 1")
            judgement = 0
            break
    
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)

def compoundizer(variantlist, family, index_sample):
    
    
    sub_pp = pprint.PrettyPrinter(indent = 5)
    
    #sub_pp.pprint(variantlist)
    
    # initialize
    ticker_dict = dict()
    for member in family.keys():
        if family[member] == '0':
            ticker_dict[member] = list()
    
    name1 = ticker_dict.keys()[0]
    name2 = ticker_dict.keys()[1]
    
    judgement = 0
    
    # check line by line, if this variant could support a compound het
    for variantline in variantlist:
        # produce a list with all the sample data
        sampledata = variantline[index_sample].split(';')
        
        # produce a dictionary to save the zygosities
        zygosities = dict()
        
        if not len(sampledata) == 3:
            print 'Something went wrong. The samples provided did not contain a trio case?'
            judgement = 0
            break
        
        for sam in sampledata:
            features = sam.split('>')
            name        = features[0]
            zygosity    = features[1]
            refcoverage = features[2] # could be numeric or .
            altcoverage = features[3] # could be numeric or .
            
            # check if, we are looking at the offspring
            if not ticker_dict.keys()[0] == name and not ticker_dict.keys()[1] == name:
                continue
            
            zygosities[name] = zygosity
        
        # check the entered values for possible compound supporters
        # denovo
        if zygosities[name1]   == '0/0' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(0.5)
            ticker_dict[name2].append(0.5)
            pass
    
        elif zygosities[name1]   == '0/1' and zygosities[name2] == '0/0':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == '0/0' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass
        
        elif zygosities[name1]   == '0/1' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == '0/1' and zygosities[name2] == './.':
            ticker_dict[name1].append(1)
            ticker_dict[name2].append(0)
            pass
        
        elif zygosities[name1]   == './.' and zygosities[name2] == '0/1':
            ticker_dict[name1].append(0)
            ticker_dict[name2].append(1)
            pass
        
        pass
    
    #sub_pp.pprint(ticker_dict)
    
    # check if there are enough supporters
    name1_sum = sum(ticker_dict[name1])
    name2_sum = sum(ticker_dict[name2])
    
    if (name1_sum + name2_sum) >= 2 and name1_sum >= 1 and name2_sum >= 1:
        judgement = 1
    else:
        judgement = 0
    
    return(judgement)

def denovo(sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    samples = sampledata.split(';')
    
    judgement = 0
    
    check_samples = dict()
    
    # create data structure for completeness check
    for sam in family.keys():
         check_samples[sam] = 0

    # go into the variant data
    for sam in samples:
        features    = sam.split('>')
        name        = features[0]
        try:
            zygosity    = features[1]
        except:
            sub_pp.pprint(samples)
        refcoverage = features[2] # could be numeric or .
        altcoverage = features[3] # could be numeric or .
        
        # sample info complete?
        if check_samples.has_key(name):
                  check_samples[name] = 1
        
        #sub_pp.pprint(family)
        
        # heterozygous in affected individual - good
        if zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue
        
        # hom ref, not affected - good
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        
        # heterozygous in non-affected - bad
        elif zygosity == '0/1' and family[name] == '0':
            #sub_pp.pprint("heterozygous in non-affected")# DEBUG
            judgement = 0
            break
        
        # hom ref in affected - bad
        elif zygosity == '0/0' and family[name] == '1':
            #sub_pp.pprint("hom ref in affected")# DEBUG
            judgement = 0
            break
        
        # homozygous can't be denovo
        elif zygosity == '1/1':
            #sub_pp.pprint("homozygous can't be denovo")# DEBUG
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue
            
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                #sub_pp.pprint("coverage")# DEBUG
                judgement = 0
                break
            
            # hom ref, non called genotype
            # poisson for low coverage and percentage for high coverage
            # poisson 10 reads (poisson average rate of success = 5) and alt reads = 0 - should get still accepted
            elif (poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007) and (altcoverage / coverage <= 0.05):
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
                judgement = 1
                continue
            
            ## het, non called genotype
            #elif int(altcoverage) >=1:
            #    judgement = 0
            #    break
            
            # not necessary to check for hom alt
            
            # coverage too low?
            else:
                #sub_pp.pprint("else")# DEBUG
                #sub_pp.pprint(poisson.cdf( float(altcoverage), float(coverage)/2 ))# DEBUG
                #sub_pp.pprint(altcoverage / coverage)# DEBUG
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':
            
            # except if vcf file was not supplemented by pileup data
            # accept variants which could not be called in the parents
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue
            
            # do not be that grateful, if there is coverage data
            # if the SNP caller could not call a genotype in an area, where there was coverage
            # we rather trust the SNP caller, than starting to call SNP based on pileup coverage
            else:
                judgement = 0
                break
    
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break
    
    return(judgement)

def dominant(sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')
        
    check_samples = dict()
    
    for samp in family.keys():
         check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
        features    = sam.split('>')
        name        = features[0]
        zygosity    = features[1]
        refcoverage = features[2] # could be numeric or .
        altcoverage = features[3] # could be numeric or .
        
        if check_samples.has_key(name):
            check_samples[name] = 1
    
        # affected family members should have the mutation (hom ref not allowed)
        if zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break
        
        # affected family members might be het
        elif zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue
        
        # affected family members might be hom alt
        # that's the major difference to de novo...
        elif zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue
        
        # non-affected family members must not have the mutation - hom ref is OK
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected family members must not have the mutation - het is bad
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break
        
        # non-affected family members must not have the mutation - hom alt is worst
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            
            # if vcf file was not supplemented
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                continue
            
            # do not do any other judgements, i.e.
            # being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            # also tolerate low coverage
            judgement = 1
            continue
        
        # accept some missing values for affected individuals
        elif zygosity == './.' and family[name] == '1':
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue
            
            # do not do any other judgements, i.e.
            # being tolerant for refcoverage >= 8 and altcoverage == 0, because it could be that the carrier is not yet sick
            # also tolerate low coverage
            judgement = 1
            continue
        
        else:
            pass
    
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return(judgement)

def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        exit("ERROR: %s column could not be identified in the annotated data" % question)
    return( int(index) )

def recessive(sampledata, family, familytype):
    #sub_pp = pprint.PrettyPrinter(indent = 8)
    
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')

    check_samples = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
        features    = sam.split('>')
        name        = features[0]
        zygosity    = features[1]
        refcoverage = features[2] # could be numeric or .
        altcoverage = features[3] # could be numeric or .
        
        if check_samples.has_key(name):
            check_samples[name] = 1
        
        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            #sub_pp.pprint(['1/1 1', name])
            continue
        
        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            #sub_pp.pprint(['0/0 0/1 1 ', name])
            break
        
        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            #sub_pp.pprint(['0/1 0', name])
            continue
        
        # non-affected individuals might be hom ref, if a family is interrogated
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'family':
            judgement = 1
            #sub_pp.pprint(['0/0 0 family', name])
            continue
        
        # non-affected individuals in a trio are the parents and have to be het
        elif zygosity == '0/0' and family[name] == '0' and familytype == 'trio':
            judgement = 0
            #sub_pp.pprint(['0/0 0 trio', name])
            break
        
        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            #sub_pp.pprint(['1/1 0', name])
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.' or refcoverage == '' or altcoverage == '':
                judgement = 1
                #sub_pp.pprint(['./. 0 . .', name])
                continue
            
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                judgement = 1
                #sub_pp.pprint(['./. 0 8 0', name])
                continue
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            elif poisson.cdf( float(refcoverage), float(coverage/2) ) <= 0.007  and refcoverage / coverage <= 0.05:
                judgement = 0
                #sub_pp.pprint(['./. 0 0 8', name])
                break
            
            # het, which is OK
            #elif int(refcoverage) >= 8 and not int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) >= 0.007 or altcoverage / coverage >= 0.05:
                judgement = 1
                #sub_pp.pprint(['./. 0 8 not 0', name])
                continue
            
            # coverage too low?
            else:
                
                # accept missing values in family interrogations
                if familytype == 'family':
                    judgement = 1
                    #sub_pp.pprint(['family', name])
                    continue
                # do not accept missing values in trio setups
                elif familytype == 'trio':
                    judgement = 0
                    #sub_pp.pprint(['trio 0', name])
                    break
                
                # for security reasons
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            
            judgement = 0
            break
    
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)

def xlinked(sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')

    check_samples = dict()
    
    # collect the two parents - one should be het (mother), one whould be hom ref (father)
    # finally should contain two samples, one het & one hom ref
    inheritance_logic = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
        features    = sam.split('>')
        name        = features[0]
        zygosity    = features[1]
        refcoverage = features[2] # could be numeric or .
        altcoverage = features[3] # could be numeric or .
        
        if check_samples.has_key(name):
            check_samples[name] = 1
        
        if family[name] == '0':
            inheritance_logic[name] = zygosity
        
        # affected individuals have to be homozygous
        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue
        
        # affected individuals should not be hom ref or het
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            break
        
        # non-affected individuals might be het
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected individuals might be hom ref
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        
        # non-affected individuals must not be hom alt
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break
        
        # now a few more complex steps, if genotype is missing (only non-affected individuals should have missing values)
        elif zygosity == './.' and family[name] == '0':
            # which chance has the current read distribution to miss out on an alt read
            # e.g. 10ref, 0alt
            # if coverage > 8, the chance is less than 0.5% to miss out on one alt read
            # http://stattrek.com/online-calculator/poisson.aspx
            
            # if vcf file was not supplemented by pileup data
            # accept variants which could not be called
            if refcoverage == '.' or altcoverage == '.':
                judgement = 1
                continue
            
            refcoverage = float(refcoverage)
            altcoverage = float(altcoverage)
            
            coverage = refcoverage + altcoverage
            
            if coverage == 0:
                judgement = 0
                break
            
            # hom ref
            #if int(refcoverage) >= 8 and int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) <= 0.007 and altcoverage / coverage <= 0.05:
                inheritance_logic[name] = '0/0'
                judgement = 1
                continue
            
            # hom alt
            #elif int(altcoverage) >=8 and int(refcoverage) == 0:
            if poisson.cdf( float(refcoverage), float(coverage)/2 ) <= 0.007 and refcoverage / coverage <= 0.05:
                inheritance_logic[name] = '1/1'
                judgement = 0
                break
            
            # het, which is OK
            #elif int(refcoverage) >= 8 and not int(altcoverage) == 0:
            elif poisson.cdf( float(altcoverage), float(coverage)/2 ) >= 0.007 or altcoverage / coverage >= 0.05:
                inheritance_logic[name] = '0/1'
                judgement = 1
                continue
            
            # coverage too low?
            else:
                
                judgement = 0
                break
        
        # do not accept missing values for affected individuals
        # they should be called hom alt by the SNP caller
        elif zygosity == './.' and family[name] == '1':
            
            judgement = 0
            break
    
    # sanity check
    het_checker = 0
    hom_checker = 0
    
    for values in inheritance_logic.values():
        if values == '0/1':
            het_checker = 1
        if values == '0/0':
            hom_checker = 1
    if het_checker == 1 and hom_checker == 1:
        judgement = 1
    else:
        judgement = 0
    
    # another sanity check
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)

main(args)
