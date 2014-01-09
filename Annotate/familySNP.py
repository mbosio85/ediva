import argparse
import csv
import pprint

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
parser.add_argument('--filteredoutfile', type=argparse.FileType('w'), dest='filteredfile', required=True, help='filtered comma separated list of SNPs annotated with inheritance pattern. [required]')
parser.add_argument('--family',  type=argparse.FileType('r'), dest='famfile', required=True, help='tab separated list of samples annotated with affection status. [required]')
parser.add_argument('--inheritance', choices=['denovo', 'dominant', 'recessive'], dest='inheritance', required=True, help='choose a inheritance model [required]')

args = parser.parse_args()

def main (args):
    ##pp = pprint.PrettyPrinter( indent=4) # DEBUG
    
    # read family relationships
    family = dict()
    for line in args.famfile:
        if line.startswith('sample'):
            continue
        line = line.rstrip('\n')
        splitline = line.split('\t')
        family[splitline[0]] = splitline[1]
    
    #pp.pprint(family) # DEBUG
    
    # read all data
    alldata = list(csv.reader(args.infile))
    
    header = alldata.pop(0)
    index_sample = identifycolumns(header, 'samples(sampleid>zygosity>Cov>AF)')
    #index_neg    = identifycolumns(header, 'neg_samples(sampleid>zygosity>Cov>AF)')
    index_MAF1k  = identifycolumns(header, 'Total1000GenomesFrequency')
    index_MAFevs = identifycolumns(header, 'TotalEVSFrequecy')
    
    
    #pp.pprint(header) # DEBUG
    out         = csv.writer(args.outfile)
    outfiltered = csv.writer(args.filteredfile)

    ############
    # both the output file headers should be consistent
    #############
    
    header.append('inheritance')
    header.append('filter')
    outfiltered.writerow(header)
    
    nonfiltheader = header
    #nonfiltheader.append('filter')
    out.writerow(nonfiltheader)
    
    ############
    # Filter out SNPs from filtered file where : 
    # kick when variant function : synonymous,unknown => line[ExonicFunction(Refseq)]
    # keep when function should be exonic,splicing => line[Function(Refseq)]
    # kick when segmentdup should be 0 => line[SegMentDup]
    # only consider Refseq annotation for performing the first criteria
    #############
    
    index_function = identifycolumns(header, 'Function(Refseq)')
    index_varfunction = identifycolumns(header, 'ExonicFunction(Refseq)')
    index_segdup = identifycolumns(header, 'SegMentDup')


    for line in alldata:
        MAF1k = line[index_MAF1k]
        MAFevs = line[index_MAFevs]
        MAF = max(float(MAF1k), float(MAFevs))
        sampledata = line[index_sample]
        #neg_sample = line[index_neg]
        
        #pp.pprint([MAF1k, MAFevs, type(MAF), sampledata]) # DEBUG
        
        
        judgement = int()
        # look for de novo variants
        if args.inheritance == 'denovo':
            judgement = denovo(MAF, sampledata, family)
            # top SNP
            if judgement == 1 and MAF <= 0.01:
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
        # look for recessive variants
        elif args.inheritance == 'recessive':
            judgement = recessive(sampledata, family)
            
            if judgement == 1 and MAF <= 0.03:
                line.append('recessive')
                line.append('pass')
                out.writerow(line)
		if (line[index_function] == 'exonic' or line[index_function] == 'exonic;splicing' or line[index_function] == 'splicing'):
                	if (line[index_varfunction] != 'synonymous SNV' and line[index_varfunction] != 'unknown' and line[index_varfunction] != 'UNKNOWN'):
                        	if (line[index_segdup] == '0'):
		                		outfiltered.writerow(line)
            else:
                line.append('NOT_' + args.inheritance)
                #line.append('bad_inheritance')
                line.append('filtered')
                out.writerow(line)
                
        elif args.inheritance == 'dominant':
            judgement = dominant(MAF, sampledata, family)
            # top SNP
            if judgement == 1 and MAF <= 0.05:
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
        

def denovo(MAF, sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')
    
    judgement = 0
    
    check_samples = dict()

    for samp in family.keys():
         check_samples[samp] = 0

    for sam in samples:
        features = sam.split('>')
        name     = features[0]
        zygosity = features[1]
        coverage = features[2]
        
        #sub_pp.pprint([name, zygosity, coverage])
        #sub_pp.pprint(family[name])
        #sub_pp.pprint(judgement)

        if check_samples.has_key(name):
                  check_samples[name] = 1


        if zygosity == '0/1' and family[name] == '1':
            judgement = 1
            continue
        elif zygosity == '0/0' and family[name] == '0':
            judgement = 1
            continue
        elif zygosity == '0/1' and family[name] == '0':
            judgement = 0
            break
        elif zygosity == '0/0' and family[name] == '1':
            judgement = 0
            break
        # homozygous can't be denovo, ATTENTION: missing values do not get thrown out
        elif zygosity == '1/1':
            judgement = 0
            break
    
    #for vals in check_samples.values():
        #if vals == 0:
            #judgement = 0
            #break

    #print 'judgement is: %d' % judgement     
    return(judgement)

def dominant(MAF, sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')
    
    #print family.keys()
        
    check_samples = dict()
    
    for samp in family.keys():
         check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
             features = sam.split('>')
             name     = features[0]
             zygosity = features[1]
             coverage = features[2]
        
             if check_samples.has_key(name):
                  check_samples[name] = 1
        
             if zygosity == '0/1' and family[name] == '1':
                  judgement = 1
                  continue
             elif zygosity == '0/0' and family[name] == '0':
                  judgement = 1
                  continue
             elif zygosity == '0/1' and family[name] == '0':
                  judgement = 0
                  break
             elif zygosity == '0/0' and family[name] == '1':
                  judgement = 0
                  break
             # homozygous can't be denovo, ATTENTION: missing values do not get thrown out
             elif zygosity == '1/1' or zygosity == './.':
                  judgement = 0
                  break
             else:
                  pass
    #sub_pp.pprint(check_samples)
    
    for vals in check_samples.values():
       if vals == 0:
            judgement = 0
            break

    return(judgement)

def identifycolumns (header, question):
    try:
        index = header.index(question)
    except:
        index = 'NA'
    return(index)

def recessive(sampledata, family):
    sub_pp = pprint.PrettyPrinter(indent = 8)
    # get samples as data structure
    all_samples = dict()
    samples = sampledata.split(';')

    check_samples = dict()
    
    for samp in family.keys():
        check_samples[samp] = 0
    
    judgement = 0
    
    for sam in samples:
        features = sam.split('>')
        # all_samples['sample'] =
        name     = features[0]
        zygosity = features[1]
        coverage = features[2]
        
        if check_samples.has_key(name):
            check_samples[name] = 1


        if zygosity == '1/1' and family[name] == '1':
            judgement = 1
            continue
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '0':
            judgement = 1
            continue
        elif zygosity == '1/1' and family[name] == '0':
            judgement = 0
            break
        elif ( zygosity == '0/0' or zygosity == '0/1' ) and family[name] == '1':
            judgement = 0
            break
    
    for vals in check_samples.values():
        if vals == 0:
            judgement = 0
            break

    return(judgement)

main(args)
