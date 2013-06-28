#!/usr/bin/python

import argparse
import csv
import sys
import tools

def main(args):
   
    #### argument reading
    try:
        snpFH = open(sys.argv[1])
    except IOError:
        sys.exit("ERROR: could not open SNP file")
    except IndexError:
        usage()
        sys.exit("ERROR: not enough arguments given")
    
    try:
        famFH = open(sys.argv[2])
    except IOError:
        sys.exit("ERROR: could not open Fam file")
    except IndexError:
        usage()
        sys.exit("ERROR: not enough arguments given")
    
    try:
        job = sys.argv[3]
    except IOError:
        usage()
        sys.exit("ERROR: not enough arguments given")
    
    if job == 'compound':
        try:
            cutoff = float(sys.argv[4])
        except IndexError:
            print 'setting MAF cutoff to default = 1'
            cutoff = float(1)
            
    
    try:
        outname = '%s.familyfiltered' % sys.argv[1]
        outFH = open(outname, 'w')
    except IOError:
        sys.exit("ERROR: could not open output file")
    
    #### read family file
    family = dict()
    affection = dict()
    
    for line in famFH:
        if line.startswith('family'):
            continue
        
        splitline = line.split()
        
        family[splitline[1]] = splitline[0]
        affection[splitline[1]] = splitline[2]
    
    #### read snp file and judge family info and SNP distribution
    csvfile = csv.reader(snpFH)
    
    index_SNPhighqual = 0
    index_SNPlowqual = 0
    index_Gene = 0
    index_HomHet = 0
    index_Refhom = 0
    index_Refhet = 0
    index_EVSMAF = 0
    index_1000GMAF = 0
    
    genesHit = dict()
    
    if job == 'denovo':
            print 'highlighting de novo variants'
    elif job == 'recfamily':
        print 'highlighting variants fitting to a recessive disease model'
    elif job == 'domfamily':
        print 'highlighting variants fitting to a dominant disease model'
    elif job == 'homozygous':
        print 'looking for variants that are het in parents and hom in offspring'
    elif job == 'compound':
        print 'looking for variants affecting the same gene, applying a cutoff of %f' % cutoff
    else:
        sys.exit('job id not understood, sorry. Select one of: denovo, domfamily, recfamily, homozygous, compound')

    
    for row in csvfile:
        # find positions of some headers
        if row[0] == 'Chr':
            header = row[0:9]
            header.append(job)
            #header.append('unfit')
            header.extend(row[9:])
            header.append('\n')
            outFH.write(','.join(header))
            
            index_SNPhighqual = row.index('SNPcall:highQual')
            index_SNPlowqual = row.index('SNPcall:lowQual')
            index_Gene = row.index('Gene')
            index_HomHet = row.index('SampleDetails(SNP[0|1|2]Ref[0|1|2])')
            index_Refhom = row.index('Refcall:hom')
            index_Refhet = row.index('Refcall:het')
            try:
                index_EVSMAF = row.index('EVS_MAF_TAC')# indels do not harbor Allele Frequencies
                index_1000GMAF = row.index('MAF_1000G')
            except:
                index_EVSMAF = 'NA'
                index_1000GMAF = 'NA'
            
            continue
        
        SNPhighqual = row[index_SNPhighqual]
        SNPlowqual = row[index_SNPlowqual]
        Refhom = row[index_Refhom]
        Refhet = row[index_Refhet]
        
        HomHet = row[index_HomHet]
        Gene = row[index_Gene]
        
        if index_EVSMAF != 'NA' and index_1000GMAF != 'NA':
            EVSMAF = float(row[index_EVSMAF])
            MAF1000G = float(row[index_1000GMAF])
            
            MAF = max(EVSMAF, MAF1000G) # use the maximum MAF of EVS or 1000Genomes
            
        else:
            MAF = 0 # indels do not harbor Allele Frequencies that's why index_MAF is 'NA'
        
        
        
        
        detailed_family = dict()
        
        if job == 'denovo':
            #print 'highlighting de novo variants'
            detailed_family = tools.denovo(SNPhighqual, SNPlowqual, Refhom, Refhet, affection, family)
        elif job == 'recfamily':
            #print 'highlighting variants fitting to a recessive disease model'
            #sys.exit('currently not implemented')
            detailed_family = tools.recfamily(SNPhighqual, SNPlowqual, Refhom, Refhet, HomHet, affection, family)
        elif job == 'domfamily':
            #print 'highlighting variants fitting to a dominant disease model'
            detailed_family = tools.domfamily(SNPhighqual, SNPlowqual, Refhom, Refhet, HomHet, affection, family)
        elif job == 'homozygous':
            #print 'looking for variants that are het in parents and hom in offspring'
            detailed_family = tools.homozygous(HomHet, affection, family)
        elif job == 'compound':
            #print 'looking for variants affecting the same gene'
            # this has to be done in another part of the script...
            genesHit = tools.compoundGenelist(HomHet, Gene, affection, family, genesHit, MAF, cutoff)
        else:
            sys.exit('job id not understood, sorry. Select one of: denovo, domfamily, recfamily, homozygous, compound')
        
        if job != 'compound':
            #write output to file
            fit = list()
            unfit = list()
            
            for element in detailed_family: # element is the family ID given in pedigree information
                if detailed_family[element] == 1:
                    fit.append(element)
                elif detailed_family[element] == 0:
                    unfit.append(element)
                else:
                    sys.stderr.write('error in assigning pattern status for %s %s \n' % (row[0], row[1]))
            
            fit = tools.fillempty(fit)
            unfit = tools.fillempty(unfit)
            
            printrow = row[0:9]
            printrow.append(' '.join(map(str, fit)))
            printrow.extend(row[9:])
            printrow.append('\n')
            outFH.write(','.join(map(str, printrow)))
    
    
    if job == 'compound':
        
        #identify potential compound hets
        potentialCompounds = tools.compound(affection, family, genesHit)
        
        ### compound het:
        # compound het is a gene that is hit at least twice in het state in the affected
        # and is het in one of the unaffected parents at these positions
        # rewind file handle
        snpFH.seek(0)
        
        all_potentials = list()
        
        for fam in potentialCompounds:
            all_potentials.extend(potentialCompounds[fam])
    
        for row in csvfile:
            # find positions of some headers
            if row[0] == 'Chr':
                
                index_SNPhighqual = row.index('SNPcall:highQual')
                index_SNPlowqual = row.index('SNPcall:lowQual')
                index_Gene = row.index('Gene')
                index_HomHet = row.index('SampleDetails(SNP[0|1|2]Ref[0|1|2])')
                index_Refhom = row.index('Refcall:hom')
                index_Refhet = row.index('Refcall:het')
                
                try:
                    index_EVSMAF = row.index('EVS_MAF_TAC')# indels do not harbor Allele Frequencies
                    index_1000GMAF = row.index('MAF_1000G')
                except:
                    index_EVSMAF = 'NA'
                    index_1000GMAF = 'NA'
                
                continue
            
            # read values
            HomHet = row[index_HomHet]
            Gene = row[index_Gene]
            
            if index_EVSMAF != 'NA' and index_1000GMAF != 'NA':
                EVSMAF = float(row[index_EVSMAF])
                MAF1000G = float(row[index_1000GMAF])
                
                MAF = max(EVSMAF, MAF1000G) # use the maximum MAF of EVS or 1000Genomes
                
            else:
                MAF = 0 # indels do not harbor Allele Frequencies that's why index_MAF is 'NA'
            
            
            if Gene in all_potentials:
                printrow = row[0:9]
                if (MAF > cutoff):
                    printrow.append('0')
                else:
                    printrow.append('1')
                printrow.extend(row[9:])
                printrow.append('\n')
                outFH.write(','.join(map(str, printrow)))
            else:
                printrow = row[0:9]
                printrow.append('0')
                printrow.extend(row[9:])
                printrow.append('\n')
                outFH.write(','.join(map(str, printrow)))
        
    snpFH.close()
    famFH.close()
    outFH.close()
    
def usage():
    print "usage: \n python SNP_applyFamilyModel_2.py file pedigree model \n\n"
    print "# file: output of SNP_comparison_1.6.pl \n"
    print "# pedigree: file containing all family members following the pattern"
    print "@ family\tsample\taffected\n@ 1\tVH001\t1\n@ 1\tVH002\t1\n"
    print "# model: denovo|recfamily|domfamily|compound"
    print "### denovo will add a column containing 1 where the script thinks it identified a denovo variant - works in trios only"
    print "### homozygous will add a column containing 1 where the script thinks it identified a variant that is heterozygous in parents, but homozygous in offspring"
    print "### homozygous - accepts a MAF cutoff to exclude SNPs in compound finding with a MAF > cutoff."
    print "### recfamily will add a column containing 1 where the script thinks it identified a homozygous variant, which is not in unaffected individuals"
    print "### domfamily will add a column highlighting variants that are present in affected family members"
    print "### compound will add a column containing 1 where both parents supply a variant in the same gene, which is heterozygous in all investigated individuals - works in trios only\n\n"

main(sys.argv)