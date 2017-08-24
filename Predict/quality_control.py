
def run_fastqc(read1, read2,fastqc,output):
    Cmd = fastqc + ' -o ' + output +" "+ read1+ " "+ read2 
    print Cmd
    subprocess.call(Cmd,shell=True)
    results = [each for each in os.listdir(output) if each.endswith('fastqc_data.txt')]
    total_reads=1

    for r in results:       
        with open('/'.join([output,r]))as rd:
            for line in rd:
                if 'Total Sequences' in line:
                    ff = line.strip().split('\t')
                    nreads = ff[1]
                    total_reads+=int(nreads)
    print 'Total reads = %d'%total_reads
    return total_reads

def run_flagstat(bamfile,samtools,output):
    
    Cmd = ' '.join([samtools,' flagstat ', bamfile ,' > ',output])
    print Cmd
    subprocess.call(Cmd,shell=True)
    
    covered_reads=0
    with open(output) as rd:
        ll = rd.readline()
        reads=ll.split('in')[0].replace(' ','').split('+')
        for r in reads:
            covered_reads += int(r)
        
    
    #print 'Covered reads = %d'%covered_reads
    return 1

def parse_gatk(report,output):
    subprocess.call("cp %s %s"%(report,output),shell=True)
    pass

def calc_coverage(bamfile,bedfile,samtools,bedtools,coverage_check,output_fig,output_path,total_reads):
    # check the bedfile and bamfile are correctly sorted and have the same chromosomes
    # otherwise trhow a warning
    ##BAM chromosome list
    cmd1 ="%s view -H %s"%(samtools,bamfile)
    p1 = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
    output = p1.communicate()[0]
    chrbam = list()
    for line in output.split('\n') :
        if line.startswith('@SQ'):
            tmp = line.split('SN:')[1].split('\t')[0]
            chrbam.append(tmp)
            
    bedchr=list()
    
    with open(bedfile) as rd:
        for line in rd:
            fields=line.split('\t')
            if not(fields[0] in bedchr):
                bedchr.append(fields[0])
    
    ## check the two lists are the same
    if chrbam != bedchr:
        print' WARNING THe chromosome lists from BAM and BED file are different '
        print ' This could lead to a wrong coverage estimation'
        print 'BAM chr: %s'%('\t'.join(chrbam))
        print 'BED chr: %s'%('\t'.join(bedchr))
    
    ## Extract stats about the covered chromosomes in common with Bed file
    
    cmd = '%s idxstats %s'%(samtools,bamfile)
    p1 = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output = p1.communicate()[0]
    bam_stats = dict()
    total_mapped_reads =0
    total_non_mapped_reads=0
    for line in output.split('\n'):
        fields= line.split('\t')
        if fields[0] in bedchr:
            # add to stats.
            bam_stats[fields[0]]= [int(fields[1]), int(fields[2])]
            total_mapped_reads+=int(fields[1])
            total_non_mapped_reads+=int(fields[2])
    ## run coverage check
    #if not(os.path.isfile(coverage_check)):
        
    cmd='%s/coverageBed -b %s -a %s -bed -mean > %s'%(bedtools, bamfile,bedfile,coverage_check)
    print cmd
    subprocess.call(cmd,shell=True) ####<<<<<<<
    
    
    # intersect bam with kit to count the total coverage to be compared later with flagstat original    
    cmd2 = "%s/intersectBed -wa -abam %s -b %s > %s/intersect_attempt.bam"%(bedtools,bamfile,bedfile,output_path)
    print cmd2
    subprocess.call(cmd2,shell=True)
    run_flagstat(output_path+'/intersect_attempt.bam',samtools,output_path+'/flagstat_kit.txt')
    print 'then here do the count between original flagstat ant this one'
    
    covered_reads=0
    with open(output_path+'/flagstat_kit.txt') as rd:
        ll = rd.readline()
        reads=ll.split('in')[0].replace(' ','').split('+')
        for r in reads:
            covered_reads += int(r)
    ratio =covered_reads/float(total_reads)
    print ratio
    
    ## load coverage information
    DP = list()
    tot_reads = 0
    with open(coverage_check) as rd:
        for line in rd:
            fields=line.strip().split('\t')
            DP.append(float(fields[3]))
            tot_reads+=round(float(fields[3])*(int(fields[2])-int(fields[1])))
    
    
    ## extract numpy basic stats and build a histogram celebrating it
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    import numpy as np
    x= np.array(DP)
    #print DP    
    ## stats
    tot = float(len(x))
    
    mu = np.mean(x)
    med = np.median(x)
    if tot ==0: tot =1
    plus1 =100*sum(x>=1)/tot
    plus10=100*sum(x>=10)/tot
    plus20=100*sum(x>=20)/tot
    print tot
    print med
    print mu
     
    title= "Coverage or EXOME bed file: percent reads captured= %.2f Mean = %f Median = %f"%(100*ratio,mu,med)
    title+='\nPercentages: >1:%f >10:%f >20:%f'%(plus1,plus10,plus20)
    
    
    bin_list = np.linspace(0, 150, 151)
    
    
    
    plt.hist(x,bins=bin_list,normed=1,facecolor='green')
    plt.title(title)
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.axvline(x.mean(), color='b', linestyle='dashed', linewidth=2)
    #plt.show()
    plt.savefig(output_fig)
    # fig = plt.gcf()
    # 
    # plot_url = py.plot_mpl(fig, filename='coverage_plot')
    return None

if __name__=='__main__':
    import sys
    import subprocess
    import os
    import argparse
    from Bio import bgzf
    import re
    import ntpath
    import struct
    import hashlib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    
    #import plotly.plotly as py
    parser = argparse.ArgumentParser()
    parser.add_argument("-bamfile"  , dest='bamfile',required=True, type=str, help="bamfile")
    parser.add_argument("-bedfile"  , dest='bedfile',required=True, type=str, help="bedfile")
    parser.add_argument("--output","-o"  , dest='o',required=True, type=str, help="output folder")
    parser.add_argument("--samtools","-s"  , dest='samtools',required=True, type=str, help="samtools_path")
    parser.add_argument("--bedtools","-b"  , dest='bedtools',required=True, type=str, help="bedtools path")
    parser.add_argument("--fastqc","-f"  , dest='fastqc',required=True, type=str, help="fastqc path")
    parser.add_argument("--read1", dest='read1',required=True, type=str, help="fastqc path")
    parser.add_argument("--read2", dest='read2',required=True, type=str, help="fastqc path")
    parser.add_argument("--gatk", dest='gatk',required=True, type=str, help="gatk report file for all variants")


    args = parser.parse_args()
    
    if os.path.isdir(args.o):
        pass
    else:
        subprocess.call('mkdir %s'%args.o,shell = True)
    
    #!/usr/bin/env python
    import subprocess
    bamfile = args.bamfile#'EX080_PN15-1005_VH065.realigned.dm.recalibrated.bam'
    bedfile = args.bedfile#'/users/so/mbosio/WES/sorted_properly.bed'
    samtools = args.samtools#'samtools'
    bedtools = args.bedtools#'bedtools'
    coverage_check=args.o+'/coverage_intermediate_report.bed'#/users/so/mbosio/intersect_with_bedfile2.bed'
    output_fig=args.o+'/coverage.pdf'
    output_flagstat = args.o+'/flagstat_report.txt'
    
    ###?step 1 Fastqc
    total_reads = run_fastqc(args.read1,args.read2,args.fastqc,args.o)
    ##?calculating total number of reads per fastqc file

    #
    #### Step 2 Flagstat    
    covered_reads= run_flagstat(bamfile,samtools,output_flagstat)


    ### Step 3 Coverage profile
    calc_coverage(bamfile,bedfile,samtools,bedtools,coverage_check,output_fig,args.o,total_reads)
    
    
    ### Step 4 GATK evaluation
    
    parse_gatk(args.gatk,args.o+'/gatk_report.txt')
    
    ###?Step 5 Put all of it in a report
    
    
    
    
