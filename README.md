
Exome Disease Variant analysis (eDiVa), is a variant annotation and prioritization pipleine which
- Integrates multiple reliable information sources to better assess variants deleteriousness
- Integrates multiple damage predictors into one single score
- Helps researchers with an easy-to-use workflow to decide for potentially disease causing genic variants
- It includes information from human genome and transcriptome sequencing, data from different omics platforms specialized in the analysis of rare and complex diseases

This project developed anintegrated platform supporting the collection and storage of genome variation, expression, omics and clinical data to provide diagnostics assistance, to elucidate pathogenicity and to identify causal coding gene mutations.

Please check the Wiki for all information about installing, requirements and running of the pipeline.




----
Content from Wiki 
----

# Requirements 

List of external tools 

* [nextflow](https://www.nextflow.io/). 
* [Samtools](http://samtools.sourceforge.net/) v1.3
* Python 2.7
* [Bedtools](http://bedtools.readthedocs.io/en/latest/) 2.25
* [Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Picard](https://broadinstitute.github.io/picard/) 1.119
* [GATK](https://software.broadinstitute.org/gatk/) 3.3+  
* [Bwa](http://bio-bwa.sourceforge.net/) 0.7.10
* [Tabix](http://www.htslib.org/download/)
* [Bcftools](https://samtools.github.io/bcftools/bcftools.html) 
* [dbNSP](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/) database split by snp and indels with
    * bcftools view -v snps dbsnp.vcf.gz | bbcftools norm  -m - > snps.vcf
    * bcftools view -v indels dbsnp.vcf.gz | bbcftools norm  -m - > indels.vcf

Python packages required:

*  Bio
*  cPickle
*  difflib
*  drmaa
*  getopt
*  gzip
*  hashlib
*  matplotlib
*  pyplot
*  MySQLdb
*  mysql-connector-python-2.0.4
*  ntpath
*  numpyasnp
*  pysam
*  readline
*  scipy
*  struct
*  subprocess
*  tokenize
*  urllib2
*  xlrd
*  xlsxwriter
*  zlib



-----

#  NextFlow Configuration 

Once Installed NextFlow, you can run eDiVA with simple one-liners, provided you configure properly NextFlow.
To do so, you will need to edit the profded nextfow.config file, adding your own parameters.

The most basic setup is to run eDiVA locally, for this you simply need to adapt the paths so they will refer to the installed tools in your machine.

**Here an example of the default nextflow.config file for local execution**
```
process {
  executor='local'
  }

env {
    REF='/users/GD/resource/human/hg19/bwa7/hg19.fasta'
    SHOREREF='/users/GD/resource/human/hg19/shore/hg19.fasta.shore'
    DBINDEL='/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.indels.vcf'
    DBSNP='/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.snps.vcf'
    BWA='/users/GD/tools/bwa/bwa-0.7.10/bwa'
    EDIVA='/users/so/mbosio/ediva/edivatools-code/'
    GATK='/users/GD/tools/GATK/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
    PICARD='/users/GD/tools/picard/picard-tools-1.119/'
    SAMTOOLS='/users/GD/tools/samtools/samtools-1.3.1/samtools'
    NOVOSORT='/users/GD/tools/novocraft/novosort/novosort'
    BEDTOOLS='/users/GD/tools/bedtools/bedtools-latest/bin/'
    CLINDEL='/users/GD/tools/clindel/bin/shore'
    EXOME='/users/GD/resource/human/probesets/merged_kits/allkit.clean.sort.merge.bed'
    BEDTOOLS='/users/GD/tools/bedtools/bedtools-latest/bin/'
    FASTQC='/users/GD/tools/FastQC/FastQC-0.11.5/fastqc'
    PYTHON='/software/so/el7.2/Python-2.7.13/bin/bin/python'
}
```

**HPC environment execution configuration **
If you plan to run eDiVA in an HPC environment, please read NextFlow [documentation](https://www.nextflow.io/docs/latest/process.html) about how to do it 

A simple editing for an SGE environment is to change the process field above with the following:
```
    memory '2 GB'
    queue 'long'
    cpus 8
    executor 'sge'
```

If you plan to run each job with a different configuration, bear in mind this is possible to do by overriding the executor parameter when calling nextflow scripts.


----


# Step 1 eDiVA-Predict

This step describes how to launch eDiVA predict on all new samples coming out of  a sequencing machine.

```
nextflow run ~/ediva/edivatools-code/Predict/eDiVA-Predict.nf \
                    --NAME sample_name \
                    --READ1 fastq.read1.gz \
                    --READ2 fastq.read2.gz \
                    --AFFECTED 1
                    --CPU 2 \
                    --OUTF output_folder/ \
                    -w work_folder/  
```

**Parameters**
* NAME : Sample name you want to be in the Bam and VCF files
* READ1 and READ2 : fastq files with the raw reads
* AFFECTED : set to 1 if the sample is affected by the disease. This is optional and default is 0 
* CPU : number of CPUs assinged to the task, 2 is the default
* OUTF: output folder where to put the Bam and VCF file. One foder per sample is recommended since the final output has alwasy the same name 'all_variants.vcf'
* -w : NextFlow work directory. It is important if you specify your own one for the manual cleanup of NextFlow or for the resuming of interrupted jobs

**Outputs**
For each sample you will obtain:
*  Bam file of aligned and filtered reads
*  all_variants.vcf  with the SNPs and INDELs  calls
*  sample_info.txt [required for the next steps]
*  Quality_control folder with the fastqc results and other quality metrics for the sample.



----

# Step 1.5 eDiVA_prepare_for_annotation

Once all samples from a family are processed and all variants have been called, this step takes care to perform the joint multisample calling of vairant positions, in order to gather all information for all samples.

```
nextflow run eDiVA_prepare_for_annotation.nf \
        --OUTF annotation/ \
        --FAMILY_INFO family_config \
        --CPU 2
```

**Parameters**

* OUTF : output folder
*  CPU : number of CPU to use for the task
*  FAMILY_INFO : it is a text file obtained by concatenating the individual *sample_info.txt* of all samples, one per line.
    *  It is important here to double check the text file to verify the AFFECTED column is correct, 0 for non-affected samples and 1 for the affected ones

**Outputs**

*  VCF file with the multisample call : combined.variants.supplement.vcf
*  Pedigree file for the next steps: pedigree.tree


----

# Step 2-3 eDiVA-Annotate + eDiVA-Score

Once you have produced the combined multisample VCF file, it is time to annotate all variants with eDiVA-Annotate and rank them  with eDiVA-Score

```
nextflow run eDiVA-Annotate.nf \
        --OUTF ./ \
        --VCF testname.vcf \
        --CPU 2 
```

This routine performs the eDiVA-Annotate part, adding annotations from several information sources for all variants in the VCF, and automatically ranks them with eDiVA-Score script. 
The two steps are joined because eDiVA-Score requires eDiVA-Annotate input columns to be calculated on the fly.
There is an available repository of pre-calculated eDiVA-Scores for all exonic variants at [???](http://404.htm)

**Output**

* CSV file with annotated variants and scores 
    * combined.variants.supplement.ranked.csv

----


# Step 4 eDiVA-Prioritize

This is the final step of the pipeline where variants are filtered and selected, depending on the pedigree information, the disease inheritance profile, and the strictness of the filter:

**Command example**
```
nextflow run eDiVA-predict.nf --OUTF ./\
        --ANNOTATED combined.variants.supplement.ranked.csv \ 
        --MODE standard \
        --PEDIGREE pedigree.tree \
        --FAMILY_TYPE trio \
        --HPO hpo_list.txt \
        --EXCLUSIONLIST edivapath/Resource/gene_exclusion_list.txt \
        --INHERITANCE choose
```


**Params**

* ANNOTATED: The ranked.csv file output from eDiVA-Annotate.nf
*  OUTF  : Output folder
*  MODE : standard or strict
*  PEDIGREE : pedigree.tree file output from Step 1.5
*  FAMILY_TYPE  : trio (for parent-child trios) or family (other analyses)
*  HPO_list= Text file with HPO terms to be used for soft filtering. It can point to an empty file such as in Resource/empty_list.txt'* .  
*  EXCLUSIONLIST : Genes excluded because normally are false postitives . A default list is present in *'//Resource/gene_exclusion_list.txt* but it can point to an empty file 
*  INHERITANCE : Inheritance mode for the disease to chose among
    *  'dominant_inherited'
    *  'recessive'
    *   'dominant_denovo'
    *   'Xlinked'
    *  'compound'
 
 
 
 **Output**
 
*  Folder with the inheritance type e.g dominant_denovo
    *  combined.variants.supplement.dominant_inherited.csv : All variants annotated with an extra column explaining if they are compatible with the inheritance pattern provided
    *  combined.variants.supplement.filtered.dominant_inherited.csv : Only the variants passing the filters are reported here, sorted by eDiVA-Score
    *  combined.variants.supplement.filtered.dominant_inherited.csv.xlsx : Excel file with the selected variants, with extra OMIM and Clinvar information where available
*  variant_prioritization_report.xlsx 
    *  This file is updated with the current information. It is an excel file with one spreadsheet for each inheritance type. In this way you can have all analyses results in one file only
