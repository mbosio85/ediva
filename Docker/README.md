## How to use Dockerfiles for eDiVA ##

eDiVA_code : 
It contains the code needed to run eDiVa.
Alone it can run easily for eDiVA-Predict, eDiVA-prepare_for_annotation, and eDiVA-prioritize
It does require eDiVA-DB to run eDiVA-annotate though.

Build Instruction:
* Go to eDiVA_code foder:
* Run :> docker build -t ediva:code ./

eDiVA-DB : 

Build Instruction:
* Download the eDiVA_public_omics.sql.gz and eDiVA_annotation.sql.gz
* Place them in the Docker/eDiVA_DB folder
* go to eDiVA_DB folder: 
* Run :> docker build -t ediva:db ./
* Wait a long time because the database is big


######
##
## Running containers
##
#####

Run instructions with nextflow [suggested mode]
* Edit nextflow.config to setup environment and docker capabilities
  * You can use the provided nexflow.config in the Docker folder to activate/modify the parts that are needed
  * It requires to have the container enabled and environment set: PATHS are those FROM THE CONTAINER
  process {
  executor='local'
  }

  process.container = 'ediva'
  process.scratch = true

  docker {
    enabled = true
    temp = '/tmp/'
    runOptions = '--network host  -p 3306 -v ${YOUR_REF_FILE_FOLDER}:/resources/ref/  -v /users/GD/resource/human/hg19/databases/dbSNP/:/resources/dbsnp/  -v /users/GD/resource/human/probesets/merged_kits/:/resources/exome'
  }

  env {
    REF='/resources/ref/hg19.fasta'
    DBINDEL='/resources/dbsnp/dbsnp_138.hg19.indels.vcf'
    DBSNP='/resources/dbsnp/dbsnp_138.hg19.snps.vcf'
    EXOME='/resources/exome/allkit.clean.sort.merge.bed'
    EDIVA='/opt/edivatools-code/'
    PYTHON='/usr/bin/python'
    SAMTOOLS='/usr/local/bin/samtools'
    PICARD='/opt/picard-tools-1.119/'
    GATK='/usr/bin/GenomeAnalysisTK.jar'
    BEDTOOLS='/usr/bin/'
    BWA='/usr/bin/bwa'
    FASTQC='/usr/bin/fastqc/'
  }


eDiVA-Predict : 

nextflow run eDiVA-Predict.nf \
--NAME sample_name \
--READ1 /users/so/nrostan/PHD_exome_course_2015/TheExomeCourse/Data/Case_1_FHHt/CD2224/CD2224.read1.fastq.gz \
--READ2 /users/so/nrostan/PHD_exome_course_2015/TheExomeCourse/Data/Case_1_FHHt/CD2224/CD2224.read2.fastq.gz \
--AFFECTED 1 \
--OUTF outfolder \
--CPU 2  -w work_folder \
-with-docker ediva:code


eDiVA-Prioritize

* First: edit the nexflow.config to mount the samples directories in the docker container file system
  * Add to runOptions: 
  * -v /path_to_sample1/:/samples/sample1/
  * -v /path_to_sample2/:/samples/sample2/
  * -v /path_to_sample3/:/samples/sample3/
* Second: prepare family info file chaging sample files paths with /samples/sampleX/
  * Local file for family_config:
    * ID status vcf bam
    * ID1 1 /path_to_sample1/all_variants.vcf /path_to_sample1/sample1.recalibrated.bam
    * ID2 0 /path_to_sample2/all_variants.vcf /path_to_sample1/sample2.recalibrated.bam
    * ID3 0 /path_to_sample3/all_variants.vcf /path_to_sample1/sample3.recalibrated.bam
  * Changed for docker:
    * ID status vcf bam
    * ID1 1 /samples/sample1/all_variants.vcf /samples/sample1/sample1.recalibrated.bam
    * ID2 0 /samples/sample2/all_variants.vcf /samples/sample2/sample2.recalibrated.bam
    * ID3 0 /samples/sample3/all_variants.vcf /samples/sample3/sample3.recalibrated.bam


nextflow run eDiVA-prioritize.nf --OUTF ./\
        --ANNOTATED testname.supplement.ranked.csv \
        --MODE standard \
        --PEDIGREE pedigree.tree \
        --FAMILY_TYPE trio \
        --HPO hpo_list.txt \
        --EXCLUSIONLIST edivapath/Resource/gene_exclusion_list.txt \
        --INHERITANCE dominant_inherited -with-docker ediva:code

######
##
##  eDiVA-Annotate  
##
######

The annotation part is more complex and require a more elaborations

* To run the annotation we need to have the eDiVA_DB container running :
* Run:
  * docker run --detach --name=ediva_db_instance --env="MYSQL_ROOT_PASSWORD=rootpwd" ediva:db
  * if you want to assign more resources to this container follow Docker instructions
* Then we need to edit nextflow.config options for docker, linking the ediva:code container to the ediva_db_instance so it can query the database
  * Add this line to runOptions '' 
  * --link ediva_db_instance:mysqlsrv-ediva.linux.crg.es
* Launch the container via nextflow:
 

nextflow run eDiVA-Annotate.nf \
        --OUTF ./ \
        --VCF testname.vcf \
        --CPU 2  -w work_folder \
        -with-docker ediva:code



