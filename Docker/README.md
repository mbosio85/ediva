# How to use Dockerfiles for eDiVA 

## eDiVA_code : 
It contains the code needed to run eDiVa.
Alone it can run easily for eDiVA-Predict, eDiVA-prepare_for_annotation, and eDiVA-prioritize
It does require eDiVA-DB to run eDiVA-annotate though.

Build Instruction:
* Go to eDiVA_code folder:
* Download GATK 3.3+ [here](https://software.broadinstitute.org/gatk/download/archive) 
* Copy there the GATK jar file 'GenomeAnalysisTK.jar' **into the eDiVA_code folder**
* Execute:
```
docker build -t ediva:code ./
```
## db:
Build instruction
* go to db folder
* Execute:
```
docker volume create --name=db
```
* Download the eDiVA_public_omics.sql.gz and eDiVA_annotation.sql.gz
```
wget https://public_docs.crg.es/sossowski/MicrobeGenomes/human/eDiVA/eDiVA_DB/eDiVa_public_omics.sql.gz
wget https://public_docs.crg.es/sossowski/MicrobeGenomes/human/eDiVA/eDiVA_DB/eDiVa_annotation.sql.gz
```
* Place them in the Docker/db folder 

## eDiVA-DB : 

Build Instruction:
* go to eDiVA_DB folder: 
* Execute:
```
docker build -t ediva:db ./ 
```
* Launch the eDiVA_DB image as a running mysql server. We will then populate the database with the information we need.
* We also mount the created data volume 'db' so it preserves the data even if we shut down the mysql server later on
```
docker run --detach --name=ediva_database --env="MYSQL_ROOT_PASSWORD=mypassword" -v db:/var/lib/mysql2/  ediva:db
```
* Next step we launch an interface to the mysql server to populate the database
* First we load the eDiVA_public database:
  * **Take care to edit the local path to your Docker/db folder**
```
docker run --rm -ti --name populate-db \ 
--link ediva_database:mysql.srv \ 
-v path_to_Docker/db/:/bin/sql ediva:code /bin/bash \
-c " zcat /bin/sql/eDiVa_public_omics.sql.gz| mysql -u edivapublic -px86d2k1B -h mysql.srv  -D eDiVa_public_omics" 

```
* Second we load the much bigger eDiVA_annotation database:
  * **Take care to edit the local path to your Docker/db folder**
```
 docker run --rm -ti --name populate-db \
 --link ediva_database:mysql.srv \
 -v path_to_Docker/db/:/bin/sql ediva:code /bin/bash \
 -c " zcat /bin/sql/eDiVa_annotation.sql.gz| mysql -u edivapublic -px86d2k1B -h 10.2.0.1 -D eDiVa_annotation"
 
 
```




# Running containers


Run instructions with nextflow: suggested mode
* Edit nextflow.config to setup environment and docker capabilities
  * You can use the provided nexflow.config in the Docker folder to activate/modify the parts that are needed
  * It requires to have the container enabled and environment set: PATHS are those FROM THE CONTAINER
  * Basically you need to edit the fields starting with `-v ` 
  * And you have to edit the REF and EXOME values in the `env `part 
  * **Example of how prepare the nextflow.config file **
~~~
  process {  executor='local' }

  process.container = 'ediva'
  process.scratch = true

  docker {
    enabled = true
    temp = '/tmp/'
    runOptions = '--network host  -p 3306 \
    -v ${YOUR_REF_FILE_FOLDER}:/resources/ref/ \
    -v ${YOUR_DBSNP_FOLDER}:/resources/dbsnp/ \
    -v ${FOLDER_OF_BED_FILE_FOR_YOUR_ANALYSIS}:/resources/exome'
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
~~~

## eDiVA-Predict : 

~~~
nextflow run eDiVA-Predict.nf \
--NAME sample_name \
--READ1 /path_to_sample1/Sample1.read1.fastq.gz \
--READ2 //path_to_sample1/Sample1.read2.fastq.gz \
--AFFECTED 1 \
--OUTF outfolder \
--CPU 2  -w work_folder \
-with-docker ediva:code
~~~

## eDiVA-Prepare for annotation
When handling a trio case for example we need to prepare a multisample VCF file for annotation by combining the individual VCF files from each sample, produced with eDiVA-Predict.
* First: edit the nexflow.config to mount the samples directories in the docker container file system
  * Add to runOptions in **nextflow.config**: 
~~~
   -v /path_to_sample1/:/samples/sample1/
   -v /path_to_sample2/:/samples/sample2/
   -v /path_to_sample3/:/samples/sample3/
~~~
* Then run the container
~~~
nextflow run eDiVA_prepare_for_annotation.nf \
        --OUTF annotation/ \
        --FAMILY_INFO family_config \
        --CPU 2
        -with-docker ediva:code
~~~

*  FAMILY_INFO : it is a text file obtained by concatenating the individual *sample_info.txt* of all samples, one per line.
    *  It is important here to double check the text file to verify the AFFECTED column is correct, 0 for non-affected
    
* Example of a local file for family_config: **fields are separated by tab**
~~~
    ID status vcf bam
    ID1 1 /path_to_sample1/all_variants.vcf /path_to_sample1/sample1.recalibrated.bam
    ID2 0 /path_to_sample2/all_variants.vcf /path_to_sample1/sample2.recalibrated.bam
    ID3 0 /path_to_sample3/all_variants.vcf /path_to_sample1/sample3.recalibrated.bam
~~~
  * Changes needed for docker execution :**fields are separated by tab**
~~~
    ID status vcf bam
    ID1 1 /samples/sample1/all_variants.vcf /samples/sample1/sample1.recalibrated.bam
    ID2 0 /samples/sample2/all_variants.vcf /samples/sample2/sample2.recalibrated.bam
    ID3 0 /samples/sample3/all_variants.vcf /samples/sample3/sample3.recalibrated.bam
~~~    
    
##  eDiVA-Annotate  


The annotation part is more complex and require a more actions

* To run the annotation we need to have the eDiVA_DB container running :
* Run:
~~~
docker run --detach --name=ediva_db_instance --env="MYSQL_ROOT_PASSWORD=rootpwd" ediva:db
~~~
  * if you want to assign more resources to this container follow Docker instructions from the official documentation
* Then we need to edit nextflow.config options for docker, linking the ediva:code container to the ediva_db_instance so it can query the database
  * Add this line to runOptions '' 
~~~
  --link ediva_db_instance:mysqlsrv-ediva.linux.crg.es
~~~
* Launch the container via nextflow:
 
~~~
nextflow run eDiVA-Annotate.nf \
        --OUTF ./ \
        --VCF testname.vcf \
        --CPU 2  -w work_folder \
        -with-docker ediva:code
~~~


## eDiVA-Prioritize

* First: edit the nexflow.config to mount the samples directories in the docker container file system
  * Add to runOptions: 
~~~
   -v /path_to_sample1/:/samples/sample1/
   -v /path_to_sample2/:/samples/sample2/
   -v /path_to_sample3/:/samples/sample3/
~~~
* Second: prepare family info file chaging sample files paths with /samples/sampleX/ as explained above. 

Run it:
~~~
nextflow run eDiVA-prioritize.nf --OUTF ./\
        --ANNOTATED testname.supplement.ranked.csv \
        --MODE standard \
        --PEDIGREE pedigree.tree \
        --FAMILY_TYPE trio \
        --HPO hpo_list.txt \
        --EXCLUSIONLIST edivapath/Resource/gene_exclusion_list.txt \
        --INHERITANCE dominant_inherited -with-docker ediva:code
~~~
 **pedigree.tree has been generated by eDiVA-Annotate execution**

