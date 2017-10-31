# specify how to test it
echo "This test is for eDiVA run with nextflow and Docker"
echo "There is no automated test for other ways of functioning"
echo "Feel free to adapt the code to run the test manually if you want"
echo "This test requires you to have build the docker containers ediva:code and eidva:db "
echo "It also requires you to have an instance of ediva_database running as mysql server"
echo "If you did not do it, launch it with this: "
echo "docker run --detach --name=ediva_database --env="MYSQL_ROOT_PASSWORD=mypassword" -v db:/var/lib/mysql2/  ediva:db"
echo "##"
echo "##"

#download test material
echo ""
echo "Downloading the test data required"
wget https://public_docs.crg.es/sossowski/MicrobeGenomes/human/eDiVA/test_material.zip
unzip test_material.zip
echo "done"


# check for ediva
echo ""
echo "Checking for EDIVAPATH in your environment"
if ! [[ -n "$EDIVAPATH" ]] ; then 
echo "set the EDIVAPATH variable with export EDIVAPATH='' where you put the path to eDiVA"; 
exit  ;
 fi
echo "done"


## check for Nextflow 
echo ""
echo "Checking for Nextflow executable"
 path_to_executable=$(which nextflow 2> /dev/null)
 if [ -x "$path_to_executable" ] ; then
    echo "nextflow found here : $path_to_executable"
 else  echo "Please put nextflow accessible via the PATH variable"; exit
 fi

## check for Docker 
echo ""
echo "Checking for docker executable"
 path_to_executable=$(which docker 2> /dev/null)
 if [ -x "$path_to_executable" ] ; then
    echo "docker is available, good :) "
 else  echo "Please install docker or make it available in PATH environment variable"; exit
 fi


## check for ediva images
echo ""
echo "Checking you build the needed images and have ediva_database up and running"
edivacode=$(docker images|grep ediva|grep code)
edivadb=$(docker images|grep ediva|grep db)
ediva_database=$(docker ps |grep ediva_database|grep "ediva:db")

echo $edivacode

if [ ! -z "$edivacode" ] ; then
    echo "ediva:code is available in the image repository "
 else  echo "Please follow build instruction on eDiVA wiki to build and run the required images (ediva:code, ediva:db ediva_database)"; exit
 fi
if [ ! -z "$edivadb" ] ; then
    echo "ediva:db is available in the image repository "
 else  echo "Please follow build instruction on eDiVA wiki to build and run the required images (ediva:code, ediva:db ediva_database)"; exit
 fi
if [ ! -z "$ediva_database" ] ; then
    echo "ediva_database as instance of ediva:db is up and running. "
 else  echo "Please follow build instruction on eDiVA wiki to build and run the required images (ediva:code, ediva:db ediva_database)"; exit
 fi



## edit the nextflow config with current paths
echo ""
echo "Setting up path variables in nextflow.config"
echo "please check it  so you can edit it manually with your own dbSNP, bedfile and reference genome information"
echo "It is important to edit the variables DBSNP, REF, BDINDEL, EXOME, and the corresponding -v paths "
sed -e "s@path_to_test_folder@$(pwd)@g" nextflow.config.model > nextflow.config


echo "This test requires you to have build the docker containers ediva:code and eidva:db "
echo "It also requires you to have an instance of ediva_database running as mysql server"
echo "If you did not do it, launch it with this: "
echo "docker run --detach --name=ediva_database --env="MYSQL_ROOT_PASSWORD=mypassword" -v db:/var/lib/mysql2/  ediva:db"


##############################################################
##
##            Predict
##
##############################################################
echo ""
echo "eDiVA-Predict on the trio"

### docker run --detach --name=ediva_database --env="MYSQL_ROOT_PASSWORD=rootpwd" ediva:db
sed -i 's/mysqldb/ediva_database/g' nextflow.config 

nextflow run $EDIVAPATH/Predict/eDiVA-Predict.nf --NAME NA12878 --READ1 NA12878/NA12878.read1.fastq.gz --READ2 NA12878/NA12878.read2.fastq.gz --AFFECTED 1 --OUTF NA12878 --CPU 1  -w work_folder -with-docker ediva:code  
nextflow run $EDIVAPATH/Predict/eDiVA-Predict.nf --NAME NA12891 --READ1 NA12891/NA12891.read1.fastq.gz --READ2 NA12891/NA12891.read2.fastq.gz --AFFECTED 0 --OUTF NA12891 --CPU 1  -w work_folder -with-docker ediva:code 
nextflow run $EDIVAPATH/Predict/eDiVA-Predict.nf --NAME NA12892 --READ1 NA12892/NA12892.read1.fastq.gz --READ2 NA12892/NA12892.read2.fastq.gz --AFFECTED 0 --OUTF NA12892 --CPU 1  -w work_folder -with-docker ediva:code

##############################################################
##
## Prepare for annotation
##
##############################################################
echo ""
echo "eDiVA- prepare for annotation on the trio"
echo "ID status vcf bam"|sed -e 's/ /\t/g' > family_config
for i in */sample_info_docker.txt ; do  cat "${i}"  >> family_config ; done

nextflow run $EDIVAPATH/Prioritize/eDiVA_prepare_for_annotation.nf   --OUTF annotation/  --FAMILY_INFO family_config   --CPU 1   -w work_folder -with-docker ediva:code -resume

############################################################## 
##
## annotate
##
##############################################################
echo ""
echo "eDiVA-Annotate"
nextflow run $EDIVAPATH/Prioritize/eDiVA-Annotate.nf  --OUTF annotation/  --VCF annotation/combined.variants.supplement.vcf  --CPU 1  -w work_folder -with-docker ediva:code



##############################################################
##
## prioritise
##
##############################################################
echo ""
echo "eDiVA-Prioritize"
#xlinked
nextflow run $EDIVAPATH/Prioritize/eDiVA-prioritize.nf --OUTF prioritization/ --ANNOTATED annotation/combined.variants.supplement.ranked.csv --MODE standard --PEDIGREE annotation/pedigree.tree --FAMILY_TYPE trio --EXCLUSIONLIST $EDIVAPATH/Resource/gene_exclusion_list.txt --INHERITANCE Xlinked -with-docker ediva:code

#recessive
nextflow run $EDIVAPATH/Prioritize/eDiVA-prioritize.nf --OUTF prioritization/ --ANNOTATED annotation/combined.variants.supplement.ranked.csv --MODE standard --PEDIGREE annotation/pedigree.tree --FAMILY_TYPE trio --EXCLUSIONLIST $EDIVAPATH/Resource/gene_exclusion_list.txt --INHERITANCE recessive -with-docker ediva:code

#compound
nextflow run $EDIVAPATH/Prioritize/eDiVA-prioritize.nf --OUTF prioritization/ --ANNOTATED annotation/combined.variants.supplement.ranked.csv --MODE standard --PEDIGREE annotation/pedigree.tree --FAMILY_TYPE trio --EXCLUSIONLIST $EDIVAPATH/Resource/gene_exclusion_list.txt --INHERITANCE compound -with-docker ediva:code


