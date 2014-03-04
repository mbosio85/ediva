#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
# eDiVaPredict toolbox 
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------


sub usage { print "\n$0 \n usage:\n",
	   "--infolder \t folder containing the raw sequence data in fastq format \n",
	   "--outfolder \t folder to write the analysis to \n",
	   "--qsubname \t name of the script you want to start lateron \n",
	   "--config \t name of the config file to be used \n",
	   "--max_coverage \t used in SNP filtering with samtools [default = 400] \n",
	   "--namestart \t  start of a substring in the read file name (first letter is numbered 1) \n",
	   "--namelength \t length of the part of the filenames which should be taken as sample (and folder) name \n",
	   "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n",
	   "--secondreadextension \t describes the name of the first read file (e.g. 2.fastq.gz or 2_sequence.txt.gz)\n",
	   "--indelcaller \t sets the indel caller to be used. Choose between gatk or clindel [default = clindel]\n",
	   "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 4]\n",
	   "--mem \t amount of Memory dedicated to your job in Gb. The amount of memory must not be bigger than available in the machine. [default = 12]\n",
	   "--help \t\t show help \n";
}


### initialize and read input variables
my $infolder;
my $outfolder;
my $qsub_name;
my $config;
my $max_cov = 400; 
my $nameStart = 'NA';
my $nameLength;
my $firstreadextension;
my $secondreadextension;
my $indelcaller = 'clindel';
my $cpu = 4;
my $mem = 12;
my $help = 0;

GetOptions("infolder=s" => \$infolder, "outfolder=s" => \$outfolder, "qsubname=s" => \$qsub_name, "config=s" => \$config, "max_coverage=i" => \$max_cov, "nameStart=s" => \$nameStart, "nameLength=s" => \$nameLength, "firstreadextension=s" => \$firstreadextension, "secondreadextension=s" => \$secondreadextension, "indelcaller=s" => \$indelcaller, "cpu=i" => \$cpu, "mem=i" => \$mem, "help=s" => \$help);

unless($infolder && $outfolder && $qsub_name && $config && $nameStart ne 'NA' && $nameLength && $firstreadextension && $secondreadextension && $help == 0) {
	usage;
	exit;
}

### read the config file
# "REFERENCE", "SHORE_REFERENCE", "DBINDEL", "DBSNP", "BWA", "GATK", "SAMTOOLS", "NOVOSORT", "PICARD", "BEDTOOLS", "EXOME", "EXOME_SHORE"]
open(FHconf, "<$config");
my ($ediva, $ref, $shore_ref, $dbindel, $dbsnp, $bwa, $gatk, $picard, $samtools, $novosort, $bedtools, $clindel, $exome, $exome_shore);
while (<FHconf>) {
	my @splitline = split('=', $_);
	
	if ($splitline[0] eq 'EDIVA') {
		$ediva = $splitline[1];
	}
	elsif ($splitline[0] eq 'REFERENCE') {
		$ref = $splitline[1];
	}
	elsif ($splitline[0] eq 'SHORE_REFERENCE') {
		$shore_ref = $splitline[1];
	}
	elsif ($splitline[0] eq 'DBINDEL') {
		$dbindel = $splitline[1];
	}
	elsif ($splitline[0] eq 'DBSNP') {
		$dbsnp = $splitline[1];
	}
	elsif ($splitline[0] eq 'BWA') {
		$bwa = $splitline[1];
	}
	elsif ($splitline[0] eq 'GATK') {
		$gatk = $splitline[1];
	}
	elsif ($splitline[0] eq 'SAMTOOLS') {
		$samtools = $splitline[1];
	}
	elsif ($splitline[0] eq 'NOVOSORT') {
		$novosort = $splitline[1];
	}
	elsif ($splitline[0] eq 'PICARD') {
		$picard = $splitline[1];
	}
	elsif ($splitline[0] eq 'BEDTOOLS') {
		$bedtools = $splitline[1];
	}
	elsif ($splitline[0] eq 'CLINDEL') {
		$clindel = $splitline[1];
	}
	elsif ($splitline[0] eq 'EXOME') {
		$exome = $splitline[1];
	}
	#elsif ($splitline[0] eq 'EXOME_SHORE') {
	#	$exome_shore = $splitline[1];
	#}
	
}
close(FHconf);


### prepare alignment and SNP calling pipeline

my @snppipe = ("


REF=$ref
SHOREREF=$shore_ref

DBINDEL=$dbindel
DBSNP=$dbsnp

BWA=$bwa
EDIVA=$ediva
GATK=$gatk
PICARD=$picard
SAMTOOLS=$samtools
NOVOSORT=$novosort
# BCFTOOLS=\$(which bcftools) # CHECK
# VCFUTILS=\$(which vcfutils.pl) # CHECK
BEDTOOLS=$bedtools

CLINDEL=$clindel
# NGSBOX=/users/GD/tools/ngsbox  # CHECK
# RSCRIPT=\$(which Rscript)  # CHECK

EXOME=$exome

### Align reads with bwa
\$BWA mem -M -t $cpu -R \"\@RG\\tID:\$NAME\\tSM:\$NAME\" \$REF \$READ1 \$READ2 | time \$SAMTOOLS view -h -b -S -F 0x900 -  > \$TMPDIR/\$NAME.noChimeric.bam


### check for Quality encoding and transform to 33 if 64 encoding is encountered
OFFSET=\$(\$SAMTOOLS view \$TMPDIR/\$NAME.noChimeric.bam | python \$EDIVA/Predict/whichQuality_bam.py)
if [[ \$OFFSET == 64 ]];
then
	echo \"fixing 64 quality encoding\"
	\$SAMTOOLS view -h \$TMPDIR/\$NAME.noChimeric.bam | python \$EDIVA/Predict/bam_rescale_quals.py - | \$SAMTOOLS view -bS - > \$TMPDIR/\$NAME.transformed.bam
	rm \$TMPDIR/\$NAME.noChimeric.bam
	mv \$TMPDIR/\$NAME.transformed.bam \$TMPDIR/\$NAME.noChimeric.bam
fi


### Sort BAM file
if [ -s \$TMPDIR/\$NAME.noChimeric.bam ];
then
   echo Sort BAM
   \$NOVOSORT --threads $cpu --tmpdir \$TMPDIR --forcesort --output \$TMPDIR/\$NAME.sort.bam -i -m ${mem}G \$TMPDIR/\$NAME.noChimeric.bam
   cp \$TMPDIR/\$NAME.sort.bam* \$OUTF/
   
   # clean up
   rm \$TMPDIR/\$NAME.noChimeric.bam
else
   echo \$TMPDIR/\$NAME.noChimeric.bam not found
   exit
fi



### Local Re-alignment
if [ -s \$TMPDIR/\$NAME.sort.bam.bai ];
then
   echo Local Re-alignment
   echo -e \"\\n #### doing Local Realignment: java -jar \$GATK -nt $cpu -T RealignerTargetCreator -R \$REF -I \$TMPDIR/\$NAME.sort.bam -o \$OUTF/\$NAME.intervals -known \$DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE \\n\"
   java -jar \$GATK -nt $cpu -T RealignerTargetCreator -R \$REF -I \$TMPDIR/\$NAME.sort.bam -o \$OUTF/\$NAME.intervals -known \$DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE
   java -jar \$GATK -T IndelRealigner -R \$REF -I \$TMPDIR/\$NAME.sort.bam -targetIntervals \$OUTF/\$NAME.intervals -o \$TMPDIR/\$NAME.realigned.bam -known \$DBINDEL --maxReadsForRealignment 10000 --consensusDeterminationModel USE_SW --downsampling_type NONE \$FMQ
else
   echo \$NAME.sort.bam.bai not found
   exit
fi



### Duplicate marking
if [ -s \$TMPDIR/\$NAME.realigned.bam ];
then
   # clean up
   rm \$TMPDIR/\$NAME.sort.bam*
   
   echo Duplicate marking
   echo -e \"\\n #### duplicate marking using: java -jar \$PICARD/MarkDuplicates.jar INPUT=\$TMPDIR/\$NAME.realigned.bam OUTPUT=\$TMPDIR/\$NAME.realigned.dm.bam METRICS_FILE=\$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true \\n\"
   java -jar \$PICARD/MarkDuplicates.jar INPUT=\$TMPDIR/\$NAME.realigned.bam OUTPUT=\$TMPDIR/\$NAME.realigned.dm.bam METRICS_FILE=\$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true 
else
   echo \$TMPDIR/\$NAME.realigned.bam not found
   exit
fi



### Base quality recalibration
if [ -s \$TMPDIR/\$NAME.realigned.dm.bam ];
then
   echo -e \" \\n #### Base quality recalibration \\n \"
   java -jar \$GATK -T BaseRecalibrator -nct $cpu --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R \$REF -I \$TMPDIR/\$NAME.realigned.dm.bam -knownSites \$DBSNP --downsampling_type NONE -o \$TMPDIR/\$NAME.recal_data.grp
   java -jar \$GATK -T PrintReads -R \$REF -I \$TMPDIR/\$NAME.realigned.dm.bam -BQSR \$TMPDIR/\$NAME.recal_data.grp -o \$OUTF/\$NAME.realigned.dm.recalibrated.bam
   
else
   echo \$NAME.realigned.dm.bam not found
   exit
fi



### Cleanup
set +e
	rm \$OUTF/\$NAME.realigned.bam
	rm \$OUTF/\$NAME.realigned.bai
	rm \$OUTF/\$NAME.realigned.dm.bam
	rm \$OUTF/\$NAME.realigned.dm.bam.bai
	
	rm \$TMPDIR/\$NAME.bam
	rm \$TMPDIR/\$NAME.realigned.bam
	rm \$TMPDIR/\$NAME.realigned.bai
	rm \$TMPDIR/\$NAME.realigned.dm.bam
	rm \$TMPDIR/\$NAME.realigned.dm.bam.bai
set -e





### GATK: Call SNPs and Indels with the GATK Unified Genotyper
if [ -s \$OUTF/\$NAME.realigned.dm.recalibrated.bam ];
then
   echo -e \"\\n #### GATK: Call SNPs and Indels with the GATK Unified Genotyper \\n\"
   java  -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.both.raw.vcf -glm BOTH --downsampling_type NONE

   echo -e \"\\n #### GATK: Split SNPs and Indels \\n\"
   java  -jar \$GATK -T SelectVariants -R \$REF --variant \$OUTF/GATK.both.raw.vcf -o \$OUTF/GATK.snps.raw.vcf -selectType SNP
   java  -jar \$GATK -T SelectVariants -R \$REF --variant \$OUTF/GATK.both.raw.vcf -o \$OUTF/GATK.indel.raw.vcf -selectType INDEL
   
else
   echo \$NAME.realigned.dm.recalibrated.bam not found
   exit
fi

if [ ! -s \$OUTF/GATK.snps.raw.vcf ];
then
   echo GATK.snps.raw.vcf not found
   exit
fi




### Filter and compare SNP calls from 2 different pipelines
# Filtering
mkdir \$OUTF/SNP_Intersection

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf --variant \$OUTF/GATK.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 25.0 \" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg


# isolate PASSed variants
grep -E \'^#|PASS\' \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf | grep -v CRGg > \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf


# Correct sample names in VFC files
# sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
# sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-SHORE/\" \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf

#if [[ ! ( -s \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf && -s \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf ) ]];
#then
#   echo GATK.snps.filtered.cleaned.vcf or SHORE.snps.filtered.cleaned.vcf not found
#   exit
#fi

## Intersecting
#  java -jar \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:SHORE \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf -V:GATK \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -priority GATK,SHORE -o \$OUTF/SNP_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING

cp \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf \$OUTF/SNP_Intersection/merged.vcf

# Evaluation
java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP  -o \$OUTF/SNP_Intersection/report.all.txt --eval \$OUTF/SNP_Intersection/merged.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK

# Annotate Enrichment
\$BEDTOOLS/intersectBed -a \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -b \$EXOME > \$OUTF/SNP_Intersection/merged.all.vcf

# borrow header from GATK vcf file
grep '^#' \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf > \$OUTF/SNP_Intersection/merged.enriched.vcf
cat \$OUTF/SNP_Intersection/merged.all.vcf >> \$OUTF/SNP_Intersection/merged.enriched.vcf

java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP -o \$OUTF/SNP_Intersection/report.enriched.txt --eval \$OUTF/SNP_Intersection/merged.enriched.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK
      
	       ");

### decide which indel caller to choose and produce the respective parts in the pipeline
my @indelpipe = ();
my @indelanalysis = ();

### use clindel calls
if ($indelcaller eq 'clindel') {
	
	@indelpipe = ("

### SHORE: Prepare format map.list
 mkdir \$OUTF/shore
 \$CLINDEL convert --sort -r \$REF -n 6 -g 2 -e 100 Alignment2Maplist \$OUTF/\$NAME.realigned.dm.recalibrated.bam \$OUTF/shore/map.list.gz

if [ ! -s \$OUTF/shore/map.list.gz ];
then
   echo shore/map.list.gz not found
   exit
fi

### Clindel call indels
 \$CLINDEL qIndel -n \$NAME -f \$SHOREREF -o \$OUTF/shore/indels -i \$OUTF/shore/map.list.gz -k \$DBINDEL -c 3 -d 3 -C 400 -r 0.3 -a 0.20 -b 5

if [ ! -s \$OUTF/shore/indels ];
then
   echo shore/indels is empty
   exit
fi

### Clean up
set +e

rm -r \$OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
rm \$OUTF/shore/Variants/ConsensusAnalysis/reference.shore
rm \$OUTF/shore/map.list.gz

set -e

mkdir \$OUTF/Indel_Intersection


# select PASSed variants
grep -E \"^#|PASS\" \$OUTF/shore/indels/indels.vcf > \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf


# check for success of the preceeding steps
if [[ ! -s \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf ]];
then
   echo CLINDEL.indel.filtered.cleaned.vcf not found
   exit
fi
 
# Evaluation
java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBINDEL -select 'set==\"CLINDEL\"' -selectName CLINDEL -o \$OUTF/Indel_Intersection/report.all.txt --eval \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -l INFO

# filter Indels for enriched regions
\$BEDTOOLS/intersectBed -a \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -b \$EXOME > \$OUTF/Indel_Intersection/merged.all.vcf

# borrow header from GATK file
grep '^#' \$OUTF/GATK.indel.raw.vcf > \$OUTF/Indel_Intersection/merged.enriched.vcf
cat \$OUTF/Indel_Intersection/merged.all.vcf >> \$OUTF/Indel_Intersection/merged.enriched.vcf

java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBINDEL -select 'set==\"CLINDEL\"' -selectName CLINDEL -o \$OUTF/Indel_Intersection/report.enriched.txt --eval \$OUTF/Indel_Intersection/merged.enriched.vcf -l INFO

		      ");
	
}

### use GATK calls only
elsif ($indelcaller eq 'gatk') {
	
	@indelpipe = ("

mkdir \$OUTF/Indel_Intersection

java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant \$OUTF/GATK.indel.raw.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName LOWQ

# select PASSed variants
grep -v \"CRG\" \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf > \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf

# check for success of the preceeding steps
if [[ ! -s \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf ]];
then
   echo GATK.indel.filtered.cleaned.vcf not found
   exit
fi

# filter Indels for enriched regions
\$BEDTOOLS/intersectBed -a \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -b \$EXOME > \$OUTF/Indel_Intersection/merged.all.vcf

grep '^#' \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf > \$OUTF/Indel_Intersection/merged.enriched.vcf
cat \$OUTF/Indel_Intersection/merged.all.vcf >> \$OUTF/Indel_Intersection/merged.enriched.vcf

java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBINDEL -select 'set==\"GATK\"' -selectName GATK -o \$OUTF/Indel_Intersection/report.enriched.txt --eval \$OUTF/Indel_Intersection/merged.enriched.vcf -l INFO


		      
		      ");


}

### use GATK and CLINDEL calls
elsif ($indelcaller eq 'both') {
	
	@indelpipe = ("
### SHORE: Prepare format map.list
 mkdir \$OUTF/shore
 \$CLINDEL convert --sort -r \$REF -n 6 -g 2 -e 100 Alignment2Maplist \$OUTF/\$NAME.realigned.dm.recalibrated.bam \$OUTF/shore/map.list.gz

if [ ! -s \$OUTF/shore/map.list.gz ];
then
   echo shore/map.list.gz not found
   exit
fi

### Clindel call indels
 \$CLINDEL qIndel -n \$NAME -f \$SHOREREF -o \$OUTF/shore/indels -i \$OUTF/shore/map.list.gz -k \$DBINDEL -c 3 -d 3 -C 400 -r 0.3 -a 0.20 -b 5

if [ ! -s \$OUTF/shore/indels ];
then
   echo shore/indels is empty
   exit
fi

### Clean up
set +e

rm -r \$OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
rm \$OUTF/shore/Variants/ConsensusAnalysis/reference.shore
rm \$OUTF/shore/map.list.gz

set -e



mkdir \$OUTF/Indel_Intersection

# filter GATK variants (Clindel variants are produced pre-filtered)
java -jar -Xmx4g \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant \$OUTF/GATK.indel.raw.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName LOWQ

# select PASSed variants
grep -v \"CRG\" \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf > \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf

grep -E \"^#|PASS\" \$OUTF/shore/indels/indels.vcf > \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf



# Correct sample names in VFC files --- not necessary for tool switching
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-CLINDEL/\" \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf

# check for success of the preceeding steps
if [[ ! ( -s \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf && -s \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf ) ]];
then
   echo GATK.indel.filtered.cleaned.vcf or CLINDEL.indel.filtered.cleaned.vcf not found
   exit
fi
 
# Intersecting
java -Xmx5g -jar \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:GATK \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -V:CLINDEL \$OUTF/Indel_Intersection/CLINDEL.indel.filtered.cleaned.vcf -priority GATK,CLINDEL -o \$OUTF/Indel_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING

# Evaluation
java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBINDEL -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"CLINDEL\"' -selectName CLINDEL -select 'set==\"GATK\"'  -selectName GATK_CLINDEL -o \$OUTF/Indel_Intersection/report.all.txt --eval \$OUTF/Indel_Intersection/merged.vcf -l INFO

# filter Indels for enriched regions
\$BEDTOOLS/intersectBed -a \$OUTF/Indel_Intersection/merged.vcf -b \$EXOME > \$OUTF/Indel_Intersection/merged.all.vcf

grep '^#' \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf > \$OUTF/Indel_Intersection/merged.enriched.vcf
cat \$OUTF/Indel_Intersection/merged.all.vcf >> \$OUTF/Indel_Intersection/merged.enriched.vcf

java -Xmx5g -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBINDEL -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"CLINDEL\"' -selectName CLINDEL -select 'set==\"GATK\"' -selectName GATK -o \$OUTF/Indel_Intersection/report.enriched.txt --eval \$OUTF/Indel_Intersection/merged.enriched.vcf -l INFO


		  ");

}
else {
	print 'Option provided to --indelcaller should be: "clindel", "gatk" or "both" ', "\n";
	exit;
}




### loop through each sample and create a respective shell script
my @files = glob($infolder . "/*$firstreadextension");

foreach my $read1 (@files) {

	# read 2 filename
	my $read2 = $read1;
	$read2 =~ s/$firstreadextension/$secondreadextension/ge;
	
	my @filepath = split("/", $read1);
	my $fileleaf = $filepath[$#filepath];


	my $name = substr($fileleaf, $nameStart - 1 , $nameLength);

	unless (-e "$outfolder") {
		mkdir "$outfolder"  or die "Cannot create output directory $outfolder";
	}
	
	if(! -e "$outfolder/$name") {
		mkdir "$outfolder/$name" or die "Cannot create output directory $outfolder/$name";
	}
	else {
		print STDERR "Runfolder $outfolder/$name already exists. Will only update qsub file\n"
	}
	
	open OUT, ">$outfolder/$name/$qsub_name" or die "Cannot create qsub file $outfolder/$name/$qsub_name";
	
	my $qsubNname; 
	if ($name =~ /^\d/) {
		$qsubNname = 's'.$name;
	}
	else {
		$qsubNname = $name;
	}# qsub doesn't like jobs starting with a digit

	my @qsub = ("#!/bin/bash

set -e

#\$ -N $qsubNname
#\$ -e $outfolder/$name/
#\$ -o $outfolder/$name/
#\$ -pe smp $cpu
#\$ -l virtual_free=${mem}G

source /etc/profile
export _JAVA_OPTIONS=\"-Djava.io.tmpdir=\$TMPDIR \$_JAVA_OPTIONS\"


NAME=$name
READ1=$read1
READ2=$read2
OUTF=$outfolder/$name


@snppipe

@indelpipe



\n");



	print OUT @qsub;

	close OUT;
};
