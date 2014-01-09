#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
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
#
#  Module: Support::qsub::qsub_batch_exome_fastq_PE.pl
#  Purpose:
#  In:
#  Out:
#

sub usage { print "\n$0 \n usage:\n",
	   "--infolder \t folder containing the raw sequence data in fastq format \n",
	   "--outfolder \t folder to write the analysis to \n",
	   "--qsubname \t name of the script you want to start lateron \n",
	   "--max_coverage \t used in SNP filtering with samtools [default = 400] \n",
	   "--namestart \t  start of a substring in the read file name (first letter is numbered 1) \n",
	   "--namelength \t length of the part of the filenames which should be taken as sample (and folder) name \n",
	   "--firstreadextension \t describes the name of the first read file (e.g. 1.fastq.gz or 1_sequence.txt.gz)\n",
	   "--secondreadextension \t describes the name of the first read file (e.g. 2.fastq.gz or 2_sequence.txt.gz)\n",
	   "--enrichment \t which enrichment kit was used? (Agilent: a35 - 35MB version 1, a50 - 50MB version 3, a75 - 75MB version 5 + UTR, Nimblegene: n35 - 35MB version 2, n50 - 50MB version 3)\n",
	   "--cpu \t number of cpu cores to be used (applicable only for a few steps) [default = 4]\n",
	   "--mem \t amount of Memory dedicated to your job in Gb. The amount of memory must not be bigger than available in the machine. [default = 12]\n",
	   ##"--aligner \t which alignment program should be used. GEM or BWA [default = BWA]\n",
	   "--minimalannotation \t use AnnoVar to get a minimal annotation containing only gene and amino acid change\n",
	   "--help \t\t show help \n";
}


my $infolder;
my $outfolder;
my $qsub_name;
my $max_cov = 400; 
my $nameStart = 'NA';
my $nameLength;
my $firstreadextension;
my $secondreadextension;
my $enrichment;
my $cpu = 4;
my $mem = 12;
##my $aligner = 'BWA';
my $minimal = 0;
my $help = 0;

my %enrichmentkits = ('a35' => '/users/GD/resource/human/probesets/agilent/35MB_standard/shore_format',
		      #'a35e' => '/users/GD/resource/human/probesets/agilent/35MB_extended/shore_format',
		      'a50' => '/users/GD/resource/human/probesets/agilent/50MB/shore_format',
		      'a75' => '/users/GD/resource/human/probesets/agilent/75MB/shore',
		      'n35' => '/users/GD/resource/human/probesets/nimblegene/v2/Target_Regions/shore_format',
		      'n50' => '/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/shore_format/');

my %coniferprobes =  ('a35' => '/users/GD/resource/human/probesets/agilent/35MB_standard/conifer/Exome_Array_plus50.bed.tiled.conifer.withoutXY',
		      #'a35e' => '',
		      'a50' => '/users/GD/resource/human/probesets/agilent/50MB/shore_format/Exome_Array_plus50.bed.tiled.conifer.noxy',
		      'a75' => '/users/GD/resource/human/probesets/agilent/75MB/S04380219_Regions.nohead.bed_plus_50_shore.tiled.conifer.noxy',
		      'n35' => '/users/GD/resource/human/probesets/nimblegene/v2/Target_Regions/SeqCap_EZ_Exome_v2_TiledRegions.bed_plus_50_shore.tiled.conifer.noxy',
		      'n50' => '/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50_shore.tiled.conifer.noxy');

GetOptions("infolder=s" => \$infolder, "outfolder=s" => \$outfolder, "qsubname=s" => \$qsub_name, "max_coverage=i" => \$max_cov, "nameStart=s" => \$nameStart, "nameLength=s" => \$nameLength, "firstreadextension=s" => \$firstreadextension, "secondreadextension=s" => \$secondreadextension, "cpu=i" => \$cpu, "enrichment=s" => \$enrichment, "mem=i" => \$mem, "minimalannotation" => \$minimal, "help=s" => \$help);

unless($infolder && $outfolder && $qsub_name && $nameStart ne 'NA' && $nameLength && $firstreadextension && $secondreadextension && $enrichment && $help == 0) {
	usage;
	exit;
}

unless(defined($enrichmentkits{$enrichment})) {
	print "\n###########################\n enrichment method unknown \n###########################\n";
	usage;
	exit;
}

unless ($minimal == 0) {
	$minimal = '--step 1';
}
else {
	$minimal = '';
}

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
# export TMPDIR=/users/GD/projects/HumanDisease/tmp # inactive for new cluster
# export _JAVA_OPTIONS=-Djava.io.tmpdir=/users/GD/projects/HumanDisease/tmp # changed for new cluster
export _JAVA_OPTIONS=-Djava.io.tmpdir=\"-Djava.io.tmpdir=\$TMPDIR \$_JAVA_OPTIONS\"
#export PATH=/users/GD/tools/annovar/annovar_2011Nov20/:\$PATH
export PATH=/users/GD/tools/annovar/annovar_2013May09/:\$PATH


NAME=$name
READ1=$read1
READ2=$read2
OUTF=$outfolder/$name

REF=/users/GD/resource/human/hg19/bwa7/hg19.fasta
SHOREREF=/users/GD/resource/human/hg19/shore/hg19.fasta.shore

EXOME=$enrichmentkits{$enrichment}
PROBE=$coniferprobes{$enrichment}
INDELPRIOR=/users/GD/tools/clindel/resources/dbindel137_121217.pseudovcf
DBINDEL=/users/GD/resource/human/hg19/databases/dbSNP/dbindel137_121217.vcf
DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp137_121217.vcf

BWA=/users/GD/tools/bwa/bwa-0.7.5a/bwa
CONIFER='/software/so/el6.3/PythonPackages-2.7.3-VirtualEnv/bin/python /users/GD/tools/conifer/conifer_v0.2.2/conifer.py'

GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.4-9-g532efad/GenomeAnalysisTK.jar
SAMTOOLS=\$(which samtools)
BCFTOOLS=\$(which bcftools)
VCFUTILS=\$(which vcfutils.pl)
ANNOVAR=/users/GD/tools/annovar/annovar_2013May09/
SHORE=/users/GD/tools/shore/shore
CLINDEL=/users/GD/tools/clindel/bin/shore
NGSBOX=/users/GD/tools/ngsbox
RSCRIPT=\$(which Rscript)
SPLITGATK='perl /users/GD/tools/exomeCRG/splitGATKoutput.pl'

ANNOVARSTEPS='$minimal'



### set a quality threshold for Variant Filtration
QD='|| QD < 4.0' # || HRun > 9' is not set any more in GATK-2.4.9

### Align reads with bwa
\$BWA mem -M -t $cpu -R \"\@RG\\tID:\$NAME\\tSM:\$NAME\" \$REF \$READ1 \$READ2 | time \$SAMTOOLS view -h -b -S -F 0x900 -  > \$TMPDIR/\$NAME.noChimeric.bam

OFFSET=\$(\$SAMTOOLS view \$TMPDIR/\$NAME.noChimeric.bam | /software/so/el6.3/PythonPackages-2.7.3-VirtualEnv/bin/python /users/GD/tools/exomeCRG/whichQuality_bam.py)
if [[ \$OFFSET == 64 ]];
then
	echo \"fixing 64 quality encoding\"
	\$SAMTOOLS view -h \$TMPDIR/\$NAME.noChimeric.bam | /software/so/el6.3/PythonPackages-2.7.3-VirtualEnv/bin/python /users/so/odrechsel/scripts/bam_rescale_quals.py - | \$SAMTOOLS view -bS - > \$TMPDIR/\$NAME.transformed.bam
	rm \$TMPDIR/\$NAME.noChimeric.bam
	mv \$TMPDIR/\$NAME.transformed.bam \$TMPDIR/\$NAME.noChimeric.bam
fi


### Sort BAM file
if [ -s \$TMPDIR/\$NAME.noChimeric.bam ];
then
   echo Sort BAM
   /users/GD/tools/novocraft/novosort/novosort --threads $cpu --tmpdir \$TMPDIR --forcesort --output \$TMPDIR/\$NAME.sort.bam -i -m ${mem}G \$TMPDIR/\$NAME.noChimeric.bam
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
   java -jar \$GATK -nt $cpu -T RealignerTargetCreator -R \$REF -I \$TMPDIR/\$NAME.sort.bam -o \$OUTF/\$NAME.intervals -known \$DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE \$FMQ
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
   echo -e \"\\n #### duplicate marking using: java -jar /users/GD/tools/picard/picard-tools-1.100/MarkDuplicates.jar INPUT=\$TMPDIR/\$NAME.realigned.bam OUTPUT=\$TMPDIR/\$NAME.realigned.dm.bam METRICS_FILE=\$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT \\n\"
   java -jar /users/GD/tools/picard/picard-tools-1.100/MarkDuplicates.jar INPUT=\$TMPDIR/\$NAME.realigned.bam OUTPUT=\$TMPDIR/\$NAME.realigned.dm.bam METRICS_FILE=\$OUTF/duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
   \$SAMTOOLS index \$TMPDIR/\$NAME.realigned.dm.bam
else
   echo \$TMPDIR/\$NAME.realigned.bam not found
   exit
fi



### Base quality recalibration
if [ -s \$TMPDIR/\$NAME.realigned.dm.bam ];
then
   echo -e \" \\n #### Base quality recalibration \\n \"
   java -jar \$GATK -T BaseRecalibrator --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R \$REF -I \$TMPDIR/\$NAME.realigned.dm.bam -knownSites \$DBSNP --downsampling_type NONE -o \$TMPDIR/\$NAME.recal_data.grp
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
	
	rm \$TMPDIR/*map
	rm \$TMPDIR/*map.gz
	rm \$TMPDIR/\$NAME.sam
	rm \$TMPDIR/\$NAME.bam
	rm \$TMPDIR/\$NAME.realigned.bam
	rm \$TMPDIR/\$NAME.realigned.bai
	rm \$TMPDIR/\$NAME.realigned.dm.bam
	rm \$TMPDIR/\$NAME.realigned.dm.bam.bai
set -e


set +e

if [[ ( -s \$OUTF/\$NAME.realigned.dm.recalibrated.bam  ) ]]
then
	echo -e \"\\n calculating rpkm \\n \"
	ln -s \$OUTF/\$NAME.realigned.dm.recalibrated.bai \$OUTF/\$NAME.realigned.dm.recalibrated.bam.bai
	\$CONIFER rpkm --probes \$PROBE --input \$OUTF/\$NAME.realigned.dm.recalibrated.bam --output \$OUTF/\$NAME.rpkm.txt
	
	if [ \$\? -ne 0 ]; then
		echo -e \"\\n rpkm calculation FAILed, but whole pipeline will keep running \\n \"
	else
		echo -n \"\$OUTF/\$NAME.rpkm.txt ### \" >> /users/GD/projects/rpkm/allRPKM.log
		echo -n \"\$(whoami)\" >> /users/GD/projects/rpkm/allRPKM.log
		date >> /users/GD/projects/rpkm/allRPKM.log
	fi
	
else
	echo -e \'\\n \$OUTF/\$NAME.realigned.dm.recalibrated.bam not found for Conifer  \\n \'
fi

set -e


### GATK: Call SNPs and Indels with the GATK Unified Genotyper
if [ -s \$OUTF/\$NAME.realigned.dm.recalibrated.bam ];
then
   echo -e \"\\n #### GATK: Call SNPs and Indels with the GATK Unified Genotyper \\n\"
   java  -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.both.raw.vcf -glm BOTH --downsampling_type NONE
   #java -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.snps.raw.vcf -glm SNP --downsampling_type NONE
   #java -jar \$GATK -T UnifiedGenotyper -nt $cpu -R \$REF -I \$OUTF/\$NAME.realigned.dm.recalibrated.bam -o \$OUTF/GATK.indel.raw.vcf -glm INDEL --downsampling_type NONE
   \$SPLITGATK \$OUTF \$OUTF/GATK.both.raw.vcf
else
   echo \$NAME.realigned.dm.recalibrated.bam not found
   exit
fi

if [ ! -s \$OUTF/GATK.snps.raw.vcf ];
then
   echo GATK.snps.raw.vcf not found
   exit
fi



### MPILEUP: Call SNPs and Indels
   echo -e \"\\n #### MPILEUP: call SNPs and InDels \\n \"
   \$SAMTOOLS mpileup -uf \$REF \$OUTF/\$NAME.realigned.dm.recalibrated.bam -L 2500 | \$BCFTOOLS view -bcg - > \$OUTF/MPILEUP.variant.raw.bcf
   if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi   
   \$BCFTOOLS view \$OUTF/MPILEUP.variant.raw.bcf | \$VCFUTILS varFilter -d5 -D$max_cov -W 20 > \$OUTF/MPILEUP.variant.raw.vcf
   if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi
   egrep \"INDEL|#\" \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.indel.raw.vcf
   grep -v INDEL \$OUTF/MPILEUP.variant.raw.vcf > \$OUTF/MPILEUP.snps.raw.vcf

if [ ! -s \$OUTF/MPILEUP.snps.raw.vcf ];
then
   echo MPILEUP.snps.raw.vcf not found
   exit
fi



### SHORE: Prepare format map.list
 mkdir \$OUTF/shore
 \$SHORE convert --sort -r \$REF -n 6 -g 1 -e 20 -s Alignment2Maplist \$OUTF/\$NAME.realigned.dm.recalibrated.bam \$OUTF/shore/map.list.gz

if [ ! -s \$OUTF/shore/map.list.gz ];
then
   echo shore/map.list.gz not found
   exit
fi



### SHORE: compute coverage plot in GFF format for browsers
 \$SHORE coverage -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/CoverageAnalysis



### SHORE: Compute enrichment
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus200 -f \$EXOME/Exome_Array_plus200.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus150 -f \$EXOME/Exome_Array_plus150.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus100 -f \$EXOME/Exome_Array_plus100.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus50 -f \$EXOME/Exome_Array_plus50.bed -H 1,1 -k
 \$SHORE count -m \$OUTF/shore/map.list.gz -o \$OUTF/shore/Count_SureSelect_plus0 -f \$EXOME/Exome_Array_plus0.bed -H 1,1 -k



### SHORE: Enrichment plots
 grep enriched \$OUTF/shore/Count_SureSelect_plus150/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus150/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_depleted.txt
 grep enriched \$OUTF/shore/Count_SureSelect_plus150/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_count_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus150/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus150/exome_count_depleted.txt
### plot data
 \$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus150/

 grep enriched \$OUTF/shore/Count_SureSelect_plus0/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus0/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_depleted.txt
 grep enriched \$OUTF/shore/Count_SureSelect_plus0/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_count_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus0/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus0/exome_count_depleted.txt
### plot data
 \$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus0/

 grep enriched \$OUTF/shore/Count_SureSelect_plus50/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus50/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_depleted.txt
 grep enriched \$OUTF/shore/Count_SureSelect_plus50/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_count_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus50/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus50/exome_count_depleted.txt
### plot data
 \$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus50/

 grep enriched \$OUTF/shore/Count_SureSelect_plus100/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus100/meancov.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_depleted.txt
 grep enriched \$OUTF/shore/Count_SureSelect_plus100/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_count_enriched.txt
 grep depleted \$OUTF/shore/Count_SureSelect_plus100/readcount.txt | cut -f5 > \$OUTF/shore/Count_SureSelect_plus100/exome_count_depleted.txt
### plot data
 \$RSCRIPT \$NGSBOX/Statistics/R_examples/exome_enrichment_stats.R \$OUTF/shore/Count_SureSelect_plus100/



### SHORE: Call SNPs and Indels
 \$SHORE qVar -n \$NAME -f \$SHOREREF -o \$OUTF/shore/Variants -i \$OUTF/shore/map.list.gz -s /users/GD/tools/shore/scoring_matrices/scoring_matrix_het.txt -E \$OUTF/shore/Count_SureSelect_plus150/meancov.txt -e -c 4 -d 4 -C $max_cov -r 3 -q 10 -Q 15 -a 0.25 -b 6 -y -v

if [ ! -s \$OUTF/shore/Variants/ConsensusAnalysis ];
then
   echo shore/Variants/ConsensusAnalysis is empty
   exit
fi


### Clindel
 \$CLINDEL qIndel -n \$NAME -f \$SHOREREF -o \$OUTF/shore/indels -i \$OUTF/shore/map.list.gz -k \$INDELPRIOR -c 2 -d 2 -C $max_cov -r 3 -a 0.20 -b 5 -y -w -v

if [ ! -s \$OUTF/shore/indels ];
then
   echo shore/indels is empty
   exit
fi


### Clean up
rm -r \$OUTF/shore/Variants/ConsensusAnalysis/supplementary_data
gzip -9 \$OUTF/shore/Variants/ConsensusAnalysis/reference.shore




### Filter and compare SNP calls from 3 different pipelines
# Filtering
mkdir \$OUTF/SNP_Intersection

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf --variant \$OUTF/GATK.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 25.0 \$QD \" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.vcf --variant \$OUTF/MPILEUP.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 15.0 || DP < 5 || DP > $max_cov\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter.pl \$OUTF/shore/Variants/ConsensusAnalysis/snp.vcf $max_cov > \$OUTF/SHORE.snps.raw.vcf

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/SNP_Intersection/SHORE.snps.filtered.vcf --variant \$OUTF/SHORE.snps.raw.vcf --mask \$OUTF/GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"QUAL < 20.0 || DP < 5 || DP > $max_cov\" --filterName CRG -U LENIENT_VCF_PROCESSING

# greping
grep -v \"CRG\" \$OUTF/SNP_Intersection/GATK.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi
grep -v \"CRG\" \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi
grep -v \"CRG\" \$OUTF/SNP_Intersection/SHORE.snps.filtered.vcf | grep -v \"SnpCluster\" > \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi

# Correct sample names in VFC files
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-MPILEUP/\" \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-SHORE/\" \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf

if [[ ! ( -s \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf && -s \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf && -s \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf ) ]];
then
   echo GATK.snps.filtered.cleaned.vcf or MPILEUP.snps.filtered.cleaned.vcf or SHORE.snps.filtered.cleaned.vcf not found
   exit
fi

# Intersecting
java -jar \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:SHORE \$OUTF/SNP_Intersection/SHORE.snps.filtered.cleaned.vcf -V:GATK \$OUTF/SNP_Intersection/GATK.snps.filtered.cleaned.vcf -V:MPILEUP \$OUTF/SNP_Intersection/MPILEUP.snps.filtered.cleaned.vcf -priority GATK,MPILEUP,SHORE -o \$OUTF/SNP_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING

# Evaluation
java -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/SNP_Intersection/report.all.txt --eval \$OUTF/SNP_Intersection/merged.vcf -l INFO

# Annotate Enrichment
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$EXOME/Exome_Array_plus150.bed \$OUTF/SNP_Intersection/merged.vcf > \$OUTF/SNP_Intersection/merged.all.vcf

# Evaluate calls on enriched regions
grep -v \"NOTENRICHED\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.enriched.vcf

java -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/SNP_Intersection/report.enriched.txt --eval \$OUTF/SNP_Intersection/merged.enriched.vcf -l INFO




### Filter and compare indel calls from 3 different pipelines
# Filtering
mkdir \$OUTF/Indel_Intersection

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf --variant \$OUTF/GATK.indel.raw.vcf --filterExpression \"MQ < 30.0 || QUAL < 20.0 || MQ0 > 5 \$QD \" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.vcf --variant \$OUTF/MPILEUP.indel.raw.vcf --filterExpression \"MQ < 30.0 || QUAL < 10.0 || DP < 5 || DP > $max_cov\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > $max_cov || GQ < 15\" --genotypeFilterName CRGg

java -jar \$GATK -T VariantFiltration -R \$REF -o \$OUTF/Indel_Intersection/SHORE.indel.filtered.vcf --variant \$OUTF/shore/Variants/ConsensusAnalysis/indels.vcf --filterExpression \"QUAL < 2.0 || DP < 4 || DP > $max_cov || RE > 1.3\" --filterName CRG  -U LENIENT_VCF_PROCESSING

# greping
grep -v \"CRG\" \$OUTF/Indel_Intersection/GATK.indel.filtered.vcf > \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi
grep -v \"CRG\" \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.vcf > \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi
grep  -v \"CRG\" \$OUTF/Indel_Intersection/SHORE.indel.filtered.vcf | grep -v \"SHOREFILTER\" > \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf
if [ \${PIPESTATUS[0]} -ne 0 -a \${PIPESTATUS[1]} -ne 0 ] ; then echo FAIL \$?; fi

# Correct sample names in VFC files
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-GATK/\" \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-MPILEUP/\" \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf
sed -i -e \"s/FORMAT\\t\$NAME/FORMAT\\t\$NAME-SHORE/\" \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf

if [[ ! ( -s \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf && -s \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf && -s \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf ) ]];
then
   echo GATK.indel.filtered.cleaned.vcf or MPILEUP.indel.filtered.cleaned.vcf or SHORE.indel.filtered.cleaned.vcf not found
   exit
fi
 
# Intersecting
java -jar \$GATK -T CombineVariants -R \$REF -genotypeMergeOptions PRIORITIZE -V:SHORE \$OUTF/Indel_Intersection/SHORE.indel.filtered.cleaned.vcf -V:GATK \$OUTF/Indel_Intersection/GATK.indel.filtered.cleaned.vcf -V:MPILEUP \$OUTF/Indel_Intersection/MPILEUP.indel.filtered.cleaned.vcf -priority GATK,MPILEUP,SHORE -o \$OUTF/Indel_Intersection/merged.vcf -U LENIENT_VCF_PROCESSING

# Evaluation
java -jar \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/Indel_Intersection/report.all.txt --eval \$OUTF/Indel_Intersection/merged.vcf -l INFO

# Annotate Enrichment
perl \$NGSBOX/Parser/VCF/vcf_filter/vcf_filter_enriched.pl \$EXOME/Exome_Array_plus150.bed \$OUTF/Indel_Intersection/merged.vcf > \$OUTF/Indel_Intersection/merged.all.vcf

# Evaluate calls on enriched regions 
grep -v \"NOTENRICHED\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.enriched.vcf

java -jar -Xmx4g \$GATK -T VariantEval -R \$REF --dbsnp \$DBSNP -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"SHORE\"' -selectName SHORE -select 'set==\"MPILEUP\"' -selectName MPILEUP -select 'set==\"GATK\"' -selectName GATK -select 'set==\"GATK-MPILEUP\"' -selectName GATK_MPILEUP -select 'set==\"GATK-SHORE\"' -selectName GATK_SHORE -select 'set==\"MPILEUP-SHORE\"' -selectName MPILEUP_SHORE -o \$OUTF/Indel_Intersection/report.enriched.txt --eval \$OUTF/Indel_Intersection/merged.enriched.vcf -l INFO




### Annotate SNPs with ANNOVAR: Intersection, all three tools predict SNP
mkdir \$OUTF/SNP_Intersection/AnnovarIntersection
egrep \"Intersection|#\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.intersection.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.intersection.vcf > \$OUTF/SNP_Intersection/AnnovarIntersection/snps.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/SNP_Intersection/AnnovarIntersection/sum \$OUTF/SNP_Intersection/AnnovarIntersection/snps.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/


### Annotate SNPs with ANNOVAR: Partial union, at least 2 tools predict SNP
mkdir \$OUTF/SNP_Intersection/AnnovarPartialUnion
egrep \"Intersection|GATK-SHORE|MPILEUP-SHORE|GATK-MPILEUP|#\" \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.union.vcf > \$OUTF/SNP_Intersection/AnnovarPartialUnion/snps.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/SNP_Intersection/AnnovarPartialUnion/sum \$OUTF/SNP_Intersection/AnnovarPartialUnion/snps.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/


### Annotate SNPs with ANNOVAR: Union, any tool predicts SNP
mkdir \$OUTF/SNP_Intersection/AnnovarUnion
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/SNP_Intersection/merged.all.vcf > \$OUTF/SNP_Intersection/AnnovarUnion/snps.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/SNP_Intersection/AnnovarUnion/sum \$OUTF/SNP_Intersection/AnnovarUnion/snps.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/



### Annotate Indels with ANNOVAR: Intersection
mkdir \$OUTF/Indel_Intersection/AnnovarIntersection
egrep \"Intersection|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.union.vcf > \$OUTF/Indel_Intersection/AnnovarIntersection/indels.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/Indel_Intersection/AnnovarIntersection/sum \$OUTF/Indel_Intersection/AnnovarIntersection/indels.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/


### Annotate Indels with ANNOVAR: Partial Union, must include SHORE or at least two tools
mkdir \$OUTF/Indel_Intersection/AnnovarUnionShore
# egrep \"Intersection|GATK-SHORE|MPILEUP-SHORE|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
egrep \"Intersection|GATK-MPILEUP|SHORE|#\" \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/merged.union.vcf
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.union.vcf > \$OUTF/Indel_Intersection/AnnovarUnionShore/indels.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/Indel_Intersection/AnnovarUnionShore/sum \$OUTF/Indel_Intersection/AnnovarUnionShore/indels.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/


### Annotate Indels with ANNOVAR: Union, indels predicted by any tool
mkdir \$OUTF/Indel_Intersection/AnnovarUnion
 \$ANNOVAR/convert2annovar.pl -format vcf4 \$OUTF/Indel_Intersection/merged.all.vcf > \$OUTF/Indel_Intersection/AnnovarUnion/indels.ann
 \$ANNOVAR/summarize_annovar.pl --buildver hg19 \$ANNOVARSTEPS --outfile \$OUTF/Indel_Intersection/AnnovarUnion/sum \$OUTF/Indel_Intersection/AnnovarUnion/indels.ann --ver1000g 1000g2012apr --verdbsnp 137 --veresp 6500si \$ANNOVAR/hg19/



### Clean up
rm \$OUTF/MPILEUP.variant.raw.bcf
rm \$OUTF/shore/map.list.gz

\n");



	print OUT @qsub;

	close OUT;
}
