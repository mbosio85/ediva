#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

sub usage { print "\n$0 \n usage:\n",
	   "--vcffile \t vcf file in the format NAME:file.vcf \n",
           "--bamfile \t bam file \n",
	   "--outfolder \t folder to write the mpileup result to \n",
	   "--qsubname \t name of the script you want to start lateron \n",
	   "--help \t\t show help \n";
}


my @infiles;
my @bamfiles;
my $outfolder;
my $qsub_name;
my $help = 0;

GetOptions("vcffile=s" => \@infiles, "bamfile=s" => \@bamfiles, "outfolder=s" => \$outfolder, "qsubname=s" => \$qsub_name, "help=s" => \$help);

unless (@infiles && @bamfiles && $outfolder && $qsub_name && $help == 0) {
    usage;
    exit;
}

# create output folder, if it doesn't exist
unless ( -e $outfolder) {
    mkdir $outfolder;
}

my @samples;
my @infile_paths;

# extract information of infiles and prepare header of output
foreach my $line (@infiles) {
    my @splitline = split(':', $line);
    push @samples, $splitline[0]."obs";
    push @samples, $splitline[0]."qual";
    push @infile_paths, $splitline[1];
}
my $header = join('\\t', @samples);

# will contain vcf --variant vcf --variant vcf ..., i.e. a --variant call in GATK is still necessary
my $vcf_files = join(' --variant ', @infile_paths); 

# write a temporary list of bam files for mpileup input
my $bam_file_list = join("\n", @bamfiles);

open FH, ">$outfolder/bam.list" or die "Could not create file containing the paths to all bam files $outfolder/bam.list";
print FH $bam_file_list;
close FH;

if ($qsub_name =~ /^\d/) {
        $qsub_name = 's'.$qsub_name;
}
else {
        $qsub_name = $qsub_name;
}# qsub doesn't like jobs starting with a digit


my @OUT = ("

#\$ -N $qsub_name
#\$ -e $outfolder/
#\$ -o $outfolder/

GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
SAMTOOLS=\$(which samtools)

REF=/users/GD/resource/human/hg19/bwa7/hg19.fasta
OUTF=$outfolder

# merge vcf files
java -Xmx2g -jar \$GATK -T CombineVariants -R \$REF --variant $vcf_files -o \$OUTF/combined.variants.vcf 

# pileup call
sed -e 's/chr//' \$OUTF/combined.variants.vcf | awk '{OFS=\"\\t\"; if (!/^#/){print \$1,\$2}}'  > \$OUTF/variant.position.bed

# header
echo -e \"chr\\tstart\\tref\\tcount\\t$header\" > \$OUTF/variant.positions.mpileup
# pileup
\$SAMTOOLS mpileup -l \$OUTF/variant.position.bed -f \$REF -b \$OUTF/bam.list >> \$OUTF/variant.positions.mpileup

");

open FH, ">$outfolder/$qsub_name" or die "Could not open $outfolder/$qsub_name for writing.";
print FH @OUT;
close FH;