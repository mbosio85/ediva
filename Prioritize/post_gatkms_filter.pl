#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


######################################
#
#	Task => Keep any type of variants with atleast 5 reads (default value) in atleast one sample from the total sample set in the VCF; others are dropped  
#	infile => GATK multi-sample VCF file with complete sample wise genoptype information
#	outfile => Filtered GATK multi-sample VCF file with complete sample wise genotype information
#	Extra outfile => GATK multi-sample VCF file with complete sample wise genotype information with the filtered out variants
#
#######################################

sub usage { print "\n$0 \n usage: \n",
	   "--infile \t multisample VCF file \n",
	   "--outfile \t post filtered VCF file \n",
	   "--value \t minimum read depth supporting the variant in atleast one sample (default : 5) \n",
	   "--help \t\t show help \n";
}

my $infile;
my $outfile;
my $value = 5;
my $help = 0;

GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "value=s" => \$value,"help=s" => \$help);

unless($infile && $outfile && $help == 0) 
{
	usage;
	exit;
}

## create path for output file and file names
my $exfile;

if ($outfile =~ m/\//) ## file with full path
{
	my @paths = split(/\//,$infile);
	my $len = scalar @paths;
	my $outpath = join("/",@paths[0..($len-2)]);

	my @files = split(/\./,$paths[$len-1]);
	$len = scalar @files;
	my $tempfile = join(".",@files[0..($len-2)]);
	$exfile = $outpath."/".$tempfile.".lowQual.vcf";

}else{ ## just filename in the same script location

	my @files = split(/\./,$outfile);
	my $len = scalar @files;
	my $tempfile = join(".",@files[0..($len-2)]);
	$exfile = $tempfile.".lowQual.vcf";
}

## open handlers
open (NEW,">>".$outfile) or die "Cant open new file \n";
open (BAD,">>".$exfile) or die "Cant open new file \n";

open (VCF, $infile) or die "Cant open $infile \n";

## start browsing data
while (<VCF>)
{
	chomp $_;
	## flag to make decision at the end for the variant 0 => NO , 1 => YES
	my $flag = 0;
	
	if ($_ =~ m/^#/) ## header line
	{
		## print the header lines for both the output VCF files
		print NEW $_."\n";
		print BAD $_."\n";
	}else ## data line
	{
		my @data = split(/\t/,$_);
		
		## iterate over all the samples for this variant in the VCF file; we start with index 9 in the data line cause thats the point
		## where the sample information starts
		for (my $i = 9; $i < @data ; $i++)
		{
			## data field format and value :: GT:AD:DP:GQ:PL => 0/0:44,0:44:99:0,117,1274
			my @sampleInfo = split(/\:/,$data[$i]);
			
			## only check when zygosity is not 0/0 or ./.
			if ($sampleInfo[0] ne '0/0' and $sampleInfo[0] ne './.')
			{
				my ($refDP,$altDP) = split(/\,/,$sampleInfo[1]);
				
				if ($altDP >= $value)
				{
					$flag = 1;
					last; ## we dont need to check anymore for this variant
				}
			}
		}
		if ($flag == 1)
		{
			print NEW $_."\n";
		}else{
			print BAD $_."\n";
		}	
	}
	
}

## close handlers
close(NEW);
close(BAD);
close(VCF);
