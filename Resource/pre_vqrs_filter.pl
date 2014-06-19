#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


######################################
#
#	Task => performs filtering on sample genotypes and on variants 
#	infile => GATK multi-sample VCF file with complete sample wise genoptype information
#	outfile => Filtered GATK multi-sample VCF file with complete sample wise genotype information
#	Extra outfile => GATK multi-sample VCF file with complete sample wise genotype information with the filtered out variants
#
#######################################

sub usage { print "\n$0 \n usage: \n",
	   "--infile \t multisample VCF file \n",
	   "--outfile \t post filtered VCF file \n",
	   "--ad_alt_allele \t minimum allelic depth supporting the variant in the sample genotype (default : 3) \n",
           "--dp \t minimum total read depth in the sample genotype (default : 8) \n",
           "--gq \t minimum genotype quality in the sample genotype (default : 10) \n",
           "--avg_gq \t  minimum average genotype quality for the variant (default : 20) \n",
           "--call_rate \t  minimum call rate in sample set (default : 80\%) \n",
           "--min_avg_af \t  minimum average allele fraction in minimum 30 sample genotypes (default : 0.32) \n",
           "--max_avg_af \t  maximum average allele fraction in minimum 30 sample genotypes (default : 0.68) \n",
	   "--help \t\t show help \n";
}

## variables
my $infile;
my $outfile;
my $ad_alt_allele = 3;
my $dp = 8;
my $gq = 10;
my $avg_gq = 20;
my $call_rate = 80;
my $min_avg_af = 0.32;
my $max_avg_af = 0.68;
my $help = 0;

## get parameters
GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "ad_alt_allele=s" => \$ad_alt_allele, "dp=s" => \$dp, "gq=s" => \$gq, "avg_gq=s" => \$avg_gq, "call_rate=s" => \$call_rate, "min_avg_af=s" => \$min_avg_af, "max_avg_af=s" => \$max_avg_af, "help=s" => \$help);

## check on necessary parameters
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

## delete existing output files and open handlers
unlink($outfile);
unlink($exfile);
open (NEW,">>".$outfile) or die "Cant open new file \n";
open (BAD,">>".$exfile) or die "Cant open new file \n";

open (VCF, $infile) or die "Cant open $infile \n";

## start browsing data
while (<VCF>)
{
	my @data = ();
        my $allelefrac = "NA";
        my $sample_count_for_allelefrac = 0;
	my $sample_count = 0;
	my $genotype_count = 0;
	my $samples = 0;
	my $seen_call_rate = "NA";
	my $seen_avg_gq = 0;	

	if ($_ =~ m/^#/) 
	{
		## vcf header line
		## print the header lines for both the output VCF files
		print NEW $_;
		print BAD $_;
	}else 
	{	
		## data line
		chomp $_;
		## split data line
		@data = split(/\t/,$_);
		

		## iterate over all the samples for this variant in the VCF file; we start with index 9 in the data line cause thats the point
		## where the sample information starts
		for (my $i = 9; $i < @data ; $i++)
		{
			## data field format and value :: GT:AD:DP:GQ:PL => 0/0:44,0:44:99:0,117,1274
			my @sampleInfo = split(/\:/,$data[$i]);
			
			## only check when zygosity is not ./. and 0/0
			if ($sampleInfo[0] ne './.')
			{
				if ($sampleInfo[0] ne '0/0')
				{
					## only variant sites
					$sample_count = $sample_count + 1;
					
					my ($refDP,$altDP) = ("NA","NA");
					if ($sampleInfo[1] =~ m/\,/)
					{
						($refDP,$altDP) = split(/\,/,$sampleInfo[1]);
						## calculate allelic fraction only for het variants
						my @gcounts = split(/\//,$sampleInfo[0]);
						if ($gcounts[0] eq '0' or $gcounts[1] eq '0')
						{
							if (($altDP + $refDP) > 0)
							{
								$allelefrac = $altDP / ($altDP + $refDP);
								$sample_count_for_allelefrac = $sample_count_for_allelefrac + 1;
							}
						}
					}

					## sample genptype filters
					# $ad_alt_allele >= 3;
					# $dp >= 8;
					# $gq >= 10;
					if ($altDP < $ad_alt_allele or $sampleInfo[2] < $dp or $sampleInfo[3] < $gq)
					{
						## set genotype to no call
						$data[$i] = "./.";
					}
				}else{
					## homozygous sites
					$sample_count = $sample_count + 1;
				}
			}
		}

		## loop over the sample genotypes again to get corrected genotype counts for calculating call rate
	        for (my $i = 9; $i < @data ; $i++)
        	{
                	if ($data[$i] =~ m/\:/)
                	{
				## data field format and value :: GT:AD:DP:GQ:PL => 0/0:44,0:44:99:0,117,1274
                        	my @sampleInfo = split(/\:/,$data[$i]);
                        	if ($sampleInfo[0] ne './.')
                        	{
                                	$genotype_count = $genotype_count + 1;
					$seen_avg_gq = $seen_avg_gq + $sampleInfo[3];
                        	}
                	}
                	$samples = $samples + 1;
        	}
	
		
		## calculate average allele fraction across samples for het variants
		if ($sample_count_for_allelefrac > 0)
		{
			$allelefrac = $allelefrac/$sample_count_for_allelefrac;
		}

		## calculate call rate percentage
		$seen_call_rate = ($genotype_count / $samples) * 100;
	
		## calculate average gq across sample genotypes
		$seen_avg_gq = $seen_avg_gq / $genotype_count;

        	## variant row wise filters
        	# $allelefrac between $min_avg_af and $max_avg_af for het variants in 30 sample genotypes
        	# call_rate >= 0.8
        	# $avg_gq >= 20
		if ($seen_call_rate >= $call_rate and $seen_avg_gq >= $avg_gq)
		{
			if ($sample_count >= 30 or $allelefrac ne "NA")
			{
				if ($allelefrac >= $min_avg_af and $allelefrac <= $max_avg_af)
				{
					## print lines in the new file
        				print NEW join("\t",@data)."\n";
				}
			}else{
				## print lines in the new file
				print NEW join("\t",@data)."\n";
			}
		}else{
			## print lines in the bad file
			print BAD join("\t",@data)."\n"; 
		}
	}	
}

## close handlers
close(NEW);
close(BAD);
close(VCF);
