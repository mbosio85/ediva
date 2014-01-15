#! /usr/bin/perl
use strict;
use warnings;

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
#  Module: Parser::VCF::vcf_filter::vcf_filter.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "\n$0 vcf_file max_depth\n\n";

my $vcf       = shift or die $usage;
my $max_depth = shift or die $usage;

open VCF, $vcf or die "Cannot open input file\n";

while( my $line = <VCF> ) {
	chomp($line);

	if($line =~ /^#/) {
		print "$line\n";
	}
	else {
		my @a = split("\t", $line);
		my @b = split(";", $a[7]);	# Annotations
		my @c = split(":", $a[8]);	# Genotypes order
		my @d = split(":", $a[9]);	# Genotypes values

		# Get annotation
		my %anno = ();
		foreach my $anno_string (@b) {
			my($type, $value) = split("=", $anno_string);
			$anno{$type} = $value;
		}

		### Get genotype
		my %geno = ();
		for(my $i = 0; $i < $#c; $i++) {
			$geno{$c[$i]} = $d[$i];
		}

		my ($ref_support, $snp_support) = split(",", $geno{AD});

		### Default: $snp_support >= 3 && $anno{DP} >=5 && $anno{DP} <= 250 && $anno{AF} >= 0.3 && $anno{RE} <= 1.3 && $a[5] >= 20
		if( $snp_support >= 3 && $anno{DP} >=5 && $anno{DP} <= $max_depth && $anno{AF} >= 0.3 && $anno{RE} <= 1.3 && $a[5] >= 20 ) { 
			print "$line\n";
		}
	}
}

close VCF;

exit(0);
