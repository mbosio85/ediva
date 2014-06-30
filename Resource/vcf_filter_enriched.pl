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
#  Module: Parser::VCF::vcf_filter::vcf_filter_enriched.pl
#  Purpose:
#  In:
#  Out:
#


####################################################################
# Overlap Genotyping SNPs with NGS SNPs
# Author: Stephan Ossowski
####################################################################


my $usage = "\n\n$0 enrichment_file SNP_VCF\n\n";
my $enriched = shift or die $usage;
my $snp_vcf  = shift or die $usage;


# Get enriched regions
my %enr = ();
open ENRICHED, $enriched or die "Cannot open $enriched file\n";
while( <ENRICHED> ) {
	chomp;
	## skip header line
	if ($_ !~ m/^#/)
	{
		my @a = split("\t", $_);
		#if($a[3] eq "enriched") 
		#{
		for( my $i = $a[1]; $i <= $a[2]; $i++) 
		{
			$enr{$a[0]}{$i} = 1;
		}
		#}
	}	
}
close ENRICHED;



# Read VCF and return SNP from enriched regions
open VCF, $snp_vcf or die "Cannot open $snp_vcf file\n";
while( <VCF> ) {

	chomp;

	if($_ =~ /#/) {
		print "$_\n";
	}
	else {
		my @a = split("\t", $_);

		if( ! exists $enr{$a[0]}{$a[1]} ) {
			$a[6] .= ",NOTENRICHED";
		}

		print join("\t",@a);
		print "\n";
		#for (my $i = 0; $i < $#a; $i++) {
		#	print "$a[$i]\t";
		#}
	}
}

close VCF;
exit(0);
