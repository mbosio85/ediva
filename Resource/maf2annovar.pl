#!/usr/bin/perl
# 
use Pod::Usage;
use Getopt::Long;
use strict;
use warnings;

my ($help, $manual);
GetOptions( 'help|h' => \$help,
	    'manual|man|m' => \$manual,
    )or pod2usage ();
@ARGV == 1 or pod2usage ();

open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my %sample;

my $title = <IN>;
while ($title =~ /^#/){
   $title =<IN>
}

my @title = split(/\t/,$title);
my ($Chromosome, $Start_Position, $Tumor_Sample_Barcode, $Reference_Allele, $Tumor_Seq_Allele1, $Tumor_Seq_Allele2);
for(0..$#title){
    $Chromosome = $_ if $title[$_] =~ "Chr";
    $Start_Position = $_ if $title[$_] =~ "Start";
  #  $Tumor_Sample_Barcode = $_ if $title[$_] eq "Tumor_Sample_Barcode";
    $Reference_Allele = $_ if $title[$_] eq "Reference_Allele";
    $Tumor_Seq_Allele1 = $_ if $title[$_] eq "Tumor_Seq_Allele1";
    $Tumor_Seq_Allele2 = $_ if $title[$_] eq "Tumor_Seq_Allele2";
}

while(<IN>){
    chomp;
    my @line = split("\t", $_);
    next unless defined $line[$Tumor_Seq_Allele1];
    my $end = $line[$Start_Position] + length($line[$Tumor_Seq_Allele1]) - 1;
    my $sample = "all_samples_converted.annovar" ; #$line[$Tumor_Sample_Barcode];
    if($line[$Tumor_Seq_Allele1] eq $line[$Reference_Allele]){
	$sample{$sample}{"$line[$Chromosome]\t$line[$Start_Position]\t$end\t$line[$Reference_Allele]\t$line[$Tumor_Seq_Allele2]"} = 1;
    }else{
	$sample{$sample}{"$line[$Chromosome]\t$line[$Start_Position]\t$end\t$line[$Reference_Allele]\t$line[$Tumor_Seq_Allele1]"} = 1;
    }
}

close IN;
foreach my $sample(sort keys %sample){
   #  open(OUT, ">$sample") or die "cannot open file $sample:$!\n";    
    foreach my $line (sort keys %{$sample{$sample}}){
        print  "$line\n";
    }
}





######################################################################################################################################
############################################################ manual page #############################################################
######################################################################################################################################

=head1 NAME                                                                                                        
maf2annovar is a script to convert maf file to a list of annovar input files
                                                              
=head1 SYNOPSIS
                                                                                                                                                                                                      
 perl maf2annovar.pl [options] <input>
 Option: 
        -h, --help                      print help message   
        -m, --manual                    print manual message
        

 Function: this script converts maf file to a list of annovar input files (each file represents a sample)

 Example: perl maf2annovar.pl /path/to/maf_file
 
=head1 OPTIONS

=over 8

=item B<--help>

 print a brief usage message and detailed explanation of options.

=item B<--manual>

 print the manual page and exit.

=back

=cut
