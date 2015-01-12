
if (!(defined($ARGV[0])) || !(defined($ARGV[1]))) {
    print "\n\n\n usage: perl splitGATKoutput.pl outfolder GATKvariant.vcf \n\n\n";
    exit(0);
}

my $path = $ARGV[0];
my $infile = $ARGV[1];

open(IN, $infile);
open(SNP, ">$path/GATK.snps.raw.vcf");
open(INDEL, ">$path/GATK.indel.raw.vcf");


while (<IN>) {
    my @splitline = split(/\t/, $_);
    
    if ($_ =~ /^\#/) {
        print SNP $_;
        print INDEL $_;
    }
    else {
        if (length($splitline[3]) == 1 && length($splitline[4]) == 1) {
            print SNP $_;
        }
        elsif (length($splitline[3]) != 1 || length($splitline[4]) != 1) {
            print INDEL $_;
        }
    }
}