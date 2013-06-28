#!/usr/bin/perl

#transform a normal bed file (with enrichment probes a la nimblegene)
#to a shore compatible format

#chr1    69114   69332
#chr1    69383   69577
#chr1    69644   69940
#chr1    69951   70028
#chr1    566170  566275

#chr     pos     end     type
#1       1       14466   depleted
#1       14467   14587   enriched
#1       14588   14638   depleted

use strict;
use warnings;
use Data::Dumper;

open(FH, $ARGV[0]) or die( &usage );
open(OUTFH, ">".$ARGV[0]."_shore");

#header
print OUTFH "chr\tpos\tend\ttype\n";

my $chr = 'NA';
my $chr_old = 'NA';
my $start = 1;
my $end_old = 1;

my %chromosomes = (1 => 249250621,
		   2 => 243199373,
		   3 => 198022430,
		   4 => 191154276,
		   5 => 180915260,
		   6 => 171115067,
		   7 => 159138663,
		   8 => 146364022,
		   9 => 141213431,
		   10 => 135534747,
		   11 => 135006516,
		   12 => 133851895,
		   13 => 115169878,
		   14 => 107349540,
		   15 => 102531392,
		   16 => 90354753,
		   17 => 81195210,
		   18 => 78077248,
		   19 => 59128983,
		   20 => 63025520,
		   21 => 48129895,
		   22 => 51304566,
		   23 => 155270560,
		   24 => 59373566,
		   25 => 16571
		  );


while(<FH>) {
    chomp();
    $_ =~ s/\s+/ /g;
    
    next if $_ !~ /^chr/;
    next if $_ =~ /random/; # there are some entries like
    next if $_ =~ /chrUn/;
    my @line_elements = split(/\s/);


    $line_elements[0] =~ /(\d+|Y|X)/;
    $chr = $1;
    $chr = translate($chr);

    if ($chr ne $chr_old) {
        $chr_old = $chr;
        $start = 1;
        
	if (defined($chromosomes{$chr - 1})) {
	    print OUTFH $chr - 1, "\t", $end_old, "\t$chromosomes{$chr - 1}\tdepleted\n"; #depleted
	}
        print OUTFH "$chr\t", $start , "\t" , $line_elements[1] - 1 , "\tdepleted\n"; #depleted
        print OUTFH "$chr\t$line_elements[1]\t$line_elements[2]\tenriched\n"; #enriched
        $start = $line_elements[2] + 1;
	
    }
    elsif ($chr eq $chr_old) {
        print OUTFH "$chr\t$start\t" , $line_elements[1] - 1 , "\tdepleted\n"; #depleted
        print OUTFH "$chr\t$line_elements[1]\t$line_elements[2]\tenriched\n"; #enriched
        $start = $line_elements[2] + 1;
	$end_old = $line_elements[2] + 1;
    }   
}

close FH;
close OUTFH;

sub translate {
	my $Chr = shift;
	return 23 if ($Chr =~ /x/i);
	return 24 if ($Chr =~ /y/i);
	return 25 if ($Chr =~ /mt/i);
	#return 'X' if ($Chr == 23);
	#return 'Y' if ($Chr == 24);
	return $Chr;
}

sub usage {
    print "\n\n\n","transform a normal bed file (with enrichment probes a la nimblegene)",
"to a shore compatible format",
"\n",
"chr1    69114   69332\n",
"chr1    69383   69577\n",
"chr1    69644   69940\n",
"chr1    69951   70028\n",
"chr1    566170  566275\n",
"\n",
"chr     pos     end     type\n",
"1       1       14466   depleted\n",
"1       14467   14587   enriched\n",
"1       14588   14638   depleted\n\n\n"
}
