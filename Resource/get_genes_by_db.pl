use warnings;
use strict;
use DBI;

### Usage
if (@ARGV != 1)
{
        print "usage:: $0 <gene_list_file> \n";
        exit 0;
}

### DB connect 
my $dbh = DBI->connect('dbi:mysql:'.'eDiVa_innoDB'.';host=www.ediva.crg.eu','rrahman','Temp1234') or die "Connection Error!!\n";


### Open input file
open (GENES, $ARGV[0]) or die "cant open file $ARGV[0] \n";


### Open output file
unlink($ARGV[0].'.db');
open (DB, ">>".$ARGV[0].'.db') or die "cant open new file \n";


### Read input file
while (<GENES>)
{
        chomp $_;
        my $gene = $_;

	### SQL query
        my $sqlQ = "SELECT geneName,name,chr,txStart,txEnd,cdsStart,cdsEnd,exonStarts,exonEnds FROM eDiVa_innoDB.Table_refGene where geneName = \'$gene\' and (chr between 1 and 22 or chr = \'X\') order by exonCount desc \;";
        my $stmt = $dbh->prepare($sqlQ);
        $stmt->execute or die "SQL Error!!\n";
        while (my @row = $stmt->fetchrow_array)
        {
                print DB join("\t", @row) . "\n";
        }
}

# Finish
$dbh->disconnect();
close(GENES);
exit(0);

