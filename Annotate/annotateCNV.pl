#!/usr/bin/perl
use warnings;
use DBI;

## process command line parameters

if (@ARGV != 1)
{
	print "USAGE :: perl annotateCNV.pl input_file_with_CNVs[bed_format]\n\n\n";
	print "NOTICE :: Currently the annotation is only based on Ensembl gene definition.\n";
	exit 0;
}

## DB parameters
my $username = 'hana';
my $database = 'eDiVa_innoDB';

## open DB connection
my $dbh = DBI->connect('dbi:mysql:'.$database.';host=www.ediva.crg.eu',$username,'hanamysql2013') or die "Connection Error!!\n";

## files to write data from eDiVa
my $outFile = "$ARGV[0].annotated"; ## output file to write annotation with input CNV definitions
unlink($outFile); ## delete the file if exists
open (ANN, ">>".$outFile) or die "Cant open new file to write annotation !!\n";

## print annotation header
print ANN "chromosome\tStart\tEnd\tGene_Overlap (format => Gene:refseqID:%ofCNV:%ofGene)\tDGV_CNV_Overlap (format => CNVlandmark:%ofCNV:%ofDGV_CNV:DGV_CNV_AF)\tAll_Genes\tProbable_same_DGV_CNV\tRest_of_Input_file_columns\n";


## open input file and process
open (INPUT, $ARGV[0]) or die "cant open $ARGV[0] file !!\n";

while(<INPUT>)
{
	chomp $_;
	if ($_ !~ m/^#/) ## skip header line
	{
		my @line = split(/\t/,$_);
		
		my $endind = @line;
		
		my $chr = $line[0];
		my $start = $line[1];
		my $end = $line[2];
		my @anns = ();
		my $ann = "";
		my $ProbSameCNV = "";
		my @genes = ();
		my @dgvs = ();
		my @Gs = ();
		my @inputrest = @line[3..$endind];
		
		## take care of chr1 or Chr1 and convert to chr1 /Chr1-> 1
		if ($chr =~ m/^c/ or $chr =~ m/^C/)
		{
			$chr = substr($chr,3);
		}
		
		## input check
		#if ($chr !~ m/[1-22,X,Y]/ or $start =~ m/[a-zA-Z]/ or $end =~ m/[a-zA-Z]/)
		#{
		#	print "Malformed input lines. Check Your input file !!\n";
		#	exit 0;
		#}
		
		## start processing
		my $sql = "select annotateCNV($chr,$start,$end);";
		my $stmt = $dbh->prepare($sql);
		$stmt->execute or die "SQL Error!!\n";
		
		##process query result
		while (@row = $stmt->fetchrow_array) 
		{
			#my $ann = $row[0];
			@anns = split(/\t/, $row[0]);
		}

		my $i = 0;
		foreach $a (@anns)
		{
			my $prev = "";
			if ($a eq "NA,NA,NA,NA")
			{
				$a = "NA";
				#$ann = $ann.$a."\t";
				if ($i == 1)
				{
					$ProbSameCNV = "NA";
				}
			}else{
				my @elements = split(/\,/,$a);
				foreach $item(@elements)
				{
					if ($item ne "NA" and $i == 0)
					{
						push (@genes, $item);
						
						my @tempG = split(/\:/, $item);
						push (@Gs, $tempG[0]);
						
					}elsif ($item ne "NA" and $i == 1)
					{
						push (@dgvs,$item);
					}
					else{}
					
					if ($item ne "NA" and $item ne $prev)
					{
						if ($i == 1)
						{
							#my @temp = split(/\(/,$item);
							my ($landmarkchr,$landmarkrange,$val1,$val2,$freq) = split(/\:/,$item);
							#$val2 = substr($val2,0,length($val2)-1);
							if ($val1+$val2 >= 120)
							{
								$ProbSameCNV = $ProbSameCNV.",".$landmarkchr.':'.$landmarkrange.",";
							}
						}
						#$ann = $ann.$item.";";
						$prev = $item;
					}
				}
				#$ann = substr($ann,0,length($ann)-1);
				#$ann = $ann."\t";
			}
			$i++;
		}
		
		## clear probably same CNV term
		if ($ProbSameCNV eq "")
		{
			$ProbSameCNV = "NA";
		}
		if($ProbSameCNV =~ m/,$/)
		{
			$ProbSameCNV = substr($ProbSameCNV,0,length($ProbSameCNV)-1);
		}
		if ($ProbSameCNV =~ m/^\,/)
		{
			$ProbSameCNV = substr($ProbSameCNV,1);
		}		

		$ProbSameCNV =~ s/,,/,/g;

		## clean out duplicate gene names from the gene list
        my %seen;
    	my @Gs_clean = ();

		for my $gene (@Gs) 
		{
			if ($seen{$gene}) 
			{
				next;
			}
			else {
				$seen{$gene} = 1;
				push (@Gs_clean, $gene);
			}
		}


		## write annotation to file
		##my $ann = join("\t",@anns); 
		print ANN $chr."\t".$start."\t".$end."\t".join(',',@genes)."\t".join(',',@dgvs)."\t".join(',',@Gs_clean)."\t".$ProbSameCNV."\t".join(',',@inputrest)."\n";
	
	}else{
	} ## end if
} ## end while

## close DB connection
$dbh->disconnect();

## close file connections
close(ANN);
close(INPUT);		

print "Annotation Complete from eDiva !!\n";