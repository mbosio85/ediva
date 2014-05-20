use warnings;
use strict;

## check command line 
if (@ARGV != 3)
{
	print "usage:: $0 <gene_file> <bp_to_extend_on_both_sides_after_merging> <output_fasta_file_name> \n";
	exit 0;
}

## vars
my $extendbp = $ARGV[1];
my %gene = ();

## browse over the generated res file from db
open (DBF, $ARGV[0].'.db') or die "cant open file \n";

while (<DBF>)
{
	chomp $_;
	my @dl = split(/\t/,$_);
	#print SAMFILE $dl[0]."\t".$dl[1]."\t".$dl[2]."\t".($dl[5] - $extendbp)."\t".($dl[6] + $extendbp)."\n";
	if (exists $gene { $dl[0] })
	{
		$gene { $dl[0] } = $gene { $dl[0] }.";".$dl[1]."\t".$dl[2]."\t".$dl[7]."\t".$dl[8];
	}else{
		$gene { $dl[0] } = $dl[1]."\t".$dl[2]."\t".$dl[7]."\t".$dl[8];
	}
}
close(DBF);

## create samtools infile
unlink($ARGV[0].'.db.samfile');
open (SAMFILE, '>>'.$ARGV[0].'.db.samfile') or die "cant open new file \n";

while ( my ($key, $value) = each(%gene) ) 
{
	if ($value =~ m/\;/)
	{
		my $minstart = "NA";
		my $maxend = "NA";
		my $chr = "NA";
		my $exon = "exon";
		#my $trs = '';
		my @dts = split(/\;/,$value);		

		## matrices of exon starts and ends for each isoform
		my @exonstarts = ();
		my @exonends = ();
		my $maxeCount = 0;
		
		for ( my $i = 0; $i < scalar @dts; $i++)
                {
                        my @line = split(/\t/,$dts[$i]);
                        my @starts = split(/\,/, $line[2]);
                        if ($i == 0)
			{
				$maxeCount = scalar @starts;
			}else{
				if (scalar @starts > $maxeCount)
				{
					$maxeCount = scalar @starts;
				}
			}
		}
                
		for ( my $i = 0; $i < scalar @dts; $i++)
                {
			for (my $j = 0; $j < $maxeCount; $j++)
			{
				$exonstarts[ $i ] [ $j ] = "NA";
				$exonends[ $i ] [ $j ] = "NA";
			}
		}

		for ( my $i = 0; $i < scalar @dts; $i++)
		{
			my @line = split(/\t/,$dts[$i]);
			$chr = $line[1];
			my @starts = split(/\,/, $line[2]);
			my @ends = split(/\,/,$line[3]);
		
			for(my $j = 0; $j < scalar @starts; $j++)
			{
				if($starts[$j] eq '')
				{
					$starts[$j] = "NA";
					$ends[$j] = "NA";
				}
				$exonstarts[ $i ] [ $j ] = $starts[$j];
				$exonends[ $i ] [ $j ] = $ends[$j];	
			}
		}	
		
		## we neeed column wise traversal
		for(my $j = 0; $j < $maxeCount; $j++)
		{
			for (my $i = 0; $i < scalar @dts; $i++)
			{
				if( $i == 0 )
				{
	                                $minstart = $exonstarts[ $i ] [ $j ];
                                        $maxend = $exonends[ $i ] [ $j ];
				}else{
					if ($exonstarts[ $i ] [ $j ] ne "NA")
					{
						if ($minstart >= $exonstarts[ $i ] [ $j ])
						{
							$minstart = $exonstarts[ $i ] [ $j ];
						}
						if ($maxend <= $exonends[ $i ] [ $j ])
						{
							$maxend = $exonends[ $i ] [ $j ];
						}
					}
				}
			}
			print SAMFILE $key."\t"."exon".($j+1)."\t".$chr."\t".($minstart - $extendbp)."\t".($maxend + $extendbp)."\n";
		}	
	

		#for ( my $i = 0; $i < scalar @dts; $i++)
		#{
		#	my @dtline = split(/\t/,$dts[$i]);
		#	$chr = $dtline[1];
		#	$trs = $trs.",".$dtline[0];
		#	if ($i == 0)
		#	{
		#		$minstart = $dtline[2];
		#		$maxend = $dtline[3];
		#	}else{
		#		if ($minstart >= $dtline[2])
		#		{
		#			$minstart = $dtline[2];
		#		}
		#		if ($maxend <= $dtline[3])
                #               {
                #                      $maxend = $dtline[3];
                #                }
		#	}
		#}
		#$trs = substr($trs,1);
		#print SAMFILE $key."\t".$trs."\t".$chr."\t".$minstart."\t".$maxend."\n";
	}else{
		my @dtline = split(/\t/,$value);
		my @es = split(/\,/,$dtline[2]);
		my @ee = split(/\,/,$dtline[3]);
		my $exon = "exon";	
	
		for(my $i = 0; $i < scalar @es; $i++)
		{
			if ($es[$i] ne '')
			{
				print SAMFILE $key."\t".$exon.($i+1)."\t".$dtline[1]."\t".($es[$i] - $extendbp)."\t".($ee[$i] + $extendbp)."\n";
			}
		}
	}
	
}

close(SAMFILE);

## get refseq from human ref genome
unlink($ARGV[2]);
open (FASTA,">>".$ARGV[2]) or die "cant open new file $ARGV[2] \n";
open(POSITIONS,$ARGV[0].'.db.samfile') or die "cant open file \n";

while(<POSITIONS>)
{
        chomp $_;
        my ($gene,$seqName,$chr,$begin,$end) = split(/\t/, $_);
        open(SAMTOOLS,"/usr/bin/samtools faidx /users/GD/resource/human/hg19/bwa7/hg19.fasta $chr:$begin-$end |");
        my $i = 0;
        while(my $line = <SAMTOOLS>)
        {
                if ($i == 0)
                {
                        my $tline = ">".$gene.",".$seqName.",";
                        $line = substr($line,1);
                        $tline = $tline.$line;
                        print FASTA $tline;
                }else{
                        print FASTA $line;
                }
                $i = 1;
        }
        close(SAMTOOLS);
}

close(POSITIONS);
close(FASTA);
