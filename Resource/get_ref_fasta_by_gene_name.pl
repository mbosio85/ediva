use warnings;
use strict;

## check command line 
if (@ARGV != 4)
{
	print "usage:: $0 <gene_file> <bp_to_extend_on_both_sides_after_merging> <output_fasta_file_name> <temp_location>\n";
	exit 0;
}

## vars
my $extendbp = $ARGV[1];
my %gene = ();
my $temploc = $ARGV[3];

## browse over the generated res file from db
open (DBF, $ARGV[0].'.db') or die "cant open file \n";

while (<DBF>)
{
	chomp $_;
	my @dl = split(/\t/,$_);
	if (exists $gene { $dl[0] })
	{
		$gene { $dl[0] } = $gene { $dl[0] }.";".$dl[1]."\t".$dl[2]."\t".$dl[7]."\t".$dl[8]."\t".$dl[5]."\t".$dl[6];
	}else{
		$gene { $dl[0] } = $dl[1]."\t".$dl[2]."\t".$dl[7]."\t".$dl[8]."\t".$dl[5]."\t".$dl[6];
	}
}
close(DBF);

## create infile for samtools
unlink($ARGV[0].'.db.samfile');
open (INFILE, '>>'.$ARGV[0].'.db.samfile') or die "cant open new file \n";

## loop over all the genes and their isoforms and exons and prepare final input file for samtools to extract reference genome in fasta format
while ( my ($key, $value) = each(%gene) ) 
{
	## process genes containing multiple transcript isoforms
	if ($value =~ m/\;/)
	{
		my $exon = "exon";
		my @dts = split(/\;/,$value);		
		
		open (GN, ">>".$temploc."/".$key);
		
		for ( my $i = 0; $i < scalar @dts; $i++)
                {
                	my @line = split(/\t/,$dts[$i]);
                       	my @starts = split(/\,/, $line[2]);
                       	my @ends = split(/\,/, $line[3]);
			my @newstarts = ();
			my @newends = ();
			my $cdsStart = $line[4];
			my $cdsEnd = $line[5];
			my $j = 0;
			my $keep = 0;
			## fix the cdsstart and cdsend in the exons and prepare new starts and ends
			for (my $i = 0; $i < scalar @starts; $i++)
			{
				if ($starts[$i] ne '')
				{
					my $nstart = $starts[$i];
					my $nend = $ends[$i];
					## check for cdsstart is inside the exon
					if ($cdsStart >= $starts[$i] and $cdsStart < $ends[$i] and $keep == 0)
					{
						$nstart = $cdsStart;
						$keep = 1;
					}
					if($cdsEnd > $starts[$i] and $cdsEnd <= $ends[$i] and $keep == 1)
					{
						$nend = $cdsEnd;
						$newstarts[$j] = $nstart;
                                                $newends[$j] = $nend;
						$keep = 0;
						## this would be the last exon to consider as this exon contains the cdsEnd
					}
					if($keep == 1)
					{
						$newstarts[$j] = $nstart;
						$newends[$j] = $nend;
						$j = $j + 1;
					}	
				}
			}
			## print the new exon starts and ends for this gene
			if(scalar @newstarts == scalar @newends)
			{
				for (my $i = 0; $i < scalar @newstarts; $i++)
				{
					if ($newstarts[$i] ne '')
					{
						print GN $line[1]."\t".$newstarts[$i]."\t".$newends[$i]."\n";
					}	
				}
			}else{
				print "debug:: ".join("\t",$dts[$i])."\n";
			}
		}
                
		close(GN);
		
		## sort the gene exon start end file and run bedtools to merge(union) over the overlapping exons
		my $genefile = $temploc."/".$key;
		my $srtCmd = "sort -n -k2,2 ".$genefile." > ".$genefile.".sorted.bed";
		system($srtCmd);
		my $bedmerCmd = "/usr/bin/bedtools merge -i ".$genefile.".sorted.bed"." > ".$genefile.".bed";
		system($bedmerCmd);
		sleep 5;
		my $gfile = $genefile.".bed";
		## read the bedtools output and print on the infile
		if (-e $gfile)
		{
			open (GN2,$gfile);
			my $i = 1;
			while (<GN2>)
			{
				chomp $_;
				my ($chr,$start,$end) = split(/\t/,$_);
				print INFILE $key."\t".$exon.$i."\t".$chr."\t".($start - $extendbp)."\t".($end + $extendbp)."\n";
				$i = $i + 1;
			}
			close(GN2);
		}
		
		## clear temploc for this gene, delete the temp files created above
		my $clCmd = "rm ".$temploc."/".$key."*";
		#system($clCmd);
		
	}else{
		my @dtline = split(/\t/,$value);
		my @es = split(/\,/,$dtline[2]);
		my @ee = split(/\,/,$dtline[3]);
		my $exon = "exon";	
		my $cdsStart = $dtline[4];
                my $cdsEnd = $dtline[5];
		my $j = 0;
		my $keep = 0;
		my @newstarts = ();
		my @newends = ();

                for (my $i = 0; $i < scalar @es; $i++)
                {
	                if ($es[$i] ne '')
                        {
         		        my $nstart = $es[$i];
                        	my $nend = $ee[$i];
                                if ($cdsStart >= $es[$i] and $cdsStart < $ee[$i] and $keep == 0)
                                {
  	                              $nstart = $cdsStart;
                                      $keep = 1;
                                }
                                if($cdsEnd > $es[$i] and $cdsEnd <= $ee[$i] and $keep == 1)
                                {
                                	$nend = $cdsEnd;
                                        $newstarts[$j] = $nstart;
                                        $newends[$j] = $nend;
                                        $keep = 0;
                                }
                                if($keep == 1)
                                {
                                        $newstarts[$j] = $nstart;
                                        $newends[$j] = $nend;
                                        $j = $j + 1;
                        	}
                	}
                }

                if(scalar @newstarts == scalar @newends)
                {
			for(my $i = 0; $i < scalar @newstarts; $i++)
			{
				if ($newstarts[$i] ne '')
				{
					print INFILE $key."\t".$exon.($i+1)."\t".$dtline[1]."\t".($newstarts[$i] - $extendbp)."\t".($newends[$i] + $extendbp)."\n";
				}
			}
		}else{
			print "debug:: ".$value."\n";
		}
	}
	
}

close(INFILE);

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

