#!/usr/bin/perl
use warnings;
use strict;
use Digest::MD5 qw(md5);
use DBI;
use threads;
use threads::shared;
#use DateTime;
use Getopt::Long;



sub usage { print "\n$0 \n usage:\n",
	   "--input \t provide the VCF file to annoate\n",
	   "--tempDir \t provide temporary location (all temp files with the current session will be removed at the end of execution) \n",
	   "--geneDef \t the gene deifnition you want to select for Annovar (ensGene,refGene,knownGene) [default: all] \n",
	   "--help \t show help \n",
	   "\n\n perl $0 --input=/path_to/inout.vcf --tempDir=/path_to_temp",
	   "\n For more details, send an email to : rubayte.rahman\@crg.eu\n";
}


##############################################################################################
## SETTINGS
##############################################################################################


## variables
my $help = 0;
my $input; ## main input file
my $geneDef = "all"; ## geneDefinition
my $sep = ",";
my %snps = (); ## hash to hold input SNPs from VCF
my @thrds = (); ## list to hold threads
our %eDiVa :shared = (); ## hash to hold input INDELs with eDiVa annotation
our %SIFT :shared = (); ## hash to hold input INDELs with SIFT annotation
our %Annovar :shared = (); ## hash to hold input INDELs with ANNOVAR annotation
our %samples :shared = (); ## hash to hold sample information
our @headers = ();

## temp location
our $templocation;##= "/users/so/rrahman/IndelAnnotationScratch";
#our $fileSuffix = DateTime->now->ymd;
#$fileSuffix = $fileSuffix."-".DateTime->now->hms;
#$fileSuffix =~ s/\:/-/g;

my $fileSuffix = localtime();
$fileSuffix =~ s/\:/-/g;
$fileSuffix =~ s/\s/-/g;

#our $fileSuffix = "2013-09-23-16-18-45";


## SIFT settings
#our $SIFTDirectory = "/users/GD/tools/sift5.0.3";
#our $siftCodingInfoGffFile = "ens.hum.ncbi37.ver66.cds.merge.gff";

## ANNOVAR settings
our $ANNOVAR = "/home/rrahman/soft/lib/Annovar"; 

## grab command line options
GetOptions("input=s" => \$input, "tempDir=s" => \$templocation, "geneDef=s" => \$geneDef, "help=s" => \$help);

## check command line parameters and take necessary actions
unless(($input and $templocation) && $help == 0)
{
	usage;
	exit 0;
	## bye bye
}

## final check of geen definition
if ($geneDef ne "refGene" && $geneDef ne "ensGene" && $geneDef ne "knownGene" && $geneDef ne "all")
{
	print "\nWARNING :: Not a valid gene definition. Please select a correct gene definition and if you are not sure then use the default settings !\n";
	usage;
	exit 0;
	## bye bye
}

##############################################################################################
## SUBROUTINES
##############################################################################################

## subroutine for finalizing annotation process
sub finalize
{
	## clear the tmp directory for this session
	my $clearCmm = "rm -r $templocation/*".$fileSuffix."*";
	system ($clearCmm);
}

## subroutine for eDiVa annotation
sub eDiVaAnnotation
{
	## DB parameters
	my $username = 'hana';
	my $database = 'eDiVa_innoDB';
	my $sql = "";
	my $stmt ;
	my @res = ();
	## open DB connection
	my $dbh = DBI->connect('dbi:mysql:'.$database.';host=localhost',$username,'hanamysql2013') or die "Connection Error!!\n";

	# extract db result
	while( my ($k, $v) = each %snps ) 
	{
		my ($chr,$pos,$ref,$alt) = split(/\;/,$k);
		$sql = "select annotateSNPGermline('$chr',$pos,'$ref','$alt');";
		$stmt = $dbh->prepare($sql);
		$stmt->execute or die "SQL Error!!\n";
	
		#process query result
		while (@res = $stmt->fetchrow_array) 
		{
			# load eDiVa hash from database
			$eDiVa{ $k } = $res[0];
		}
    }
	## close DB connection
	$dbh->disconnect();
}


## subroutine for Annovar annotation
sub AnnovarAnnotation
{	
	my $annCmm;
	## prepare Annovar input
	#my $annInCmm = "perl $ANNOVAR/convert2annovar.pl --includeinfo -format vcf4 $input > $templocation/annInfile".$fileSuffix."   2> ".$input.".annovar.log";
	#print "MESSAGE :: Running Annovar command \> $annInCmm\n";
	#system ($annInCmm);

	my $annFile = "$templocation/annInfile".$fileSuffix."";

	open (D, $input);
	open (N, ">>".$annFile);
	while (<D>)
	{
		chomp $_;
		if ($_ !~ m/^#/)
		{
			my @lines = split(/\t/,$_);

			#my @line = split(/\;/,$lines[1]);

			my $chr = $lines[0];
			my $position = $lines[1];
			my $ref = $lines[3];
			my $alt = $lines[4];
		
			if ($alt =~ m/\,/)
			{
				my @als = split(/\,/,$alt);
			
				foreach my $as (@als)
				{
					print N $chr."\t".$position."\t".$position."\t".$ref."\t".$as."\n";
				}
			}else
			{
				print N $chr."\t".$position."\t".$position."\t".$ref."\t".$alt."\n";
			}
		}	
	}

	#close(D);
	#close(N);

	## load annovar hash
	open (LOAD, $annFile) or die "Cant open $annFile file \n";
	while (<LOAD>)
	{
		chomp $_;
		my @dt = split(/\t/,$_);
		#if ($dt[9] =~ m/\,/)
		#{
		#	my @alss = split(/\,/,$dt[9]);
		#	foreach my $alsss (@alss)
		#	{
		#		my $newKey = $dt[5].";".$dt[6].";".$dt[8].";".$alsss;
				#my $newVal = $dt[0].";".$dt[1].";".$dt[2].";".$dt[3].";".$dt[4];
		#		$Annovar { $newKey } = "NA";
		#	}
		#}else{
			#my $newKey = $dt[5].";".$dt[6].";".$dt[8].";".$dt[9]; ## for VCF
			my $newKey = $dt[0].";".$dt[1].";".$dt[3].";".$dt[4];
			#my $newVal = $dt[0].";".$dt[1].";".$dt[2].";".$dt[3].";".$dt[4];
			$Annovar { $newKey } = "";
		#}
	
	}
	close(LOAD);

	## run Annovar annotation
	if ($geneDef eq 'ensGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);		
	}elsif($geneDef eq 'refGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}elsif($geneDef eq 'knownGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}elsif($geneDef eq 'all')
	{
		print "MESSAGE :: No sepicific gene definition selected, hence Annovar is going to run on all definitions !\n";
		## refgene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## ensgene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## knowngene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}else
	{
		print "MESSAGE :: Not a valid gene definition ! Exiting ..";
		&finalize;
		exit 0;
	}

	## read annotation from annovar
	my $annFianlAnnE = "$templocation/Ensembl".$fileSuffix.".genome_summary.csv";
	my $annFianlAnnR = "$templocation/Refseq".$fileSuffix.".genome_summary.csv";
	my $annFianlAnnK = "$templocation/Known".$fileSuffix.".genome_summary.csv";
	
	if (-e $annFianlAnnR)
	{
		open (ANNR, $annFianlAnnR) or die "Cant open $annFianlAnnR file \n";
		## read refseq annotation
		while (<ANNR>)
		{
			if ($_ !~ m/^Func/)
			{
				chomp $_;
				my @dt = split(/\,/,$_);
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[24].";".$dt[25];
				##print $valueTOmatch."\n";
				
				## fix missing values
				if ($dt[0] eq '')
				{
					$dt[0] = 'NA';
				}
				if ($dt[1] eq '')
				{
					$dt[1] = 'NA';
				}
				if ($dt[2] eq '')
				{
					$dt[2] = 'NA';
				}

				if ($dt[3] eq '')
				{
					$dt[3] = 'NA';
				}
				
				
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];
				#if ($geneDef eq 'refGene')
				#{
		        #	$Annovar { $valueTOmatch } = $annToPass;
				#}else{
					$Annovar { $valueTOmatch } = $annToPass;
				#}
			}
		}
		close(ANNR);
	}

	if (-e $annFianlAnnE)
	{
		open (ANNE, $annFianlAnnE) or die "Cant open $annFianlAnnE file \n";
		## read ens annotation
		while (<ANNE>)
		{
			if ($_ !~ m/^Func/)
			{
				chomp $_;
				my @dt = split(/\,/,$_);
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[24].";".$dt[25];
				#print $valueTOmatch."\n";
				
				## fix missing values
				if ($dt[0] eq '')
				{
					$dt[0] = 'NA';
				}
				if ($dt[1] eq '')
				{
					$dt[1] = 'NA';
				}
				if ($dt[2] eq '')
				{
					$dt[2] = 'NA';
				}

				if ($dt[3] eq '')
				{
					$dt[3] = 'NA';
				}				
				
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];

				if ($geneDef eq 'ensGene')
				{
		        	$Annovar { $valueTOmatch } = $annToPass;
				}else{
		        	$Annovar { $valueTOmatch } = $Annovar { $valueTOmatch }.",".$annToPass;
				}

			}
		}
	close(ANNE);
	}

	if (-e $annFianlAnnK)
	{
		open (ANNK, $annFianlAnnK) or die "Cant open $annFianlAnnK file \n";		
		## read known annotation
		while (<ANNK>)
		{
			if ($_ !~ m/^Func/)
			{
				chomp $_;
				my @dt = split(/\,/,$_);
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[24].";".$dt[25];

				
				## fix missing values
				if ($dt[0] eq '')
				{
					$dt[0] = 'NA';
				}
				if ($dt[1] eq '')
				{
					$dt[1] = 'NA';
				}
				if ($dt[2] eq '')
				{
					$dt[2] = 'NA';
				}

				if ($dt[3] eq '')
				{
					$dt[3] = 'NA';
				}
				
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];
				
				if ($geneDef eq 'knownGene')
				{
		        	$Annovar { $valueTOmatch } = $annToPass;
				}else{
		        	$Annovar { $valueTOmatch } = $Annovar { $valueTOmatch }.",".$annToPass;
				}
				        #$Annovar { $valueTOmatch } = $Annovar { $valueTOmatch }.",".$annToPass;
		        #	}
				#}
			}
		}
	close(ANNK);
	}
}

## subroutine for fetching sift score from sift sqlite3 database
sub siftScorefetch
{
	my $c = shift;
	my $pos = shift;
	my $r = shift;
	my $a = shift;
	my $trs = shift;
	my $score;
	if ($trs =~ m/^ENST/)
	{

	}elsif($trs =~ m/^NM/)
	{

	}elsif($trs =~ m/^uc/)
	{

	}else
	{
		print "\nWARNING :: Sift can't recognize ANNOVAR transcript $trs !\n";
	}
	return $score;
}

## subroutnine por providing header to the main annotation output file
sub getHeader
{
	my $stringTOreturn; ## header to return
	## check for gene definiton and construct header according to that
	if($geneDef eq 'ensGene')
	{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
		AminoAcidChange(Ensembl),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>Cov>AF)";
	}elsif($geneDef eq 'refGene')
	{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),
		dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>Cov>AF)";
	}elsif($geneDef eq 'knownGene')
	{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,
		TotalEVSFrequecy,Eur1000GenomesFrequency,Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>Cov>AF)";
	}else{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
		AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>Cov>AF)";
	}
	## replace newlines with nothing at header line
	$stringTOreturn =~ s/\n|\s+//g;
	return $stringTOreturn;
}

########################################################################################################################
## MAIN starts
########################################################################################################################

## start annotation process
print "MESSAGE :: Intializing annotation process .. ";

## open input file and load input hash
open (INPUT, $input) or die "\ncant open $input file !!\n";

while(<INPUT>)
{
	chomp $_;
	
	if ($_ =~ m/CHROM/)
	{
		@headers = split(/\t/,$_);
	}
	
	if ($_ !~ m/^#/) ## skip header line
	{
		my @line = split(/\t/,$_);
		my $chr = $line[0];
		my $position = $line[1];
		my $ref = $line[3];
		my $alt = $line[4];
		my @infos = split(/\;/,$line[7]);
		my $AF;
		
		foreach my $info (@infos)
		{
			if ($info =~ m/^AF=/)
			{
				$AF = substr($info,3);
			}
		}
		

		## take care of chr1 or Chr1 and convert to chr1 /Chr1-> 1
		if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
		{
			$chr = substr($chr,3);
		}
		## process based on alteration
		if ($alt =~ m/\,/)
		{
			my @alts = split(/\,/,$alt);
			my @afs = split(/\,/,$AF);
			
			for (my $j = 0; $j < @alts; $j++)
			#foreach my $al(@alts)
			{
				my $al = $alts[$j];
				my $alfr = $afs[$j];
				
				$snps{ "$chr;$position;$ref;$al" } = "NA";
				for(my $i = 9; $i < @line; $i++)
				{
					my @gts = split(/\:/,$line[$i]);
					if ($gts[0] ne './.')
					{
						if (exists $samples { "$chr;$position;$ref;$al" })
						{
							$samples { "$chr;$position;$ref;$al" } = $samples { "$chr;$position;$ref;$al" }.";".$headers[$i].">".$gts[0].">".$gts[2].">".$alfr; 
						}else{
							$samples { "$chr;$position;$ref;$al" } = $headers[$i].">".$gts[0].">".$gts[2].">".$alfr;
						}
					}
				}			
			}
		}else{
			$snps{ "$chr;$position;$ref;$alt" } = "NA";
				for(my $i = 9; $i < @line; $i++)
				{
					my @gts = split(/\:/,$line[$i]);
					if ($gts[0] ne './.')
					{
						if (exists $samples { "$chr;$position;$ref;$alt" })
						{
							$samples { "$chr;$position;$ref;$alt" } = $samples { "$chr;$position;$ref;$alt" }.";".$headers[$i].">".$gts[0].">".$gts[2].">".$AF; 
						}else{
							$samples { "$chr;$position;$ref;$alt" } = $headers[$i].">".$gts[0].">".$gts[2].">".$AF;
						}
					}
				}
		}

	} ## end if
} ## end while
close(INPUT);



#open (INPUT, $input) or die "\ncant open $input file !!\n";

#while(<INPUT>)
#{
#	chomp $_;
#	if ($_ !~ m/^#/) ## skip header line
#	{
#		my @lines = split(/\t/,$_);

#		my @line = split(/\;/,$lines[1]);

#		my $chr = $line[0];
#		my $position = $line[1];
#		my $ref = $line[2];
#		my $alt = $line[3];

		## take care of chr1 or Chr1 and convert to chr1 /Chr1-> 1
#		if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
#		{
#			$chr = substr($chr,3);
#		}
		## process based on altearation
#		if ($alt =~ m/\,/)
#		{
#			my @alts = split(/\,/,$alt);
#			foreach my $al(@alts)
#			{
#				$snps{ "$chr;$position;$ref;$al" } = "NA";
#			}
#		}else{
#			$snps{ "$chr;$position;$ref;$alt" } = "NA";
#		}

#	} ## end if
#} ## end while
#close(INPUT);





## Initialization completed
print "\nMESSAGE :: Intialization completed .. ";

## start threading and annotating
print "\nMESSAGE :: Annotation starting .. This may take a few moments ..\n";

## spawn threds and push to list
push @thrds, threads -> new(\&eDiVaAnnotation); ## spawn a thread for eDiVa annotation
push @thrds, threads -> new(\&AnnovarAnnotation); ## spawn a thread for Annovar annotation

## join spawned threads
foreach my $thr (@thrds)
{
        $thr -> join();
}

#&AnnovarAnnotation;
#&eDiVaAnnotation;


## annotation completed
print "\nMESSAGE :: Annotation completed .. ";

## write annotaton in file 
print "\nMESSAGE :: Writing annotation to output file .. ";


## file to write data eDiVa
my $outFile = "$input.annotated"; ## output file to write annotation with input CNV definitions
unlink($outFile); ## delete the file if exists
open (ANN, ">>".$outFile) or die "Cant open new file to write annotation !!\n";

## write header to main output file
my $headerOutputFile = &getHeader;
print ANN  $headerOutputFile."\n";
## write data lines now
while (my($key, $value) = each(%snps)) 
{
	my ($edivaannotationtoprint,$annovarannotationtoprint,$sftScr) = ("NA","NA","NA");
	my ($chr,$position,$ref,$alt) = split(/\;/, $key);
	$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
	$annovarannotationtoprint = $Annovar{$key} if $Annovar{$key};

	## clean annovar annotation
	#$annovarannotationtoprint =~ s/\"//g;
	#$annovarannotationtoprint =~ s/,,/,NA,/g;
	#$annovarannotationtoprint =~ s/NA,/NA,NA/g;
	
	## clean annovar annotation
	$annovarannotationtoprint =~ s/\"//g;
	#$annovarannotationtoprint =~ s/,,/,NA,/g;

	## grab annovar transcript for quering SIFT score 
	#my @ans = split(/\,/,$annovarannotationtoprint);
	#if ($ans[2] ne "NA")
	#{
	#	sftScr = &siftScorefetch($chr,$position,$ref,$alt,$ans[2]);
	#}
	## write annotation to file 
	print ANN $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$samples { $key }."\n";
}

close(ANN);

## sort the file
#my $srtCmm = "sort -n -k1,1 -n -k2,2 --field-separator=, $outFile > $outFile.sorted";
#system($srtCmm);

## writing completed
print "\nMESSAGE :: Writing annotation completed .. ";
#print "\nMESSAGE :: Your annotated files are $outFile \& $outFile.sorted";
print "\nMESSAGE :: Your annotated file is $outFile";

## Finalize everything
print "\nMESSAGE :: Finalizing annotation process ..";
&finalize;
print "\nMESSAGE :: Finalization completed !\n";
