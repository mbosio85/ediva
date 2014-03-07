#!/usr/bin/perl
use warnings;
use strict;
use Digest::MD5 qw(md5);
use DBI;
use threads;
use threads::shared;
use Getopt::Long;


######################################
#
#	Task => Annotate variants 
#	infile => VCF file with complete sample wise genotype information
#	outfile => text file (csv) with complete annotation with sample wise genotype information
#	Extra outfile => text file (csv) without genic annotation with sample wise information for variants that are not bi-allelic (e.g tri-allelic) 
#
#######################################


## subroutine for usage of the tool
sub usage { print "\n$0 \n usage:\n",
	   "--input \t\t VCF file containing the variants to annoate \n",
	   "--tempDir \t\t Temporary scratch location for temp files (all temp files with the current session will be removed at the end of execution) [default: input VCF location] \n",
	   "--geneDef \t\t Gene deifnition you want to select for genic annotation (ensGene,refGene,knownGene,all) [default: refGene] \n",
	   "--variantType \t\t Type of variants to annotate from input VCF file (SNP,INDEL,all) [default: all] \n",
	   "--sampleGenotypeMode \t complete: reports all genotypes from the input VCF file\n\t\t\t compact: reports only heterozygous and homozygous alteration genotypes from the input VCF file [default: compact] \n",
	   "--forceNewFileCreate \t If set, then it will over-write existing output annotation file with the same name \n",
	   "--help \t show help \n\n"
}


##############################################################################################
## SETTINGS
##############################################################################################


## variables
my $help = 0;
my $input; ## main input vcf file
my $geneDef = "refGene"; ## gene Definition
my $sep = ","; ## separator for annotation outfile; currently comma (,) is default;
my $type = "all"; ## type of variants to annotate from input vcf file
my $gtMode = "compact"; ## type of variants to annotate from input vcf file
my $forceDel; ## varibale for force deleting the output annotation file (if exists)
our $templocation = "INPATH"; ## scratch place for creating the temp files while annotating

my %variants = (); ## hash to hold input variants from VCF
my %not_biallelic_variants = (); ## hash to hold input variants from VCF where the site is not bi-allelic and annovar annotation is absent
my @thrds = (); ## list to hold threads

our %eDiVa :shared = (); ## hash to hold input INDELs with eDiVa annotation
our %Annovar :shared = (); ## hash to hold input INDELs with ANNOVAR annotation
our %samples :shared = (); ## hash to hold sample information
our @headers = ();

my $fileSuffix = localtime();
$fileSuffix =~ s/\:/-/g;
$fileSuffix =~ s/\s/-/g;

## ANNOVAR settings
our $ANNOVAR = "/users/GD/tools/eDiVaCommandLine/lib/Annovar"; 

##############################################################################################
## COMMAND LINE OPTIONS
##############################################################################################


## grab command line options
GetOptions("input=s" => \$input, "tempDir=s" => \$templocation, "geneDef=s" => \$geneDef, "variantType=s" => \$type, "forceNewFileCreate" => \$forceDel,"sampleGenotypeMode=s" => \$gtMode, "help=s" => \$help);

## check mandatory command line parameters and take necessary actions
unless(($input) && $help == 0)
{
	usage;
	exit 0;
}

## final check of geen definition
if ($geneDef ne "refGene" && $geneDef ne "ensGene" && $geneDef ne "knownGene" && $geneDef ne "all")
{
	print "\nWARNING :: Not a valid gene definition. Please select a correct gene definition and if you are not sure then use the default settings !\n";
	usage;
	exit 0;
}


## final check of variant type
if ($type ne "SNP" && $type ne "INDEL" && $type ne "CNV" && $type ne "all")
{
	print "\nWARNING :: Not a valid variant type. Please select a correct variant type and if you are not sure then use the default settings !\n";
	usage;
	exit 0;
}

## final check of sample genotype mode type
if ($gtMode ne "complete" && $gtMode ne "compact")
{
	print "\nWARNING :: Not a valid sample genptype mode value. Please select a correct value and if you are not sure then use the default settings !\n";
	usage;
	exit 0;
}


## final check of temp location
if ($templocation ne "INPATH" and !(-d $templocation))
{
	print "\nWARNING :: Not a valid temporary location or the program does not have access to this location. If you are not sure then use the default settings !\n";
	usage;
	exit 0;
}



##############################################################################################
## CONFIGURATION VARIFY
##############################################################################################

## check of annovar settings
my $vcftoAnn = "$ANNOVAR/convert2annovar.pl";
my $genicAnn = "$ANNOVAR/ediva_summarize_annovar.pl";


if (!(-d $ANNOVAR))
{
	print "\nERROR :: Annovar library location is not found. Please set the Annovar library path correctly inside the program \n";
	exit 0;
}

if (!(-e $vcftoAnn))
{
	print "\nERROR :: Program for converting VCF to Annovar format is not found in the Annovar library location \n";
	exit 0;
}

if (!(-e $genicAnn))
{
	print "\nERROR :: Program for genic annotation is not found in the Annovar library location \n";
	exit 0;
}

## check DB connection
my $dbConn = DBI->connect('dbi:mysql:'.'eDiVa_innoDB'.';host=www.ediva.crg.eu','hana','hanamysql2013') or die "\nERROR :: Could not connect to eDiVa Database Server \n";


##############################################################################################
## OUTPUT FILE(s)
##############################################################################################

my ($outFile,$SortedoutFile,$outFileIns); ## files to write annotation

if ($input =~ m/\//) ## input file with full path
{
	my @paths = split(/\//,$input);
	my $len = scalar @paths;
	my $mainOutpath = join("/",@paths[0..($len-2)]);

	my @files = split(/\./,$paths[$len-1]);
	$len = scalar @files;
	my $tempfile = join(".",@files[0..($len-2)]);
	$outFile = $mainOutpath."/".$tempfile.".annotated";
	$SortedoutFile = $mainOutpath."/".$tempfile.".sorted.annotated";
	$outFileIns = $mainOutpath."/".$tempfile.".inconsistent.annotated";
	
	## set temp location
	if ($templocation eq "INPATH")
	{
		$templocation = $mainOutpath;
	}

}else{ ## just filename

	my @files = split(/\./,$input);
	my $len = scalar @files;
	my $tempfile = join(".",@files[0..($len-2)]);
	$outFile = $tempfile.".annotated";
	$SortedoutFile = $tempfile.".sorted.annotated";
	$outFileIns = $tempfile.".inconsistent.annotated";

	## set temp location
	if ($templocation eq "INPATH")
	{
		$templocation = ".";
	}

}


## check on output file existence
if (-e $outFile or -e $SortedoutFile)
{
	## check for new file creation flag
	if ($forceDel == 1)
	{
		## delete the files if exist
		unlink($outFile);
		unlink($SortedoutFile); 
		print "MASSAGE :: Target output file(s) already exists. Removing them now \n";
	}
	else{
		print "\nWARNING :: Target output file(s) already exists. Either rename them, remove them or set the --forceNewFileCreate variable \n";
		exit 0;
	}
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
	my $dbh = DBI->connect('dbi:mysql:'.$database.';host=www.ediva.crg.eu',$username,'hanamysql2013') or die "Connection Error!!\n";

	# extract db result
	while( my ($k, $v) = each %variants ) 
	{
		my ($chr,$pos,$ref,$alt) = split(/\;/,$k);

		## decide for variant type
		my $lenref = length $ref;
		my $lenalt = length $alt;
			
		if (($lenref + $lenalt) > 2)## INDEL
		{
			$sql = "select annotateINDEL('$k','$chr',$pos);";
		}else{ ##SNP
			$sql = "select annotateSNPGermline('$chr',$pos,'$ref','$alt');";		
		}
		
		## prepare statement and query
		$stmt = $dbh->prepare($sql);
		$stmt->execute or die "SQL Error!!\n";
	
		#process query result
		while (@res = $stmt->fetchrow_array) 
		{
			# load eDiVa hash from database
			$eDiVa{ $k } = $res[0];
		}
    }

	# extract db result
	while( my ($k, $v) = each %not_biallelic_variants ) 
	{
		my ($chr,$pos,$ref,$alt) = split(/\;/,$k);

		## decide for variant type
		my $lenref = length $ref;
		my $lenalt = length $alt;
			
		if (($lenref + $lenalt) > 2)## INDEL
		{
			$sql = "select annotateINDEL('$k','$chr',$pos);";
		}else{ ##SNP
			$sql = "select annotateSNPGermline('$chr',$pos,'$ref','$alt');";		
		}
		
		## prepare statement and query
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
	my $annInCmm = "perl $ANNOVAR/convert2annovar.pl --includeinfo -format vcf4 $input > $templocation/annInfile".$fileSuffix."   2> ".$input.".annovar.log";
	##print "MESSAGE :: Running Annovar command \> $annInCmm\n";
	system ($annInCmm);

	my $annFile = "$templocation/annInfile".$fileSuffix."";


	## run Annovar annotation
	if ($geneDef eq 'ensGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);		
	}elsif($geneDef eq 'refGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}elsif($geneDef eq 'knownGene')
	{
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}elsif($geneDef eq 'all')
	{
		##print "MESSAGE :: No sepicific gene definition selected, hence Annovar is going to run on all definitions !\n";
		## refgene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## ensgene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## knowngene
		$annCmm = "$ANNOVAR/ediva_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		##print "MESSAGE :: Running Annovar command \> $annCmm\n";
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
				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
				
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\,/)
				{
					my @annalts = split(/\,/,$dt[30]);
					$dt[30] = $annalts[0];
				}
				$valueTOmatch =~ s/\"//g;
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
				$annToPass =~ s/\"//g;
				$Annovar { $valueTOmatch } = $annToPass;
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
				
				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
				
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\,/)
				{
					my @annalts = split(/\,/,$dt[30]);
					$dt[30] = $annalts[0];
				}
				$valueTOmatch =~ s/\"//g;

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
				$annToPass =~ s/\"//g;
				
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

				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
				
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\,/)
				{
					my @annalts = split(/\,/,$dt[30]);
					$dt[30] = $annalts[0];
				}
				$valueTOmatch =~ s/\"//g;
				
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
				$annToPass =~ s/\"//g;
				
				if ($geneDef eq 'knownGene')
				{
		        	$Annovar { $valueTOmatch } = $annToPass;
				}else{
		        	$Annovar { $valueTOmatch } = $Annovar { $valueTOmatch }.",".$annToPass;
				}
			}
		}
	close(ANNK);
	}
}


## subroutnine por providing header to the main annotation output file
sub getHeader
{
    my $stringTOreturn; ## header to return
        
    ## check for gene definiton and construct header according to that
    if($geneDef eq 'ensGene')
    {
        $stringTOreturn = "Chr,Position,Reference,Alteration,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
        AminoAcidChange(Ensembl),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
        Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>DPRef>DPAlt>AF)";
    }elsif($geneDef eq 'refGene')
    {
        $stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),
        dbsnpIdentifier,dbSNPfrequency,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
        Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>DPRef>DPAlt>AF)";
    }elsif($geneDef eq 'knownGene')
    {
        $stringTOreturn = "Chr,Position,Reference,Alteration,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequency,AfrEVSFrequency,
        TotalEVSFrequency,Eur1000GenomesFrequency,Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>DPRef>DPAlt>AF)";
    }else{
        $stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
        AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
        Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>DPRef>DPAlt>AF)";
    }
    ## replace newlines with nothing at header line
    $stringTOreturn =~ s/\n|\s+//g;
    return $stringTOreturn;
}



## subroutnine por providing header to the inconsistent annotation output file
sub getHeaderIns
{
    my $stringTOreturn; ## header to return

    $stringTOreturn = "Chr,Position,Reference,Alteration,GenicAnnotation,dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
    Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
    PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,samples(sampleid>zygosity>DPRef>DPAlt>AF)";

    ## replace newlines with nothing at header line
    $stringTOreturn =~ s/\n|\s+//g;
    return $stringTOreturn;
}



########################################################################################################################
## MAIN starts
########################################################################################################################

## start processing input VCF file
print "MESSAGE :: Processing input VCF file - $input \n";

## open input file and load input hash
open (INPUT, $input) or die "\n cant open $input file \n";

## start browsing over the VCF file
while(<INPUT>)
{
	chomp $_;

	if ($_ =~ m/^##/) ## skip initial header lines
	{
		## nothing to do 
	}
	elsif ($_ =~ m/^#CHROM/) ## grab sample names with other columns from the header line
	{
		@headers = split(/\t/,$_);
	}
	else ## data lines
	{
		my @line = split(/\t/,$_);
		my $chr = $line[0];
		my $position = $line[1];
		my $ref = $line[3];
		my $alt = $line[4];
		my @infos = split(/\;/,$line[7]);
		my $AF;
		
		## always test for complete genotype format field consistency in the VCF; if abnormal report for that variant
		if ($line[8] =~ m/\:/)
		{
			my @gtcheck = split(/\:/,$line[8]);
			if (@gtcheck < 5 or @gtcheck >= 7)
			{
				print "WARNING:: weird genotype format found in $input at => $chr and $position \n";
			}
		}else{
			print "WARNING:: weird genotype format found in $input at => $chr and $position \n";
		}
			
		## take care of chr1 or Chr1 and convert to chr1 /Chr1-> 1
		if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
		{
			$chr = substr($chr,3);
		}

		## grab the AF from the INFO field in the VCF
		foreach my $info (@infos)
		{
			if ($info =~ m/^AF=/)
			{
				$AF = substr($info,3);
				last;
			}
		}

		## process based on alteration
		if ($alt =~ m/\,/) ## section for all sites other than bi-allelic
		{
			my @alts = split(/\,/,$alt);
			my @afs = split(/\,/,$AF);
			
			## process each alteration allele
			for (my $j = 0; $j < @alts; $j++)
			{
				my $al = $alts[$j];
				my $alfr = $afs[$j];
				
				## decide for variant type
				my $lenref = length $ref;
				my $lenalt = length $al;
				
				if (($lenref + $lenalt) > 2)## INDEL
				{
					## make indelID
					my $token_ref = unpack('L', md5($ref));
					my $token_obs = unpack('L', md5($al));
                	
                	if ($type eq 'INDEL' or $type eq 'all')
                	{
        	        	## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
    	            	## e.g A,C in the alternate column in VCF will report only A in the main annotation file
	                	## we are doing this because we want to keep the annotation main file consistent

	                	if ($j == 0)
    	            	{
							$variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$al";
						}else{
							$not_biallelic_variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$al";
						}

						for(my $i = 9; $i < @line; $i++)
						{
			    			my ($dpref,$dpalt,$samAf);
							my @gts = ();
							
							if ($line[$i] =~ m/\:/)
							{
								@gts = split(/\:/,$line[$i]);
	                		    my @ads = split(/\,/,$gts[1]);
    	            	 		$dpref = $ads[0];
        	       		 		@ads = @ads[1..(scalar @ads -1 )];
            	    		    $dpalt = $ads[$j];
            				
            					## for missing genotype or homozygous reference genotype set the AF to 0
            					if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            					{
            						$samAf = "0";
        	    				}else{
    	            				$samAf = $alfr;
								}

    						}else{
				    			$gts[0] = $line[$i];
								$dpref = ".";
								$dpalt = "."; 

            					## for missing genotype or homozygous reference genotype set the AF to 0
            					if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            					{
            						$samAf = "0";
        	    				}else{
    	            				$samAf = $alfr;
								}
    						
    						}        				
						
							## check for sample genotype mode
							if ($gtMode eq "complete")
							{					
               			 		if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
                				{
	            	    			$samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
    							}else
        		        		{    
        	   	    	   			$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
   			         			}
   			         		}else ## compact
   			         		{
   			         			## kick out genotypes of '0/0','./.','0|0' and '.|.'
    		            		if ($gts[0] ne '0/0' and $gts[0] ne './.' and $gts[0] ne '.|.' and $gts[0] ne '0|0')
	    	            		{
	    	            			if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
                					{
	            	    				$samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
    								}else
        		        			{    
        	   	    	   				$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
   			         				}
								}
   			         		}	
						}	

					}

				}else{ ## SNP
					
                	## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
                	## e.g A,C in the alternate column in VCF will report only A in the main annotation file
                	## we are doing this because we want to keep the annotation main file consistent
                	if ($type eq 'SNP' or $type eq 'all')
                	{
	                	if ($j == 0)
    	            	{
							$variants{ "$chr;$position;$ref;$al" } = "$chr;$position;$ref;$al";
						}else{
							$not_biallelic_variants{ "$chr;$position;$ref;$al" } = "$chr;$position;$ref;$al";
						}
					
						for(my $i = 9; $i < @line; $i++)
						{
							my ($dpref,$dpalt,$samAf);
							my @gts = ();
							
							if ($line[$i] =~ m/\:/)
							{
								@gts = split(/\:/,$line[$i]);
    	            		    my @ads = split(/\,/,$gts[1]);
        	        	 		$dpref = $ads[0];
            	   		 		@ads = @ads[1..(scalar @ads -1 )];
                			    $dpalt = $ads[$j];
							    
							    ## for missing genotype or homozygous reference genotype set the AF to 0
            					if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            					{
            						$samAf = "0";
        	    				}else{
			        	        	$samAf = $alfr;
								}
							}else{
				    			$gts[0] = $line[$i];
								$dpref = ".";
								$dpalt = "."; 
							    
							    ## for missing genotype or homozygous reference genotype set the AF to 0
            					if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            					{
            						$samAf = "0";
        		    			}else{
				        	        $samAf = $alfr;
								}
							}							
                	    
							
							## check for sample genotype mode
							if ($gtMode eq "complete")
							{								
								if (exists $samples { "$chr;$position;$ref;$al" })
								{
									$samples { "$chr;$position;$ref;$al" } = $samples { "$chr;$position;$ref;$al" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf; 
								}else{
									$samples { "$chr;$position;$ref;$al" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
								}
							}else ## compact
							{
		   	            		## kick out genotypes of '0/0','./.','0|0' and '.|.'
    		            		if ($gts[0] ne '0/0' and $gts[0] ne './.' and $gts[0] ne '.|.' and $gts[0] ne '0|0')
	    	            		{
	    	            			if (exists $samples { "$chr;$position;$ref;$al" })
									{
										$samples { "$chr;$position;$ref;$al" } = $samples { "$chr;$position;$ref;$al" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf; 
									}else{
										$samples { "$chr;$position;$ref;$al" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
									}
								}
							}
						}

					}
				}
			}
		}else ## section for bi-allelic sites
		{
			## decide for variant type
			my $lenref = length $ref;
			my $lenalt = length $alt;
			
			if (($lenref + $lenalt) > 2)## INDEL
			{
				## make indelID
				my $token_ref = unpack('L', md5($ref));
				my $token_obs = unpack('L', md5($alt));
				
				if ($type eq 'INDEL' or $type eq 'all')
                {
					$variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$alt";
				
					for(my $i = 9; $i < @line; $i++)
					{
			    		my ($dpref,$dpalt,$samAf);
			    		my @gts = ();
			    		
			    		if ($line[$i] =~ m/\:/)
			    		{
							@gts = split(/\:/,$line[$i]);
							if ($gts[1] =~ m/\,/)
							{
	 		           			($dpref,$dpalt) = split(/\,/,$gts[1]);
							}else{
								$dpref = $gts[1];
								$dpalt = ".";
							}		    		
			    		    ## for missing genotype or homozygous reference genotype set the AF to 0
	            			if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
    	        			{
        	    				$samAf = "0";
            				}else{
                				$samAf = $AF;
							}

			    		}else{
			    			$gts[0] = $line[$i];
							$dpref = ".";
							$dpalt = "."; 

			    		    ## for missing genotype or homozygous reference genotype set the AF to 0
	            			if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
    	        			{
        	    				$samAf = "0";
            				}else{
                				$samAf = $AF;
							}

			    		}                		
						
						## check for sample genotype mode
						if ($gtMode eq "complete")
						{
	                		if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
    	            		{
        	        		    $samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
            		    	}else
        	    	    	{        
    	            		   	$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
   		         			}
   		         		}else ## compact
   		         		{
    	            		## kick out genotypes of '0/0','./.','0|0' and '.|.'
    	            		if ($gts[0] ne '0/0' and $gts[0] ne './.' and $gts[0] ne '.|.' and $gts[0] ne '0|0')
    	            		{
	                			if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
    	            			{
        	        			    $samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
            		    		}else
        	    	    		{        
    	            			   	$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
   		         				}
							}
   		         		}		
					}	
				
				}

			}else{ ## SNP
                if ($type eq 'SNP' or $type eq 'all')
                {
					$variants{ "$chr;$position;$ref;$alt" } = "$chr;$position;$ref;$alt";
					
					for(my $i = 9; $i < @line; $i++)
					{
						my ($dpref,$dpalt,$samAf);
						my @gts = ();
						
						if ($line[$i] =~ m/\:/)
						{
							@gts = split(/\:/,$line[$i]);
        	    			
        	    			if ($gts[1] =~ m/\,/)
        	    			{
	        	    			($dpref,$dpalt) = split(/\,/,$gts[1]);
						    }else{
								$dpref = $gts[1];
								$dpalt = "."; 
						    }
							
							## for missing genotype or homozygous reference genotype set the AF to 0
        	    			if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            				{
            					$samAf = "0";
        	    			}else{
    	            			$samAf = $AF;
							}

						}else{
							$gts[0] = $line[$i];
							$dpref = ".";
							$dpalt = "."; 

						    ## for missing genotype or homozygous reference genotype set the AF to 0
            				if ($gts[0] eq './.' or $gts[0] eq '0/0' or $gts[0] eq '.|.' or $gts[0] eq '0|0')
            				{
            					$samAf = "0";
        	    			}else{
    	            			$samAf = $AF;
							}

						}						
						
						## check for sample genotype mode
						if ($gtMode eq "complete")
						{
	                		if (exists $samples { "$chr;$position;$ref;$alt" })
    	            		{
        	        		    $samples { "$chr;$position;$ref;$alt" } = $samples { "$chr;$position;$ref;$alt" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
            	   			}else
            		    	{        
        	        		   	$samples { "$chr;$position;$ref;$alt" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
    	            		}
    	            	}else ## compact
    	            	{
    	            		## kick out genotypes of '0/0','./.','0|0' and '.|.'
    	            		if ($gts[0] ne '0/0' and $gts[0] ne './.' and $gts[0] ne '.|.' and $gts[0] ne '0|0')
    	            		{
	                			if (exists $samples { "$chr;$position;$ref;$alt" })
    	            			{
        	    	    		    $samples { "$chr;$position;$ref;$alt" } = $samples { "$chr;$position;$ref;$alt" }.";".$headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
            		   			}else
        	    		    	{        
    	    	        		   	$samples { "$chr;$position;$ref;$alt" } = $headers[$i].">".$gts[0].">".$dpref.">".$dpalt.">".$samAf;
	    	            		}	
    	            		}
    	            	}	
					}
				
				}	
			}
			
		} ## end of if-else on process based on alteration

	} ## end if
} ## end while
close(INPUT);






## Initialization completed
print "MESSAGE :: Finished processing input VCF file - $input \n";

## start threading and annotating
print "MESSAGE :: Annotation starting \n";

## spawn threds and push to list
push @thrds, threads -> new(\&eDiVaAnnotation); ## spawn a thread for eDiVa annotation
push @thrds, threads -> new(\&AnnovarAnnotation); ## spawn a thread for Annovar annotation

## join spawned threads
foreach my $thr (@thrds)
{
        $thr -> join();
}

## write annotaton in file 
print "MESSAGE :: Writing annotation to output file \n";


## open file handler
open (ANN, ">>".$outFile) or die "Cant open new file to write annotation \n";

## write header to output file
my $headerOutputFile = &getHeader;
print ANN  $headerOutputFile."\n";

## write data lines to main output file
while (my($key, $value) = each(%variants)) 
{
	my ($edivaannotationtoprint,$annovarannotationtoprint) = ("NA","NA");
	my ($chr,$position,$ref,$alt) = split(/\;/, $value);
	$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
	$annovarannotationtoprint = $Annovar{$value} if $Annovar{$value};
	
	## write annotation to file 
	print ANN $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$samples { $key }."\n";
}

## close the handler
close(ANN);

## check for inconsistent file presence, take action and open file handler
if (-e 	$outFileIns)
{
	unlink($outFileIns);
}
open (ANNINS, ">>".$outFileIns) or die "Cant open new file to write annotation \n";

## write header for inconsistent file
$headerOutputFile = &getHeaderIns;
print ANNINS  $headerOutputFile."\n";

## write data lines to main output file
while (my($key, $value) = each(%not_biallelic_variants)) 
{
	my ($edivaannotationtoprint,$annovarannotationtoprint) = ("NA","NA");
	my ($chr,$position,$ref,$alt) = split(/\;/, $value);
	$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};

	## write annotation to file 
	print ANNINS $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$samples { $key }."\n";	
}

## close the handler
close(ANNINS);


## sort the file
my $srtCmm = "sort -n -k1,1 -n -k2,2 --field-separator=, $outFile > $SortedoutFile ";
system($srtCmm);

## writing completed
print "MESSAGE :: Writing annotation completed \n";
print "MESSAGE :: Your annotated file is $outFile \n";
print "MESSAGE :: Your sorted annotated file is $SortedoutFile \n";
print "MESSAGE :: Reported sites which are not bi-allelic are in $outFileIns \n";

## Finalize everything
print "MESSAGE :: Finalizing annotation process \n";
&finalize;
print "MESSAGE :: Finalization completed \n";