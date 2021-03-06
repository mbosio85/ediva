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
sub usage { print "\n$0 usage:\n",
	   "--input,-i \t\t VCF file containing the variants to annoate \n\n",
	   "--quicklookup,-q \t Option : Command line variant argument. Variant format => chr\:position\:reference\:alteration \n",
	   "\t\t\t Option : File containing variants (one variant per each line). Variant format => chr\:position\:reference\:alteration \n\n",
	   "--tempDir,-t \t\t Temporary scratch location for temp files (all temp files with the current session will be removed at the end of execution) \n\t\t\t default: input VCF location \n\n",
	   "--geneDef,-g \t\t Gene deifnition you want to select for genic annotation (ensGene,refGene,knownGene,all) \n\t\t\t default: refGene \n\n",
	   "--variantType,-v \t Type of variants to annotate from input VCF file (SNP,INDEL,all) \n\t\t\t default: all \n\n",
	   "--sampleGenotypeMode,-s  complete: reports all genotypes from the input VCF file\n\t\t\t compact: reports only heterozygous and homozygous alteration genotypes from the input VCF file \n\t\t\t none: exclude sample wise genotype information in annotation \n\t\t\t default: compact \n\n",
	   "--onlyGenicAnnotation,-o  If set, then only genic annotation will be performed \n\n",
	   "--forceNewFileCreate,-f  If set, then it will over-write existing output annotation file with the same name \n\n",
	   "--help,-h \t\t show help \n\n"
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
my $onlygenic = 0; ## variable for only genic annotation 
my $forceDel = 0; ## varibale for force deleting the output annotation file (if exists)
my $qlookup = "NA"; ## varibale for enabling quick lookup mode of the program
our $templocation = "INPATH"; ## scratch place for creating the temp files while annotating

my %variants = (); ## hash to hold input variants from VCF
my %not_biallelic_variants = (); ## hash to hold input variants from VCF where the site is not bi-allelic and annovar annotation is absent
my @thrds = (); ## list to hold threads

our %eDiVa :shared = (); ## hash to hold input INDELs with eDiVa annotation
our %Annovar :shared = (); ## hash to hold input INDELs with ANNOVAR annotation
our %samples :shared = (); ## hash to hold sample information
our %edivaStr :shared = (); ## hash to hold simple tandem repeat data from ediva public omics database
our @headers = ();

## for missing db annotation for snps and indels
our $missanndb :shared = "NA";
our $missanndbindel :shared = "NA";
our $missanndb_coordinate :shared = "0";

my $fileSuffix = localtime();
$fileSuffix =~ s/\:/-/g;
$fileSuffix =~ s/\s/-/g;

## ANNOVAR settings
our $ANNOVAR = "/users/GD/tools/eDiVaCommandLine/lib/Annovar"; 

##############################################################################################
## COMMAND LINE OPTIONS
##############################################################################################


## grab command line options
unknownArguments() if (!GetOptions("input=s" => \$input, "tempDir=s" => \$templocation, "geneDef=s" => \$geneDef, "variantType=s" => \$type, "onlyGenicAnnotation" => \$onlygenic, "forceNewFileCreate" => \$forceDel, "quicklookup=s" => \$qlookup, "sampleGenotypeMode=s" => \$gtMode, "help" => \$help));


## check mandatory command line parameters and take necessary actions
unless(($input || $qlookup ne "NA") && $help == 0)
{
	usage;
	exit 0;
}

## check on both input file and quick lookup mode
if ($input and $qlookup ne "NA")
{
	print "\nERROR :: Either run the program with a VCF file or by using quick lookup command line arguments \n";
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
if ($gtMode ne "complete" && $gtMode ne "compact" && $gtMode ne "none")
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

## final check on quick lookup mode parameter
if ($qlookup ne "NA")
{
	if (-e $qlookup)
	{
		## nothing to do ! its a file !! varinat format inside the file is checked later on !!
	}else{
		if ($qlookup !~ m/\:/)
		{
			print "\nERROR :: Not a valid format. Correct format is chr:position:reference:alteration \n";
			exit 0;		
		}		
	}
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
## my $dbConn = DBI->connect('dbi:mysql:'.'eDiVa_innoDB'.';host=www.ediva.crg.eu','hana','hanamysql2013') or die "\nERROR :: Could not connect to eDiVa Database Server \n";


##############################################################################################
## OUTPUT FILE(s)
##############################################################################################

my ($outFile,$SortedoutFile,$outFileIns); ## files to write annotation

if ($qlookup eq "NA")
{
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

}else{
	
	if (-e $qlookup)
	{
		if ($qlookup =~ m/\//) ## input file with full path
		{
			my @paths = split(/\//,$qlookup);
			my $len = scalar @paths;
			my $mainOutpath = join("/",@paths[0..($len-2)]);

			my @files = split(/\./,$paths[$len-1]);
			$len = scalar @files;
			my $tempfile = join(".",@files[0..($len-2)]);
			$outFile = $mainOutpath."/".$tempfile.".annotated";
			$SortedoutFile = $mainOutpath."/".$tempfile.".sorted.annotated";
			#$outFileIns = $mainOutpath."/".$tempfile.".inconsistent.annotated";
	
			## set temp location
			#if ($templocation eq "INPATH")
			#{
			#	$templocation = $mainOutpath;
			#}

		}else{ ## just filename

			my @files = split(/\./,$qlookup);
			my $len = scalar @files;
			my $tempfile = join(".",@files[0..($len-2)]);
			$outFile = $tempfile.".annotated";
			$SortedoutFile = $tempfile.".sorted.annotated";
			#$outFileIns = $tempfile.".inconsistent.annotated";

			## set temp location
			#if ($templocation eq "INPATH")
			#{
			#	$templocation = ".";
			#}
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
	
	}

}
	

##############################################################################################
## SUBROUTINES
##############################################################################################

## subroutine for unknown command line arguments
sub unknownArguments
{
	print "\nERROR :: Unknown command line arguments were found \nPlease see the following usage \n";
	usage;
	exit 0;

}

## subroutine for finalizing annotation process
sub finalize
{
	## clear the tmp directory for this session
	my $clearCmm = "rm -r $templocation/*".$fileSuffix."*";
	system ($clearCmm);
}


## sub for preparing missing db annotation
sub preparemissdb
{
	## for snps
	for(my $i = 0; $i<9; $i++)
	{
		$missanndb = $missanndb.$sep."0";
	}

	for(my $i = 0; $i<14; $i++)
	{
		$missanndb = $missanndb.$sep."NA";
	}
	
	## for genomic co-ordinates for indels
	for(my $i = 0; $i<8; $i++)
	{
		$missanndb_coordinate = $missanndb_coordinate.$sep."NA";
	}
	
	## for indels
	for(my $i = 0; $i<8; $i++)
	{
		$missanndbindel = $missanndbindel.$sep."0";
	}
	
}

## sub for replacing commas inside double qoutes for annovar genic annotation lines
sub replaceCommainQoute
{
	my $newstr;
	my $str = shift;
	my @stra = split(//,$str);
	my $cnt = 0;

	foreach my $char (@stra)
	{
	        if ($char eq "\"")
	        {
	                $cnt = $cnt + 1;
	        }
	        if (($cnt%2) != 0)
	        {
	                $char =~ s/\,/\;/g;
	        }
	        $newstr = $newstr.$char;
	}
	
	return $newstr;
} 


## subroutine for eDiVa public omics data fetch
sub eDiVaPublicOmics
{
	## DB parameters
	my $username = 'edivacrg';
	my $database = 'eDiVa_public_omics';
	my $dbhost = 'mysqlsrv-ediva.linux.crg.es';
	my $pass = 'FD5KrT3q';
		
	## open DB connection
	my $dbh = DBI->connect('dbi:mysql:'.$database.';host='.$dbhost.'',$username,$pass) or die "Connection Error!!\n";

	my $sql = "select chr,pos,lengthofrepeat,copyNum,region from eDiVa_public_omics.Table_simpleRepeat\;";
	
	## prepare statement and query
	my $stmt = $dbh->prepare($sql);
	$stmt->execute or die "SQL Error!!\n";
	
	#process query result
	while (my @res = $stmt->fetchrow_array) 
	{
		$edivaStr{ $res[0].";".$res[1] } = $res[4].",".$res[2]; 
	}

   	## close DB connection
	$dbh->disconnect();
}


## subroutine for eDiVa annotation
sub eDiVaAnnotation
{
	## DB parameters
	my $username = 'edivacrg';
	my $database = 'eDiVa_annotation';
	my $dbhost = 'mysqlsrv-ediva.linux.crg.es';
	my $pass = 'FD5KrT3q';
		
	## open DB connection
	my $dbh = DBI->connect('dbi:mysql:'.$database.';host='.$dbhost.'',$username,$pass) or die "Connection Error!!\n";

	# extract db result
	while( my ($k, $v) = each %variants ) 
	{
		my ($chr,$pos,$ref,$alt) = split(/\;/,$k);

		## decide for variant type
		my $lenref = length $ref;
		my $lenalt = length $alt;
			
		if (($lenref + $lenalt) > 2)## INDEL
		{
			my $sql = "";
			my $sql2 = "";
			my $stmt ;
			my $stmt2 ;
			my @res = ();
			my @res2 = ();

			## prepare query 
			#$sql = "select annotateINDEL('$k','$chr',$pos);";
			$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0') from Table_Chr$chr\_indel where indelid = '$k' limit 1;";
			
			$sql2 = "select ifnull(SDIndel,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA') 
			from Table_Chr$chr where position = $pos limit 1;";

			$sql =~ s/\n//g;
			$sql2 =~ s/\n//g;
			
			## prepare statement and query
			$stmt = $dbh->prepare($sql);
			$stmt2 = $dbh->prepare($sql2);

			$stmt->execute or die "SQL Error!!\n";
			$stmt2->execute or die "SQL Error!!\n";	
			
			#process query result
			@res = $stmt->fetchrow_array;
			@res2 = $stmt2->fetchrow_array;
			
			if (scalar @res > 1 and scalar @res2 > 1) ## both returned database rows
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}				
				for my $reselement (@res2)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
				## add NAs for damage potential scores and cadd scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";


			}elsif(scalar @res > 1 and scalar @res2 < 1) ## only indel table returned database row
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
				
				# take care of missing positional values			
				$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndb_coordinate;
				## add NAs for damage potential scores and cadd scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";
			
			}elsif(scalar @res < 1 and scalar @res2 > 1) ## only snp table returned database row
			{
				# take care of missing positional values			
				if (exists $eDiVa{ $k })
				{
					## this should never happen !!
					$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndbindel;
					
				}else{
					$eDiVa{ $k } = $missanndbindel;
				}		
				
				for my $reselement (@res2)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}

				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";
							
			} 
			else{
				
				## no entry in the database for both the queries
				if (exists $eDiVa{ $k })
				{
					## this should never happen !!
					$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndbindel;
					
				}else{
					$eDiVa{ $k } = $missanndbindel;
				}		
				
				# take care of missing positional values			
				$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndb_coordinate;
				
				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";				
			}		


		}else{ ##SNP

			my $sql = "";
			my $stmt ;
			my @res = ();

			#$sql = "select annotateSNPGermline('$chr',$pos,'$ref','$alt');";					
			$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
			ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
			from Table_Chr$chr where position = $pos and Reference = '$ref' and Alt = '$alt' limit 1;";
			
			$sql =~ s/\n//g;
			
			## prepare statement and query
			$stmt = $dbh->prepare($sql);
			$stmt->execute or die "SQL Error!!\n";
	
			#process query result
			@res = $stmt->fetchrow_array;
			
			if (scalar @res > 1)
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
			}else{
				$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
				ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
				ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
				ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') from Table_Chr$chr where position = $pos limit 1;";
			
				$sql =~ s/\n//g;
			
				## prepare statement and query
				$stmt = $dbh->prepare($sql);
				$stmt->execute or die "SQL Error!!\n";

				#process query result
				@res = $stmt->fetchrow_array;
				
				if (scalar @res > 1)
				{
					# load eDiVa hash from database
					for my $reselement (@res)
					{
						if (exists $eDiVa{ $k })
						{
							$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
						}else{
							$eDiVa{ $k } = $reselement;
						}
					}
				}else{
					## handle missing database annotation entry
					$eDiVa{ $k } = $missanndb;
				}
			}
		} ## end of if-else for variant decision
		
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
			my $sql = "";
			my $sql2 = "";
			my $stmt ;
			my $stmt2 ;
			my @res = ();
			my @res2 = ();

			## prepare query 
			#$sql = "select annotateINDEL('$k','$chr',$pos);";
			$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0') from Table_Chr$chr\_indel where indelid = '$k' limit 1;";
			
			$sql2 = "select ifnull(SDIndel,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA') 
			from Table_Chr$chr where position = $pos limit 1;";

			$sql =~ s/\n//g;
			$sql2 =~ s/\n//g;
			
			## prepare statement and query
			$stmt = $dbh->prepare($sql);
			$stmt2 = $dbh->prepare($sql2);

			$stmt->execute or die "SQL Error!!\n";
			$stmt2->execute or die "SQL Error!!\n";	
			
			#process query result
			@res = $stmt->fetchrow_array;
			@res2 = $stmt2->fetchrow_array;
			
			if (scalar @res > 1 and scalar @res2 > 1) ## both returned database rows
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}				
				for my $reselement (@res2)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";


			}elsif(scalar @res > 1 and scalar @res2 < 1) ## only indel table returned database row
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
				
				# take care of missing positional values			
				$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndb_coordinate;
				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";
			
			}elsif(scalar @res < 1 and scalar @res2 > 1) ## only snp table returned database row
			{
				# take care of missing positional values			
				if (exists $eDiVa{ $k })
				{
					## this should never happen !!
					$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndbindel;
					
				}else{
					$eDiVa{ $k } = $missanndbindel;
				}		
				
				for my $reselement (@res2)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}

				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";
							
			} 
			else{
				
				## no entry in the database for both the queries
				if (exists $eDiVa{ $k })
				{
					## this should never happen !!
					$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndbindel;
					
				}else{
					$eDiVa{ $k } = $missanndbindel;
				}		
				
				# take care of missing positional values			
				$eDiVa{ $k } = $eDiVa{ $k }.$sep.$missanndb_coordinate;
				
				## add NAs for damage potential scores for indels
				$eDiVa{ $k } = $eDiVa{ $k }.$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA".$sep."NA";				
			}		

		}else{ ##SNP
		
			my $sql = "";
			my $stmt ;
			my @res = ();

			#$sql = "select annotateSNPGermline('$chr',$pos,'$ref','$alt');";					
			$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
			ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
			ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
			ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') 
			from Table_Chr$chr where position = $pos and Reference = '$ref' and Alt = '$alt' limit 1;";
			
			$sql =~ s/\n//g;
			
			## prepare statement and query
			$stmt = $dbh->prepare($sql);
			$stmt->execute or die "SQL Error!!\n";
	
			#process query result
			@res = $stmt->fetchrow_array;
			
			if (scalar @res > 1)
			{
				# load eDiVa hash from database
				for my $reselement (@res)
				{
					if (exists $eDiVa{ $k })
					{
						$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
					}else{
						$eDiVa{ $k } = $reselement;
					}
				}
			}else{
				$sql = "select ifnull(dbsnpid,'NA'),ifnull(EurEVSFreq,'0'),ifnull(AfrEVSFreq,'0'),ifnull(TotalEVSFreq,'0'),ifnull(EurAFKG,'0'),ifnull(AfrAFKG,'0'),
				ifnull(AsaAFKG,'0'),ifnull(AmrAFKG,'0'),ifnull(AFKG,'0'),ifnull(SDSnp,'0'),ifnull(`placentalMammal.phyloP`,'NA'),ifnull(`primates.phyloP`,'NA'),ifnull(`vertebrates.phyloP`,'NA'),
				ifnull(`placentalMammal.phastCons`,'NA'),ifnull(`primates.phastCons`,'NA'),ifnull(`vertebrates.phastCons`,'NA'),ifnull(gerp1,'NA'),ifnull(gerp2,'NA'),
				ifnull(sift,'NA'),ifnull(polyphen2,'NA'),ifnull(mutationassessor,'NA'),ifnull(condel,'NA'),ifnull(cadd1,'NA'),ifnull(cadd2,'NA') from Table_Chr$chr where position = $pos limit 1;";
			
				$sql =~ s/\n//g;
			
				## prepare statement and query
				$stmt = $dbh->prepare($sql);
				$stmt->execute or die "SQL Error!!\n";

				#process query result
				@res = $stmt->fetchrow_array;
				
				if (scalar @res > 1)
				{
					# load eDiVa hash from database
					for my $reselement (@res)
					{
						if (exists $eDiVa{ $k })
						{
							$eDiVa{ $k } = $eDiVa{ $k }.$sep.$reselement;
						}else{
							$eDiVa{ $k } = $reselement;
						}
					}
				}else{
					## handle missing database annotation entry
					$eDiVa{ $k } = $missanndb;
				}
			}		
		}
		
    } ## end of while on not_biallelic_variants

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
				my $newAnnovarLine = &replaceCommainQoute($_);
				my @dt = split(/\,/,$newAnnovarLine);
								
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\;/)
				{
					my @annalts = split(/\;/,$dt[30]);
					$dt[30] = $annalts[0];
				}
                                $dt[26] =~ s/chr//g;
                                $dt[26] =~ s/Chr//g;
				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
				$valueTOmatch =~ s/\"//g;

				## fix missing values
				if ($dt[0] eq '')
				{
					$dt[0] = 'NA';
				}
						## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
				
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
				my $newAnnovarLine = &replaceCommainQoute($_);
				my @dt = split(/\,/,$newAnnovarLine);
				
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\;/)
				{
					my @annalts = split(/\;/,$dt[30]);
					$dt[30] = $annalts[0];
				}
                                $dt[26] =~ s/chr//g;
                                $dt[26] =~ s/Chr//g;
				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
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
				my $newAnnovarLine = &replaceCommainQoute($_);
				my @dt = split(/\,/,$newAnnovarLine);
				
				$dt[30] =~ s/\"//g;
				if ($dt[30] =~ m/\;/)
				{
					my @annalts = split(/\;/,$dt[30]);
					$dt[30] = $annalts[0];
				}
                                $dt[26] =~ s/chr//g;
                                $dt[26] =~ s/Chr//g;
				my $valueTOmatch = $dt[26].";".$dt[27].";".$dt[29].";".$dt[30];
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
    
	if ($onlygenic == 1)
	{
		## only genic annotation header
		## check for gene definiton and construct header according to that
	    if($geneDef eq 'ensGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
	        AminoAcidChange(Ensembl),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }elsif($geneDef eq 'refGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),
	        NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }elsif($geneDef eq 'knownGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }else{
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
	        AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),NA,NA,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }
		
	}else{
	    ## normal annotation header
		## check for gene definiton and construct header according to that
	    if($geneDef eq 'ensGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
	        AminoAcidChange(Ensembl),dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
	        Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
	        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }elsif($geneDef eq 'refGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),
	        dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
	        Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
	        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }elsif($geneDef eq 'knownGene')
	    {
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,
	        TotalEVSFrequency,Eur1000GenomesFrequency,Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
	        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }else{
	        $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
	        AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,EurEVSFrequency,AfrEVSFrequency,TotalEVSFrequency,Eur1000GenomesFrequency,
	        Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
	        PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";
	    }		
	}    
	
    ## replace newlines with nothing at header line
    $stringTOreturn =~ s/\n|\s+//g;
    return $stringTOreturn;
}



## subroutnine por providing header to the inconsistent annotation output file
sub getHeaderIns
{
    my $stringTOreturn; ## header to return

    $stringTOreturn = "Chr,Position,Reference,Alteration,AlleleFrequency,GenicAnnotation,dbsnpIdentifier,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
    Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
    PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength,samples(sampleid>zygosity>DPRef>DPAlt>AlleleFraction)";

    ## replace newlines with nothing at header line
    $stringTOreturn =~ s/\n|\s+//g;
    return $stringTOreturn;
}


## subroutnine por providing header to the quick look up mode annotation output file
sub getHeaderQlookup
{
    my $stringTOreturn; ## header to return

    $stringTOreturn = "Chr,Position,Reference,Alteration,dbsnpIdentifier,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
    Afr1000GenomesFrequency,Asia1000GenomesFrequency,Amr1000GenomesFrequency,Total1000GenomesFrequency,SegMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
    PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTScore,polyphen2,MutAss,Condel,Cadd1,Cadd2,SimpleTandemRepeatRegion,SimpleTandemRepeatLength";

    ## replace newlines with nothing at header line
    $stringTOreturn =~ s/\n|\s+//g;
    return $stringTOreturn;
}



########################################################################################################################
## MAIN starts
########################################################################################################################

if ($qlookup eq "NA")
{
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
		
		## check for malformed VCF file and take action
		if (scalar @headers < 8)
		{
			print "ERROR:: Not a valid VCF format \n";
			exit 0;			
		}
		
		## check for genotype format and sample information in the VCF and genotype mode parameter value and take action
		if (scalar @headers < 9 and $gtMode ne "none")
		{
			print "ERROR:: No Genptype format column and no sample genotype information column present in the VCF. Please run the tool with the --sampleGenotypeMode parameter set to \"none\" \n";
			exit 0;						
		}
	}
	else ## data lines
	{
		my @line = split(/\t/,$_);

		## check for malformed VCF file and take action
		if (scalar @line < 8)
		{
			print "ERROR:: Not a valid VCF format \n";
			exit 0;			
		}

		my $chr = $line[0];
		my $position = $line[1];
		my $ref = $line[3];
		my $alt = $line[4];
		my @infos = split(/\;/,$line[7]);
		my $AF = ".";			
		my $gtindex = "NF";
		my $adindex = "NF";
		print $alt."\n";
		print $alt !~ m/[ACGT]/ . "\n";
		## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
		if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
		{
			$chr = substr($chr,3);
		}

		## take care of chr 23 or 24 and convert to X or Y
		if ($chr eq "23")
		{
			$chr = "X";
		}

		if ($chr eq "24")
		{
			$chr = "Y";
		}

		if ($chr eq "25")
		{
			$chr = "MT";
		}

		if ($chr =~ m/[MTmt]/ or $chr =~ m/\_/)
		{
			print "WARNING:: $0 does not support chromosome $chr currently. This variant line will be skipped in the final annotation output \n";			
		}
		elsif($alt !~ m/[ACGT]/)
		{
			print "WARNING:: Unknown alternate allele detected at $chr and $position. This variant line will be skipped in the final annotation output \n";
		}
		else
		{
		## grab the AF from the INFO field in the VCF
		foreach my $info (@infos)
		{
			if ($info =~ m/^AF=/)
			{
				$AF = substr($info,3);
				last;
			}
		}

		## confirm AF extraction from the info tags
		if ($AF eq ".")
		{
			print "WARNING:: AF tag not found in the INFO column at chromosome $chr and position $position. AF will be set to \".\" for this variant \n";
		}

		## always test for complete genotype format field consistency in the VCF; if abnormal report for that variant
		if ($gtMode ne "none")
		{	
			if (scalar @line > 8 and $line[8] =~ m/\:/)
			{
				my @gtcheck = split(/\:/,$line[8]);
				
				for(my $gti = 0; $gti < @gtcheck; $gti++)
				{
					if ($gtcheck[$gti] eq "GT")
					{
						$gtindex = $gti;
					}
					if ($gtcheck[$gti] eq "AD")
					{
						$adindex = $gti;
					}
				}
				
				## check for the GT and AD fields
				if ($gtindex eq "NF")
				{
					print "WARNING:: Weird genotype format $line[8] found in $input at chromosome $chr and position $position. No GT field present in $line[8] \n";
				}
				if ($adindex eq "NF")
				{
					print "WARNING:: Weird genotype format $line[8] found in $input at chromosome $chr and position $position. No AD field present in $line[8] \n";
				}
				
			}else{
				print "WARNING:: Weird genotype format $line[8] found in $input at chromosome $chr and position $position. Expected genotype format field separator \":\" \n";
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
					my $token_ref = "NA";
					my $token_obs = "NA";
					
					## make indelID
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
					if ($ref !~ m/[Nn]/ and $al !~ m/[Nn]/)
					{
						$token_ref = unpack('L', md5($ref));
						$token_obs = unpack('L', md5($al));
					}
					
                	if ($type eq 'INDEL' or $type eq 'all')
                	{
        	        	## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
    	            	## e.g A,C in the alternate column in VCF will report only A in the main annotation file
	                	## we are doing this because we want to keep the annotation main file consistent

	                	if ($j == 0 and $token_ref ne "NA" and $token_obs ne "NA")
    	            	{
							$variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$al;$alfr";
						}else{
							$not_biallelic_variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$al;$alfr";
						}

						## if sample wise information is present in the VCF then process ; otherwise skip
						## also check for none value in genotype mode parameter
						if (scalar @line > 8 and $gtMode ne "none")
						{
							for(my $i = 9; $i < @line; $i++)
							{
								my ($genotype,$dpref,$dpalt,$samAf);
								my @gts = ();
							
								if ($line[$i] =~ m/\:/)
								{
									@gts = split(/\:/,$line[$i]);
								
									if ($adindex ne "NF" and $gts[$adindex] =~ m/\,/)
									{
	    	            			    my @ads = split(/\,/,$gts[$adindex]);
    	    	        	 			$dpref = $ads[0];
        	    	   		 			@ads = @ads[1..(scalar @ads -1 )];
            	    				   	$dpalt = $ads[$j];

		        	    				if ($dpref+$dpalt > 0)
		        	    				{
			        	    				$samAf = $dpalt/($dpref+$dpalt);
		    	    	    				$samAf = substr($samAf,0,5);
										}else{
											$samAf = "0.000";
										}

            	    			    }else{
										$dpref = ".";
										$dpalt = ".";
										$samAf = "."; 
									}
            	    				
            	    				if ($gtindex ne "NF")
            	    				{
            	    					$genotype = $gts[$gtindex];
									}else{
									    $genotype = ".";
									}			    									    
								}else{
				    				$genotype = $line[$i];
									$dpref = ".";
									$dpalt = "."; 
							    	$samAf = ".";   
								}							
						
								## check for sample genotype mode
								if ($gtMode eq "complete")
								{					
               				 		if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
                					{
	            	    					$samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
    								}else
        		        			{    
        	   	    	  	 			$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
   			         				}
   			         			}else ## compact
   			         			{
   			         				## kick out genotypes of '0/0','./.','0|0' and '.|.'
    		            			if ($genotype ne '0/0' and $genotype ne './.' and $genotype ne '.|.' and $genotype ne '0|0')
	    	            			{
	    	            				if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
                						{
	            	    					$samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
    									}else
        		        				{    
        	   	    	   					$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
   			         					}
									}
   			         			}	
							}
						} ## end of sample wise information checking IF	

					}

				}else{ ## SNP
					
                	## we are only going to report the first alternate allele in the cases where the site is more than bi-allelic
                	## e.g A,C in the alternate column in VCF will report only A in the main annotation file
                	## we are doing this because we want to keep the annotation main file consistent
                	if ($type eq 'SNP' or $type eq 'all')
                	{
						## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
	                	if ($j == 0 and $ref !~ m/[Nn]/ and $al !~ m/[Nn]/)
    	            	{
							$variants{ "$chr;$position;$ref;$al" } = "$chr;$position;$ref;$al;$alfr";
						}else{
							$not_biallelic_variants{ "$chr;$position;$ref;$al" } = "$chr;$position;$ref;$al;$alfr";
						}
					
						## if sample wise information is present in the VCF then process ; otherwise skip
						## also check for none value in genotype mode parameter
						if (scalar @line > 8 and $gtMode ne "none")
						{					
							for(my $i = 9; $i < @line; $i++)
							{
								my ($genotype,$dpref,$dpalt,$samAf);
								my @gts = ();

								if ($line[$i] =~ m/\:/)
								{
									@gts = split(/\:/,$line[$i]);
								
									if ($adindex ne "NF" and $gts[$adindex] =~ m/\,/)
									{
	    	            			    my @ads = split(/\,/,$gts[$adindex]);
    	    	        	 			$dpref = $ads[0];
        	    	   		 			@ads = @ads[1..(scalar @ads -1 )];
            	    				   	$dpalt = $ads[$j];

		        	    				if ($dpref+$dpalt > 0)
		        	    				{
			        	    				$samAf = $dpalt/($dpref+$dpalt);
		    	    	    				$samAf = substr($samAf,0,5);
										}else{
											$samAf = "0.000";
										}
										            	    				   	
            	    			    }else{
										$dpref = ".";
										$dpalt = "."; 
										$samAf = ".";								
									}
            	    				
            	    				if ($gtindex ne "NF")
            	    				{
            	    					$genotype = $gts[$gtindex];
									}else{
									    $genotype = ".";
									}			    									    
								}else{
				    				$genotype = $line[$i];
									$dpref = ".";
									$dpalt = "."; 
							    	$samAf = ".";   
								}							
                	    
							
								## check for sample genotype mode
								if ($gtMode eq "complete")
								{								
									if (exists $samples { "$chr;$position;$ref;$al" })
									{
										$samples { "$chr;$position;$ref;$al" } = $samples { "$chr;$position;$ref;$al" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf; 
									}else{
										$samples { "$chr;$position;$ref;$al" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
									}
								}else ## compact
								{
		   	            			## kick out genotypes of '0/0','./.','0|0' and '.|.'
    		            			if ($genotype ne '0/0' and $genotype ne './.' and $genotype ne '.|.' and $genotype ne '0|0')
	    	            			{
	    	            				if (exists $samples { "$chr;$position;$ref;$al" })
										{
											$samples { "$chr;$position;$ref;$al" } = $samples { "$chr;$position;$ref;$al" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf; 
										}else{
											$samples { "$chr;$position;$ref;$al" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
										}
									}
								}
							}
						} ## end of sample wise information checking of IF
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
				my $token_ref = "NA";
				my $token_obs = "NA";
				
				## make indelID
				## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section				
				if ($ref !~ m/[Nn]/ and $alt !~ m/[Nn]/)
				{
					$token_ref = unpack('L', md5($ref));
					$token_obs = unpack('L', md5($alt));
				}
				
				if ($type eq 'INDEL' or $type eq 'all')
                {
					
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
                	if ($token_ref ne "NA" and $token_obs ne "NA")
	            	{
						$variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$alt;$AF";
					}else{
						$not_biallelic_variants{ "$chr;$position;$token_ref;$token_obs" } = "$chr;$position;$ref;$alt;$AF";
					}
				
					## if sample wise information is present in the VCF then process ; otherwise skip
					## also check for none value in genotype mode parameter
					if (scalar @line > 8 and $gtMode ne "none")
					{						
						for(my $i = 9; $i < @line; $i++)
						{
							my ($genotype,$dpref,$dpalt,$samAf);
			    			my @gts = ();
			    		
							if ($line[$i] =~ m/\:/)
							{
								@gts = split(/\:/,$line[$i]);
        	    			
        	    				if ($adindex ne "NF" and $gts[$adindex] =~ m/\,/)
        	    				{
	        	    				($dpref,$dpalt) = split(/\,/,$gts[$adindex]);

	        	    				if ($dpref+$dpalt > 0)
	        	    				{
		        	    				$samAf = $dpalt/($dpref+$dpalt);
		        	    				$samAf = substr($samAf,0,5);
									}else{
										$samAf = "0.000";
									}
							    }else{
									$dpref = "."; #$gts[1];
									$dpalt = "."; 
									$samAf = ".";								
							    }
							
								## for missing genotype or homozygous reference genotype set the AF to 0
								if ($gtindex ne "NF")
								{
									$genotype = $gts[$gtindex];
								}else{
									$genotype = ".";
								}
							}else{
								$genotype = $line[$i];
								$dpref = ".";
								$dpalt = "."; 
								$samAf = ".";								
							}						
						
							## check for sample genotype mode
							if ($gtMode eq "complete")
							{
	                			if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
    	            			{
        	        			    $samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
            			    	}else
        	    		    	{        
    	            			   	$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
   		         				}
   		         			}else ## compact
   		         			{
    	            			## kick out genotypes of '0/0','./.','0|0' and '.|.'
    	        	    		if ($genotype ne '0/0' and $genotype ne './.' and $genotype ne '.|.' and $genotype ne '0|0')
    	        	    		{
	                				if (exists $samples{ "$chr;$position;$token_ref;$token_obs" })
    	            				{
        	        				    $samples{ "$chr;$position;$token_ref;$token_obs" } = $samples{ "$chr;$position;$token_ref;$token_obs" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
            			    		}else
        	    		    		{        
    	            				   	$samples{ "$chr;$position;$token_ref;$token_obs" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
   		         					}
								}
   		         			}		
						}
					} ## end of sample wise information checking of IF	
				
				}

			}else{ ## SNP
                if ($type eq 'SNP' or $type eq 'all')
                {					
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section
                	if ($ref !~ m/[Nn]/ and $alt !~ m/[Nn]/)
	            	{
						$variants{ "$chr;$position;$ref;$alt" } = "$chr;$position;$ref;$alt;$AF";
					}else{
						$not_biallelic_variants{ "$chr;$position;$ref;$alt" } = "$chr;$position;$ref;$alt;$AF";
					}
					
					## if sample wise information is present in the VCF then process ; otherwise skip
					## also check for none value in genotype mode parameter
					if (scalar @line > 8 and $gtMode ne "none")
					{										
						for(my $i = 9; $i < @line; $i++)
						{
							my ($genotype,$dpref,$dpalt,$samAf);
							my @gts = ();
							
							if ($line[$i] =~ m/\:/)
							{
								@gts = split(/\:/,$line[$i]);
        	    			
        	    				if ($adindex ne "NF" and $gts[$adindex] =~ m/\,/)
        	    				{
	        	    				($dpref,$dpalt) = split(/\,/,$gts[$adindex]);
	        	    				
	        	    				if ($dpref+$dpalt > 0)
	        	    				{
		        	    				$samAf = $dpalt/($dpref+$dpalt);
		        	    				$samAf = substr($samAf,0,5);
									}else{
										$samAf = "0.000";
									}
							    }else{
									$dpref = "."; #$gts[1];
									$dpalt = "."; 
									$samAf = ".";								
							    }
							
								## for missing genotype or homozygous reference genotype set the AF to 0
								if ($gtindex ne "NF")
								{
									$genotype = $gts[$gtindex];
								}else{
									$genotype = ".";
								}
							}else{
								$genotype = $line[$i];
								$dpref = ".";
								$dpalt = "."; 
								$samAf = ".";								
							}						
						
							## check for sample genotype mode
							if ($gtMode eq "complete")
							{
	                			if (exists $samples { "$chr;$position;$ref;$alt" })
    	            			{
        	        			    $samples { "$chr;$position;$ref;$alt" } = $samples { "$chr;$position;$ref;$alt" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
            	   				}else
            			    	{        
        	        			   	$samples { "$chr;$position;$ref;$alt" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
    	            			}
    	            		}else ## compact
    	            		{
    	        	    		## kick out genotypes of '0/0','./.','0|0' and '.|.'
    	        	    		if ($genotype ne '0/0' and $genotype ne './.' and $genotype ne '.|.' and $genotype ne '0|0')
    	        	    		{
	            	    			if (exists $samples { "$chr;$position;$ref;$alt" })
    	        	    			{
        	    		    		    $samples { "$chr;$position;$ref;$alt" } = $samples { "$chr;$position;$ref;$alt" }.";".$headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
            			   			}else
        	    			    	{        
    	    	    	    		   	$samples { "$chr;$position;$ref;$alt" } = $headers[$i].">".$genotype.">".$dpref.">".$dpalt.">".$samAf;
	    	        	    		}	
    	            			}
    	            		}	
						}
					} ## end of sample wise information checking of IF
				}	
			}
			
		} ## end of if-else on process based on alteration

	} ## end if-else on chrMt and others
	} ## end if on data lines
} ## end while
close(INPUT);


## Initialization completed
print "MESSAGE :: Finished processing input VCF file - $input \n";

}else{
	## Quick lookup mode 

	## check for input type in -q parameter
	## decide for file or query string
	if (-e $qlookup)
	{
		## its a file
		print "MESSAGE :: Processing input file - $qlookup \n";
		
		open (FL, $qlookup) or die "cant open input file $qlookup \n";
		
		while(<FL>)
		{
			chomp $_;
			my ($chr,$pos,$ref,$alt);
			my @var = split(/\:/,$_);
	
			## check for quick lookup data field consistency
			if (scalar @var != 4)
			{
				print "\nERROR :: Not a valid format. Correct format is chr:position:reference:alteration \n";
				exit 0;
			}else{
				## check for simple checking of the quick loookup data fields
				if ($var[1] =~ m/[a-zA-Z]/)
				{
					print "\nERROR :: Not a valid position value \n";
					exit 0;			
				}
				elsif ($var[2] =~ m/[0-9Nn]/)
				{
					print "\nERROR :: Not a valid reference allele value \n";
					exit 0;			
				}
				elsif ($var[3] =~ m/[0-9Nn]/)
				{
					print "\nERROR :: Not a valid alternate allele value \n";
					exit 0;			
				}else{
					## assign variation from quick look up format
					$chr = $var[0];
					$pos = $var[1];
					$ref = $var[2];
					$alt = $var[3];
	
					## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
					if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
					{
						$chr = substr($chr,3);
					}

					## take care of chr 23 or 24 and convert to X or Y
					if ($chr eq "23")
					{
						$chr = "X";
					}

					if ($chr eq "24")
					{
						$chr = "Y";
					}

					if ($chr eq "25")
					{
						$chr = "MT";
					}
			
					### prepare ediva hash
					## decide for variant type
					my $lenref = length $ref;
					my $lenalt = length $alt;
			
					if (($lenref + $lenalt) > 2) ## INDEL
					{
						my $token_ref = "NA";
						my $token_obs = "NA";
				
						## make indelID
						## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section				
						if ($ref !~ m/[Nn]/ and $alt !~ m/[Nn]/)
						{
							$token_ref = unpack('L', md5($ref));
							$token_obs = unpack('L', md5($alt));
						}
						$variants{ "$chr;$pos;$token_ref;$token_obs" } = "$chr;$pos;$ref;$alt;.";
					}else{ ## SNP
						$variants{ "$chr;$pos;$ref;$alt" } = "$chr;$pos;$ref;$alt;.";				
					}					
				} ## end of inner if - elsif -else

			}## end of outer if -else
		}
		
		close(FL);
		
	}else
	{
		## query string
		my ($chr,$pos,$ref,$alt);
		my @var = split(/\:/,$qlookup);
	
		## check for quick lookup data field consistency
		if (scalar @var != 4)
		{
			print "\nERROR :: Not a valid format. Correct format is chr:position:reference:alteration \n";
			exit 0;
		}else{
			## check for simple checking of the quick loookup data fields
			if ($var[1] =~ m/[a-zA-Z]/)
			{
				print "\nERROR :: Not a valid position value \n";
				exit 0;			
			}
			elsif ($var[2] =~ m/[0-9Nn]/)
			{
				print "\nERROR :: Not a valid reference allele value \n";
				exit 0;			
			}
			elsif ($var[3] =~ m/[0-9Nn]/)
			{
				print "\nERROR :: Not a valid alternate allele value \n";
				exit 0;			
			}else{
				## assign variation from quick look up format
				$chr = $var[0];
				$pos = $var[1];
				$ref = $var[2];
				$alt = $var[3];
	
				## take care of chr1 or Chr1 and convert to chr1/Chr1-> 1
				if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
				{
					$chr = substr($chr,3);
				}

				## take care of chr 23 or 24 and convert to X or Y
				if ($chr eq "23")
				{
					$chr = "X";
				}

				if ($chr eq "24")
				{
					$chr = "Y";
				}

				if ($chr eq "25")
				{
					$chr = "MT";
				}
			
				### prepare ediva hash
				## decide for variant type
				my $lenref = length $ref;
				my $lenalt = length $alt;
			
				if (($lenref + $lenalt) > 2) ## INDEL
				{
					my $token_ref = "NA";
					my $token_obs = "NA";
				
					## make indelID
					## if ref or alternate allele is N, then annovar fails to make genic annotation for them; so put them in inconsistent section				
					if ($ref !~ m/[Nn]/ and $alt !~ m/[Nn]/)
					{
						$token_ref = unpack('L', md5($ref));
						$token_obs = unpack('L', md5($alt));
					}
					$variants{ "$chr;$pos;$token_ref;$token_obs" } = "$chr;$pos;$ref;$alt;.";
				}else{ ## SNP
					$variants{ "$chr;$pos;$ref;$alt" } = "$chr;$pos;$ref;$alt;.";				
				}					
			} ## end of inner if - elsif -else

		}## end of outer if -else

	} ## end of input parameter deciding if-else
}

## prepare missing data handler for db annotation
&preparemissdb;

## spawn threds and push to list
if ($qlookup eq "NA")
{
	if ($onlygenic == 1)
	{
		## start a sigle thread for annovar genic annotation
		print "MESSAGE :: Annotation starting \n";
		push @thrds, threads -> new(\&AnnovarAnnotation); ## spawn a thread for Annovar annotation
	}else{
		## start threading and annotating
		print "MESSAGE :: Annotation starting \n";

		push @thrds, threads -> new(\&eDiVaAnnotation); ## spawn a thread for eDiVa annotation
		push @thrds, threads -> new(\&AnnovarAnnotation); ## spawn a thread for Annovar annotation
		push @thrds, threads -> new(\&eDiVaPublicOmics); ## spawn a thread for eDiVa public omics	
	}
}else{
	
	push @thrds, threads -> new(\&eDiVaAnnotation); ## spawn a thread for eDiVa annotation
	push @thrds, threads -> new(\&eDiVaPublicOmics); ## spawn a thread for eDiVa public omics		
	
}

## join spawned threads
foreach my $thr (@thrds)
{
        $thr -> join();
}

## write annotation to file or ender output
if ($qlookup eq "NA")
{
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
		my ($edivaannotationtoprint,$annovarannotationtoprint,$samplewiseinfortoprint) = ("NA","NA","NA");
		my $edivapublicanntoprint = "NA,NA";
		my ($chr,$position,$ref,$alt,$aftoprint) = split(/\;/, $value);
		my $annovarValueToMatch = $chr.";".$position.";".$ref.";".$alt;

		$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
		$annovarannotationtoprint = $Annovar{$annovarValueToMatch} if $Annovar{$annovarValueToMatch};
		$samplewiseinfortoprint = $samples { $key } if $samples { $key };
		$edivapublicanntoprint = $edivaStr { $chr.";".$position} if $edivaStr { $chr.";".$position}; 
	
		## write annotation to file 
		print ANN $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$aftoprint.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$edivapublicanntoprint.$sep.$samplewiseinfortoprint."\n";
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
		my ($edivaannotationtoprint,$annovarannotationtoprint,$samplewiseinfortoprint) = ("NA","NA","NA");
		my $edivapublicanntoprint = "NA,NA";
		my ($chr,$position,$ref,$alt,$aftoprint) = split(/\;/, $value);

		$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
		$samplewiseinfortoprint = $samples { $key } if $samples { $key };
		$edivapublicanntoprint = $edivaStr { $chr.";".$position} if $edivaStr { $chr.";".$position}; 

		## write annotation to file 
		print ANNINS $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$aftoprint.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$edivapublicanntoprint.$sep.$samplewiseinfortoprint."\n";	
	}

	## close the handler
	close(ANNINS);


	## sort the file
	my $srtCmm = "sort -k1,1 -n -k2,2 --field-separator=, $outFile > $SortedoutFile ";
	system($srtCmm);

	## writing completed
	print "MESSAGE :: Writing annotation completed \n";
	print "MESSAGE :: Your annotated file is $outFile \n";
	print "MESSAGE :: Your sorted annotated file is $SortedoutFile \n";
	print "MESSAGE :: Reported non bi-allelic sites are in $outFileIns \n";
	
	## Finalize everything
	print "MESSAGE :: Finalizing annotation process \n";
	&finalize;
	print "MESSAGE :: Finalization completed \n";

}else{

	## check for input type in -q parameter
	## decide for output file or command line output
	if (-e $qlookup)
	{	
		## render annotation to output file
		open (ANN,">>".$outFile) or die "Cant open new file \n";
		
		## get header 
		my $header = &getHeaderQlookup;
		
		## print header
		print ANN $header."\n";
		
		while (my($key, $value) = each(%variants)) 
		{
			my ($annovarannotationtoprint,$edivaannotationtoprint) = ("NA","NA");
			my $edivapublicanntoprint = "NA,NA";
			my ($chr,$position,$ref,$alt,$aftoprint) = split(/\;/, $value);
			#my $annovarValueToMatch = $chr.";".$position.";".$ref.";".$alt;

			$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
			$edivapublicanntoprint = $edivaStr { $chr.";".$position} if $edivaStr { $chr.";".$position}; 
			#$annovarannotationtoprint = $Annovar{$annovarValueToMatch} if $Annovar{$annovarValueToMatch};
		
			## write to file
			print ANN $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$edivaannotationtoprint.$sep.$edivapublicanntoprint."\n";
		}
		
		close(ANN);
		
		## sort the file
		my $srtCmm = "sort -k1,1 -n -k2,2 --field-separator=, $outFile > $SortedoutFile ";
		system($srtCmm);

		## writing completed
		print "MESSAGE :: Writing annotation completed \n";
		print "MESSAGE :: Your annotated file is $outFile \n";
		print "MESSAGE :: Your sorted annotated file is $SortedoutFile \n";
		
	}else{
		## render annotation to output
		while (my($key, $value) = each(%variants)) 
		{
			my ($edivaannotationtoprint,$annovarannotationtoprint) = ("NA","NA");
			my $edivapublicanntoprint = "NA,NA";
			my ($chr,$position,$ref,$alt,$aftoprint) = split(/\;/, $value);
			#my $annovarValueToMatch = $chr.";".$position.";".$ref.";".$alt;

			$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
			$edivapublicanntoprint = $edivaStr { $chr.";".$position} if $edivaStr { $chr.";".$position}; 
			#$annovarannotationtoprint = $Annovar{$annovarValueToMatch} if $Annovar{$annovarValueToMatch};
		
			my @edivavals = split(/\,/, $edivaannotationtoprint);
			my @edivapublicvals = split(/\,/, $edivapublicanntoprint);

			## render to command line output
			print "chromosome: $chr \n";
			print "position: $position \n";
			print "Reference: $ref \n";
			print "Alteration: $alt \n";
			print "dbSNP identifier: $edivavals[0] \n";
			print "EVS european frequency: $edivavals[1] \n";
			print "EVS african frequency: $edivavals[2] \n";
			print "EVS total frequency: $edivavals[3] \n";
			print "1000genomes european frequency: $edivavals[4] \n";
			print "1000genomes african frequency: $edivavals[5] \n";
			print "1000genomes american frequency: $edivavals[6] \n";
			print "1000genomes asian frequency: $edivavals[7] \n";
			print "1000genomes total frequency: $edivavals[8] \n";
			print "Segment duplication: $edivavals[9] \n";
			print "Placental mammal phyloP: $edivavals[10] \n";
			print "Primates phyloP: $edivavals[11] \n";
			print "Vertebrates phyloP: $edivavals[12] \n";
			print "Placental mammal phastcons: $edivavals[13] \n";
			print "Primates phastcons: $edivavals[14] \n";
			print "Vertebrates phastcons: $edivavals[15] \n";
			print "Gerp score1: $edivavals[16] \n";			
			print "Gerp score2: $edivavals[17] \n";
			print "Sift: $edivavals[18] \n";			
			print "polyphen2: $edivavals[19] \n";
			print "Mutationassessor: $edivavals[20] \n";			
			print "Condel: $edivavals[21] \n";
			print "Cadd score1: $edivavals[22] \n";			
			print "Cadd score2: $edivavals[23] \n";
			print "Simple tandem repeat region: $edivapublicvals[0] \n";			
			print "Simple tandem repeat length: $edivapublicvals[1] \n";
		}	## end while
					
	} ## end of deciding if-else for -q paramter 		
}