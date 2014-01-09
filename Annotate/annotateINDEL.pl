#!/usr/bin/perl
use warnings;
use strict;
use Digest::MD5 qw(md5);
use DBI;
use threads;
use threads::shared;
use DateTime;
use Getopt::Long;



sub usage { print "\n$0 \n usage:\n",
	   "--input \t provide the VCF file to annoate\n",
	   "--tempDir \t provide temporary location (all temp files with the current session will be removed at the end of execution) \n",
	   "--geneDefinition \t the gene definition you want to select for Annovar (ensGene,refGene,knownGene) [default: all] \n",
	   "--help \t show help \n",
	   "\n\n perl $0 --input=/path_to/inout.vcf --tempDir=/path_to_temp",
	   "\n For more details, send an email to : rubayte.rahman\@crg.eu\n";
}

## variables
my $help = 0;
my $input; ## main input file
my $geneDef = "all"; ## geneDefinition
my $sep = ",";
my %indels = (); ## hash to hold input INDELs from VCF
my @thrds = (); ## list to hold threads
our %eDiVa :shared = (); ## hash to hold input INDELs with eDiVa annotation
our %SIFT :shared = (); ## hash to hold input INDELs with SIFT annotation
our %Annovar :shared = (); ## hash to hold input INDELs with ANNOVAR annotation
our %AnnovarRes :shared = (); ## hash to hold input INDELs with ANNOVAR annotation

## temp location
our $templocation;##= "/users/so/rrahman/IndelAnnotationScratch";
our $fileSuffix = DateTime->now->ymd;
$fileSuffix = $fileSuffix."-".DateTime->now->hms;
$fileSuffix =~ s/\:/-/g;

## SIFT settings
our $SIFTDirectory = "/users/GD/tools/sift5.0.3";
our $siftCodingInfoGffFile = "ens.hum.ncbi37.ver66.cds.merge.gff";

## ANNOVAR settings
our $ANNOVAR = "/users/so/rrahman/Annovar_INDEL/annovar_2012Oct23"; 

## grab command line options
GetOptions("input=s" => \$input, "tempDir=s" => \$templocation, "geneDefinition=s" => \$geneDef, "help=s" => \$help);

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

## subroutine for converting VCF indel to SIFT format
sub prepareSIFTformat
{
	open (SIFTIN, ">>$templocation"."/$fileSuffix"."/siftInfile".$fileSuffix) or die "Cant open new file in $templocation \n";
	#open (SIFTIN, ">>$templocation"."/siftInfile".$fileSuffix) or die "Cant open new file in $templocation \n";
	while( my ($k, $v) = each %indels ) 
	{ 
		## variant format :: $chr;$position;$ref;$alt
		my @dt = split(/\;/,$v);
	    if (length($dt[2])> length($dt[3])) ## deletion
	    {
	            my $endpos = $dt[1]+length($dt[2]) - 1;
	            my $refs = substr($dt[2],1); ## remove the initial nucleotide, only deleted alteration
	            #$SIFT { $v } = "$dt[0],$dt[1],$endpos,1,$refs/";
	            print SIFTIN $dt[0].",".$dt[1].",".$endpos.",1,".$refs."\/,$v\n";
	    }else{## insertion
	            my $alts = substr($dt[3],1); ## remove the initial nucleotide, only inserted alteration
	           	#$SIFT { $v } = "$dt[0],$dt[1],$dt[1],1,$alts";
	            print SIFTIN $dt[0].",".$dt[1].",".$dt[1].",1,".$alts.",$v\n";
	    }
	}
	close(SIFTIN);
}

## subroutine for converting SIFT input to IndelPrediction input
sub reformat_chrfile_032811
{
	my $fileTOreformat = shift;

	open (REF, "$templocation"."/siftInfile".$fileSuffix) or die "Cant open $templocation"."/siftInfile".$fileSuffix." file \n";
	open (ORGFILE, ">>$templocation"."/siftInfile".$fileSuffix."_indelPred") or die "Cant open new file to write IndelPrediction input \n";
	while(<REF>)
	{
		chomp $_;

		$_ =~ s/^\s+|\s+$//g;
        my $chr = "";
        my $start = "";
        my $stop = "";
        my $orn = "";
        my $allele = "";
        my $comment = "";
        my @elts = split /,/, $_;
        $chr = $elts[0];
        $start = $elts[1];
        $stop = $elts[2];
        $orn = $elts[3];
        $allele = $elts[4];
        $comment = $elts[5];

        if (!defined($comment)) { $comment = ""; } # Comment can be empty if user does not have one

        next if ($chr eq "");
        #rare bad case if stop < start
        if ($start > $stop){
                next;
        }       

        #if insertion
        if ($start == $stop){
                if ($allele =~ /(\w+)/){        #if already in snp_classifier format (no slash): ATGGC
                        $allele = uc ($1);
                }
                elsif($allele =~ /^\-*\/(\w+)/){  #if of the form -/ATTGCA -> ATTGCA
                        $allele = uc ($1)
                }
                elsif ($allele =~ /^(\w+)\/\-*/){  #if of the form ATTG/- -> ATTG
                        $allele = uc ($1);
                }
        }

         #if deletion
        if ($stop - $start >= 1){
                $allele = "-/";
        }
        if ($comment !~ /\w|\d/){
                #print TEMPFILE "$chr,$start,$stop,$orn,$allele\n";
                 if ($start == $stop) {
                        print ORGFILE "$chr,$start,$stop,$orn,$allele\n";
                        }
                else{
                        print ORGFILE "$chr,$start,$stop,$orn,/\n";
                        }

        }
        else{
               #print TEMPFILE "$chr,$start,$stop,$orn,$allele,$comment\n";
               if ($start == $stop) {
                        print ORGFILE "$chr,$start,$stop,$orn,$allele,$comment\n";
               }
                else{
                        print ORGFILE "$chr,$start,$stop,$orn,/,$comment\n";
                }

        }

	} ## end While

	close(ORGFILE);
	close(REF);
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
	my $dbh = DBI->connect('dbi:mysql:'.$database.';host=localhost',$username,'') or die "Connection Error!!\n";

	# extract db result
	while( my ($k, $v) = each %indels ) 
	{
		my ($chr,$pos,$ref,$alt) = split(/\;/,$v);
		$sql = "select annotateINDEL('$k',$chr,$pos);";
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

## subroutine for SIFT annotation
sub SIFTAnnotationNew
{
	my $crashLoc = $templocation."/$fileSuffix";
	my $mkdirCmm = "mkdir ".$crashLoc;
	system ($mkdirCmm);
	## prepare sift input 
	&prepareSIFTformat;
	## sift input file
	my $FileIndelPred = $templocation."/$fileSuffix/siftInfile".$fileSuffix;
	##start SIFT with the created sift input
	my $siftCmm = "perl $SIFTDirectory/bin/SIFT_exome_indels_corrected.pl -i $FileIndelPred -c $SIFTDirectory/coding_info/Coding_info_37 -d /no_backup/GD/projects/eDiVa/eDiVaSupport/sift_new_institute_sg/May31_2013/data -o $crashLoc";	
	system($siftCmm);
	my $findCmm = "find $crashLoc/ -name '*condensed_classification_file' -exec scp {} $crashLoc/condensed_classification_file$fileSuffix \\;";
	system ($findCmm);
	open (SIFT, "$crashLoc/condensed_classification_file$fileSuffix") or die "Cant open $crashLoc\/condensed_classification_file$fileSuffix file \n";

	while(<SIFT>)
	{
		chomp $_;
		my @dt = split(/\t/,$_);
		my @ids = split(/\,/,$dt[0]);
		my $newindelformat = $ids[5];

		my $confidence = $dt[27];
		my $effect = $dt[26];
		my $score;
		my $war;
		if ($_ =~ m/NCBI/)
		{
			$war = "Yes";
		}else{
			$war = "No";
		}

		if ($confidence and $effect)
		{
			$score = "$confidence($effect)";
		}
		else{
			$score = "NA";
		}
		if (! $dt[7])
		{
			$dt[7] = "NA";
		}
		if (! $dt[13])
		{
			$dt[13] = "NA";
		}
		my $valToPass = "$dt[7],$dt[13],$score,$war";
		$SIFT { $newindelformat } = $valToPass;
	}

	close(SIFT);

}

## subroutine for SIFT annotation
sub SIFTAnnotation
{	
	## prepare sift input 
	&prepareSIFTformat;
	## prepare sift IndelPrediction input
	&reformat_chrfile_032811;
	## NOTE TO MYSELF :: In future this two conversions will be merged to one
	my $createDir = "mkdir ".$templocation."/".$fileSuffix."";
	my $filetoPredict = "$templocation"."/siftInfile".$fileSuffix."_indelPred";
	open (LOAD, $filetoPredict) or die "Cant open $filetoPredict file \n";
	while (<LOAD>)
	{
		chomp $_;
		my @dt = split(/\,/,$_);
		$dt[5] =~ s/^\s*(.*?)\s*$/$1/;
		$SIFT { $dt[5] } = "$dt[0],$dt[1],$dt[2],$dt[3],$dt[4]";
	}
	close(LOAD);
	my $tmpIndelPred = $templocation."/".$fileSuffix."";
	my $resPrediction = $filetoPredict."_res";
	system ($createDir);

	system("python $SIFTDirectory/bin/non3nindels/indelPredPipeline.py $filetoPredict $tmpIndelPred $SIFTDirectory/coding_info/Coding_info_37");

    open( PREDFILE, "$filetoPredict.predictions" ) or die "Cant open filetoPredict.predictions file \n";
    open (RES, ">>$resPrediction") or die "Cant open new file $resPrediction \n";
    while (<PREDFILE>)
    {
		chomp;
		my @predictions = split /\s+/, $_;
		my $indel = $predictions[0];

		my @indelParts = split /,/, $indel;
		my $newindelformat;
		if ( $indelParts[3] eq "/"){
		    #$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,-".$indelParts[3];
		   	$newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[3];
		}
		else{
		    $newindelformat = $indelParts[0].",".$indelParts[1].",".$indelParts[2].",1,".$indelParts[3];
		}

		#print $newindelformat."\n";
		my $affectedEnsgID = $predictions[1];
		my $effect  = $predictions[2];
		my $confidence = $predictions[3];
		my $classificationRule = $predictions[4];
		my @array = ($effect, $confidence, $classificationRule);
		#@{$indelClassificationResults{$newindelformat}} = join(';',@array);
		#print $newindelformat."\t".$effect."\n";
		print RES $newindelformat.",",join (';', @array)."\n";
		while ( my ($key, $value) = each(%SIFT) ) 
		{
        	if ($value eq $newindelformat)
        	{
        		$SIFT { $key } = "$confidence($effect)";
        	}
		}
    }
    close(RES);
    close(PREDFILE);
	## we dont care about the other two files, cause they dont have any SIFT score
}

## subroutine for Annovar annotation
sub AnnovarAnnotation
{	
	my $annCmm;
	## prepare Annovar input
	my $annInCmm = "perl $ANNOVAR/convert2annovar.pl --includeinfo -format vcf4 $input > $templocation/annInfile".$fileSuffix."   2> ".$input.".annovar.log";
	print "MESSAGE :: Running Annovar command \> $annInCmm\n";
	system ($annInCmm);

	my $annFile = "$templocation/annInfile".$fileSuffix."";
	## load annovar hash
	open (LOAD, $annFile) or die "Cant open $annFile file \n";
	while (<LOAD>)
	{
		chomp $_;
		my @dt = split(/\t/,$_);
		if ($dt[9] =~ m/\,/)
		{
			my @alss = split(/\,/,$dt[9]);
			foreach my $alsss (@alss)
			{
				my $newKey = $dt[5].";".$dt[6].";".$dt[8].";".$alsss;
				my $newVal = $dt[0].";".$dt[1].";".$dt[2].";".$dt[3].";".$dt[4];
				$Annovar { $newKey } = $newVal;
			}
		}else{
			my $newKey = $dt[5].";".$dt[6].";".$dt[8].";".$dt[9];
			my $newVal = $dt[0].";".$dt[1].";".$dt[2].";".$dt[3].";".$dt[4];
			$Annovar { $newKey } = $newVal;
		}
	
	}
	close(LOAD);

	## run Annovar annotation
	if ($geneDef eq 'ensGene')
	{
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);		
	}elsif($geneDef eq 'refGene')
	{
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}elsif($geneDef eq 'knownGene')
	{
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
	}else{
		print "MESSAGE :: No sepicific gene definition selected, hence Annovar is going to run on all definitions !\n";
		## refgene
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --step 1 --outfile $templocation/Refseq$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## ensgene
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --genetype ensgene --step 1 --outfile $templocation/Ensembl$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
		## knowngene
		$annCmm = "$ANNOVAR/eDiVa_indels_summarize_annovar.pl --buildver hg19  --genetype knowngene --step 1 --outfile $templocation/Known$fileSuffix $templocation/annInfile".$fileSuffix." $ANNOVAR/hg19/ 2> ".$input.".annovar.log";
		print "MESSAGE :: Running Annovar command \> $annCmm\n";
		system ($annCmm);
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
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[23].";".$dt[24].";".$dt[25];
				##print $valueTOmatch."\n";
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];
				while ( my ($key, $value) = each(%Annovar) ) 
				{
		        	if ($value eq $valueTOmatch)
		        	{
		        		$AnnovarRes { $key } = $annToPass;
		        	}
				}
			}
		}		
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
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[23].";".$dt[24].";".$dt[25];
				##print $valueTOmatch."\n";
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];
				while ( my ($key, $value) = each(%Annovar) ) 
				{
		        	if ($value eq $valueTOmatch)
		        	{
		        		if ($geneDef eq "all")
		        		{
		        			$AnnovarRes { $key } = $AnnovarRes { $key }.",".$annToPass;
		        		}else
		        		{
		        			$AnnovarRes { $key } = $annToPass;
		        		}
		        	}
				}
			}
		}		
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
				my $valueTOmatch = $dt[21].";".$dt[22].";".$dt[23].";".$dt[24].";".$dt[25];
				##print $valueTOmatch."\n";
				my $annToPass = $dt[0].",".$dt[1].",".$dt[2].",".$dt[3];
				while ( my ($key, $value) = each(%Annovar) ) 
				{
		        	if ($value eq $valueTOmatch)
		        	{
		        		if ($geneDef eq "all")
		        		{
		        			$AnnovarRes { $key } = $AnnovarRes { $key }.",".$annToPass;
		        		}else
		        		{
		        			$AnnovarRes { $key } = $annToPass;
		        		}
		        	}
				}
			}
		}		
	}
	
	close(ANNE);
	close(ANNR);
	close(ANNK);
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
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SugMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTAminoAcidChange,NMD,SIFTscore,NCBIWarning";
	}elsif($geneDef eq 'refGene')
	{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),
		dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SugMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTAminoAcidChange,NMD,SIFTscore,NCBIWarning";
	}elsif($geneDef eq 'knownGene')
	{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,
		TotalEVSFrequecy,Eur1000GenomesFrequency,Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SugMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTAminoAcidChange,NMD,SIFTscore,NCBIWarning";
	}else{
		$stringTOreturn = "Chr,Position,Reference,Alteration,Function(Refseq),Gene(Refseq),ExonicFunction(Refseq),AminoAcidChange(Refseq),Function(Ensembl),Gene(Ensembl),ExonicFunction(Ensembl),
		AminoAcidChange(Ensembl),Function(Known),Gene(Known),ExonicFunction(Known),AminoAcidChange(Known),dbsnpIdentifier,dbSNPfrequency,EurEVSFrequecy,AfrEVSFrequecy,TotalEVSFrequecy,Eur1000GenomesFrequency,
		Afr1000GenomesFrequency,Amr1000GenomesFrequency,Asia1000GenomesFrequency,Total1000GenomesFrequency,SugMentDup,PlacentalMammalPhyloP,PrimatesPhyloP,VertebratesPhyloP,PlacentalMammalPhastCons,
		PrimatesPhastCons,VertebratesPhastCons,Score1GERP++,Score2GERP++,SIFTAminoAcidChange,NMD,SIFTscore,NCBIWarning";
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
	if ($_ !~ m/^#/) ## skip header line
	{
		my @line = split(/\t/,$_);
		my $chr = $line[0];
		my $position = $line[1];
		my $ref = $line[3];
		my $alt = $line[4];

		## take care of chr1 or Chr1 and convert to chr1 /Chr1-> 1
		if ($chr =~ m/^chr/ or $chr =~ m/^Chr/)
		{
			$chr = substr($chr,3);
		}
		## process based on altearation
		if ($alt =~ m/\,/)
		{
			my @alts = split(/\,/,$alt);
			foreach my $al(@alts)
			{
				## make indelID
				my $token_ref = unpack('L', md5($ref));
				my $token_obs = unpack('L', md5($al));
				my $indel = $chr.';'.$position.';'.$token_ref.';'.$token_obs;
				$indels{ $indel } = "$chr;$position;$ref;$al";
			}
		}else{
			## make indelID
			my $token_ref = unpack('L', md5($ref));
			my $token_obs = unpack('L', md5($alt));
			my $indel = $chr.';'.$position.';'.$token_ref.';'.$token_obs;
			$indels{ $indel } = "$chr;$position;$ref;$alt";
		}

	} ## end if
} ## end while
close(INPUT);
## Initialization completed
print "\nMESSAGE :: Intialization completed .. ";

## start threading and annotating
print "\nMESSAGE :: Annotation starting .. This may take a few moments ..";

## spawn threds and push to list
push @thrds, threads -> new(\&eDiVaAnnotation); ## spawn a thread for eDiVa annotation
push @thrds, threads -> new(\&SIFTAnnotationNew); ## spawn a thread for SIFT annotation
push @thrds, threads -> new(\&AnnovarAnnotation); ## spawn a thread for Annovar annotation

## join spawned threads
foreach my $thr (@thrds)
{
        $thr -> join();
}

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
while (my($key, $value) = each(%indels)) 
{
	my $newKey = $value;
	my ($edivaannotationtoprint,$annovarannotationtoprint) = ("NA","NA");
	my ($chr,$position,$ref,$alt) = split(/\;/, $value);
	$edivaannotationtoprint = $eDiVa{$key} if $eDiVa{$key};
	my $siftannotationtoprint = $SIFT{$newKey};

	$annovarannotationtoprint = $AnnovarRes{$newKey} if $AnnovarRes{$newKey};

	## clean SIFT annotation 
	if (! $siftannotationtoprint)
	{
		$siftannotationtoprint = "NA,NA,NA,NA";
	}
	## clean annovar annotation
	$annovarannotationtoprint =~ s/\"//g;
	$annovarannotationtoprint =~ s/,,/,NA,/g;
	$annovarannotationtoprint =~ s/NA,/NA,NA/g;
	## write annotation to file 
	print ANN $chr.$sep.$position.$sep.$ref.$sep.$alt.$sep.$annovarannotationtoprint.$sep.$edivaannotationtoprint.$sep.$siftannotationtoprint."\n";
}


close(ANN);

my $srtCmm = "sort -n -k1,1 -n -k2,2 --field-separator=, $outFile > $outFile.sorted";
system($srtCmm);

## writing completed
print "\nMESSAGE :: Writing annotation completed .. ";
print "\nMESSAGE :: Your annotated files are $outFile \& $outFile.sorted";

## Finalize everything
print "\nMESSAGE :: Finalizing annotation process ..";
&finalize;
print "\nMESSAGE :: Finalization completed !\n";
