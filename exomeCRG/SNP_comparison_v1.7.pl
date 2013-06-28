#!/usr/bin/perl

#######################################################################
# SNP_comparison.pl --infile/--listinfile file --outfile [--delimiter --frequencycutoff --genecutoff --MAFcutoff --FUNCfilter --EXONICFUNCfilter] --help


### will be changed to Indel_comparison


# --infile may be invoked several times
# --outfile specifies the name of the result file
# --delimiter may be set if not tab, must be a unique letter  in your file. [default 'tab']
# --cutoff report only SNPs that are found in at least as many files as cutoff [default = 0]
# --nodescription print no description of each SNP once?
# --help show usage
# compares two or more SNP lists and prints out SNPs that are found in 
# multiple lists
# different comparisons can be thought of:
# - which files contain the same SNP hit
# - which files contain the same Gene hit
# -- collate the function of the hit genes
# VERSION
# 1.6
# include exclusionlist from ftp://ftp.nhgri.nih.gov/pub/NIHUDP/ADAMS_METHODS/Exclusion_lists_and_linkage_file/Gene%20exclusion%20list.txt (tr '\015\032' '\n' < gene_exclusion_list.txt > gene_exclusion_list.txt_unix)
# 1.5
# change the way SNP_comparison handles files completely
# files should be unpacked now all and all files should be read line by line in parallel
# ranking improved a lot
# 1.44
# check input read from files better
# 1.43
# calculate allelfrequencies
# 1.42
# use threads while going through all reference, coverage & allSNPfiles
#  threading is built with the help of http://chicken.genouest.org/perl/multi-threading-with-perl/
# enter data from EVS Server
# 1.31 & 1.311 
# try to access result files directly
# 1.3
# investigation of additional fields integrated
# genecutoff, MAFcutoff, funcfilter, exonicfuncfilter
# additionally whole rewrite of the filters for SNP/GENE compilation
# the output is now based on the frequency of GENEs in samples not SNPs
# introduced ranking
# 1.2
# nodescription
# one line one SNP position containing all file names and data
# listinfile digests a file containing a list of paths to files
# 1.1
# cutoff
# one line of data for all SNP entries
# 1.0
# first stable release comparing multiple files
#

###TODO
#fits family model (SNP_apply family model)
#
#exclude IGSF3 DONE FAM104B DONE CACNA1B (9:140777306) DONE (whole gene)
### end TODO


#######################################################################

use warnings;
use strict;
#use threads;
# use IPC::Shareable;
use Getopt::Long; # read command line arguments easily
use File::Spec; # split path names into location/path/filename
use File::Path; # find out file locations
use Data::Dumper; # facilitate debugging
use IO::File; # use object oriented way to store FileHandles (is not possible in Hashes with standardFH)
use List::Util; #easily find min and max of a list of numbers


#============== initialize =================

my $EVSfile = '/users/GD/resource/human/hg19/databases/EVS/ESP5400.chrall.snps.vcf';
my $exclusionList = '/users/GD/resource/human/gene_exclusion_list/gene_exclusion_list.txt';
my $condelcall = '/users/so/odrechsel/scripts/condel/bin/was-modified.pl /users/so/odrechsel/scripts/condel/config/';
my $mutationAssessor = '/users/GD/resource/human/hg19/databases/mutationAssessor/mutationassessor_exonicregions_120425';
my $siftcall = '/users/GD/tools/sift/SIFTindel.sh';

my $dirlistinfile = '';
my $outfile = '';
my $collectionfile = '';
my $mailaddress = '';

my %allGENEdata;
my $delimiter = ",";
my $nodescription = 0;
my $help = 0;
my $useexclusionlist = 0;
my %samplefiles;
my %excludedgenes;
my $allyes = 0;
my $unlink = 0;
my $debug = 0;
my @errormessages;

#===filter variables===
my $frequencycutoff = 0;
my $genecutoff = 0;
my $MAFcutoff = 1;
my $funcfilter = '';
my $exonicfuncfilter = '';
use constant REFERENCE_HOM => 0.7;
use constant REFERENCE_HET => 0.3;

#============== read parameters ============
GetOptions("dirlistinfile=s" => \$dirlistinfile, "outfile=s" => \$outfile, "collectionfile=s" => \$collectionfile, "mailaddress=s" => \$mailaddress,
	   "evsfile=s" => \$EVSfile, "useexclusionlist" => \$useexclusionlist, "exclusionlist=s" => \$exclusionList, "mutationassessor=s" => \$mutationAssessor,
	   "help" => \$help, # "nodescription" => \$nodescription,
	   #"genecutoff=i" => \$genecutoff, "MAFcutoff=f" => \$MAFcutoff, "FUNCfilter=s" => \$funcfilter, "EXONICFUNCfilter=s" => \$exonicfuncfilter, 
	   "allyes" => \$allyes, "unlink" => \$unlink, "debug" => \$debug); # GetOptions("command line parameter=string" =>(saved in) \variable defined earlier)
	   #unlink and allyes are only for debugging use

#============== check parameters ===========
if (( ($outfile eq '' || $dirlistinfile eq '') && $collectionfile eq '' && $mailaddress eq '') || $help == 1 ) { # belch the usage line and terminate if not all essential parameters are given.  
	&usage;
}
else {
	$samplefiles{'allsamples'}{'mailaddress'} = $mailaddress;
}
#============== start work =================
if ($collectionfile eq '') {
	#============== read files =================
	if ($dirlistinfile) {
		open (DIRLISTINFILE, $dirlistinfile) or die "could not open $dirlistinfile \n";
		while(<DIRLISTINFILE>) {
			chomp($_);
			my $directory = $_;
			$directory =~ /(.+)\/(\S+?)$/;
			my $samplename = $2;
			my $snpfile = $directory."/SNP_Intersection/AnnovarIntersection/sum.exome_summary.csv";
			$snpfile = checkannovarversion($snpfile); #path could be either AnnovarIntersection or Annovar
			my $coveragefile = $directory."/shore/CoverageAnalysis/coverage.0.gff.gz";
			my $SNPintersection = $directory."/SNP_Intersection/merged.all.vcf";
			my $referencecalls = $directory."/shore/Variants/ConsensusAnalysis/reference.shore.gz";
			
			my $indelfile = $directory."/Indel_Intersection/AnnovarIntersection/sum.exome_summary.csv";
			my $indelintersection = $directory."/Indel_Intersection/merged.all.vcf";
			
			
			next unless (&filetest($snpfile, \@errormessages));
			next unless (&filetest($coveragefile, \@errormessages));
			next unless (&filetest($SNPintersection, \@errormessages));
			next unless (&filetest($referencecalls, \@errormessages));
			next unless (&filetest($indelfile, \@errormessages));
			next unless (&filetest($indelintersection, \@errormessages));
			
			
			$samplefiles{$samplename}{'SNP'}{'name'} = $snpfile;
			$samplefiles{$samplename}{'coverage'}{'name'} = $coveragefile;
			$samplefiles{$samplename}{'allSNP'}{'name'} = $SNPintersection;
			$samplefiles{$samplename}{'reference'}{'name'} = $referencecalls;
			$samplefiles{$samplename}{'Indel'}{'name'} = $indelfile;
			$samplefiles{$samplename}{'allIndel'}{'name'} = $indelintersection;
			
			print "found all necessary files for $samplename\n" # MESSAGE
		}
		close(DIRLISTINFILE);
	}
	
	# unzip files
	foreach my $sample (keys %samplefiles) {
		next if ($sample eq 'allsamples');
		$samplefiles{$sample}{'coverage'}{'unzipped'} = unzipfile($sample, 'coverage', \%samplefiles, $allyes);
		$samplefiles{$sample}{'reference'}{'unzipped'} = unzipfile($sample, 'reference', \%samplefiles, $allyes);
	}
	
	# open all file handles
	my $evsfilefound = 0;
	my $mutassfilefound = 0;
	openAllFiles(\%samplefiles, $outfile, $EVSfile, \$evsfilefound, $mutationAssessor, \$mutassfilefound);
	
	# print header line in outputfile
	writeHeader(\%samplefiles, $delimiter, 'collectionfilehandle'); #REAL
	
	# run through all files to obtain all information
	my $flag = 1;
	
	# initialize the SNP/Indel information
	fillLastSNPline(\%samplefiles);
	$samplefiles{'allsamples'}{'currentposition'} = identifyCurrentPosition(\%samplefiles, $delimiter);
	
	#print Dumper(%samplefiles); #DEBUG
	#exit; #DEBUG
	
	# run through SNPs
	print "compiling information \n";
	while ($flag) {
		my (@positiveSamples, @negativeSamples);
		#identify which samples show a high quality SNP/Indel at the currentposition and which don't 
		identifySamples(\%samplefiles, $delimiter, \@positiveSamples, \@negativeSamples);
		# check if all files have reached eof
		unless (nothingToDo(\%samplefiles)) { 
			$flag = 0;
		}
		#enterDataForPositveSamples
		enterDataForPositveSamples(\%samplefiles, \@positiveSamples, \%allGENEdata, $delimiter);
		
		#print "enterDataForNegativeSamples\n";
		enterDataForNegativeSamples(\%samplefiles, \@negativeSamples, $delimiter);
		
		#identifyQualityMedian
		identifyQualityMedian(\%samplefiles, $delimiter);
		calculateAverageCoverage(\%samplefiles);
		
		# enter EVS values or NA values if no values seen
		# enterEVSdata(\%samplefiles, $evsfilefound);
		
		# enter Mutation Assessor values or NA values if no values seen
		# enterMutationAssessor(\%samplefiles, $mutassfilefound, $delimiter);
		
		# determine allel frequencies (high & low qual SNPs GATK, MPILEUP, SHORE - so no need to look into het/hom from Annovar)
		determineAllelFrequencies(\%samplefiles, $funcfilter, $exonicfuncfilter); # count allel and gene numbers but exclude genes that get filtered
	
		# write results
		
		writeResults(\%samplefiles, $delimiter, $frequencycutoff, $genecutoff, $MAFcutoff, $funcfilter, $exonicfuncfilter);
		$samplefiles{'allsamples'}{'currentposition'} = identifyNextPosition(\%samplefiles, $delimiter, \$flag, \@positiveSamples);
	}
	
	# close file handles
	closeAllFiles(\%samplefiles);

} # end if collectionfile is given
else {
	$collectionfile =~ /(.+?)\.collection/;
	$outfile = $1;
}


# re-read the output
print "re-read the output\n";
my %alldata;
openOutFile (\%samplefiles, $outfile, $exclusionList, $useexclusionlist);
readInitialOutput(\%samplefiles, \%alldata, $delimiter, \%allGENEdata, $genecutoff, $MAFcutoff, $funcfilter, $exonicfuncfilter);

# do condel
#print "do condel\n";
#doCondel(\%samplefiles, $condelcall);
#readCondelOutput(\%samplefiles, \%alldata);

# TODO: check usefulness of this section!
# do ranking
#print "do ranking\n";
#groupvalues(\%alldata);
#rank('MAF_bin', 'ascending', 'MAF_Rank', \%alldata);
#rank('SegDup_bin', 'ascending', 'SegDup_Rank', \%alldata);
##rank('Condel_bin', 'descending', 'Condel_Rank', \%alldata);
#rank('Conserved_bin', 'ascending', 'Conserved_Rank', \%alldata);
#rank('LJB_PhyloP_bin', 'descending', 'LJB_PhyloP_Rank', \%alldata);

#rankproduct(\%alldata);

#rank('RankProduct', 'ascending', 'RankProduct_Sorted', \%alldata);

## apply exclusion list
#print "read exclusion list\n";
#if ($useexclusionlist == 1) {
#	readExclusionList(\%samplefiles, \%excludedgenes);
#}

writeSIFTinput(\%samplefiles, \%alldata);
writeSIFTconfig(\%samplefiles, $outfile);
doSIFT(\%samplefiles, $siftcall); ### DEBUG - REAL
readSIFT(\%samplefiles, \%alldata); ### DEBUG - REAL

writeFinalResults(\%samplefiles, \%alldata, $delimiter, $debug, \%allGENEdata, $genecutoff, \%excludedgenes, $useexclusionlist);

# delete unzipped files
print "delete unzipped files\n";
if ($collectionfile eq '') {
	deleteUnzipped(\%samplefiles); # REAL
	deleteSIFT(\%samplefiles)
}
else {
	deleteCondel(\%samplefiles); # REAL
	deleteSIFT(\%samplefiles)
}

###### subroutines ######

sub calculateAverageCoverage{
	my $samplefiles = shift;
	
	my $coverage_sum = 0;
	foreach my $sample (keys %{$samplefiles}) {
		next if $sample eq 'allsamples';
		
		if ($$samplefiles{$sample}{'coverage'}{'value'} ne 'NA') {
			$coverage_sum += $$samplefiles{$sample}{'coverage'}{'value'};
		}
	}
	
	my $coverage_average = $coverage_sum / (scalar(keys %{$samplefiles}) - 1);
	$$samplefiles{'allsamples'}{'coverage_average'} = sprintf("%.2f",$coverage_average); 
	
};
sub checkannovarversion{
	my $snpfile = shift;
	if(-e $snpfile) {
		return($snpfile);
	}
	$snpfile =~ s/AnnovarIntersection/Annovar/;
	if(-e $snpfile) {
		return($snpfile);
	}
	else {
		filenotfoundexception($snpfile);
	}
}

sub checkLowQualSNPs {
	my $samplefiles = shift;
	my $sample = shift;
	my $currentposition = shift;
	
	my $flag = 1;
	while ($flag) {
		my $line;
		
		# check if a line was ever read and fill it if not
		if (exists $$samplefiles{$sample}{'allIndel'}{'line'} ) {
			$line = $$samplefiles{$sample}{'allIndel'}{'line'};
		}
		else {
			$line = $$samplefiles{$sample}{'allIndel'}{'filehandle'} -> getline();
			$$samplefiles{$sample}{'allIndel'}{'line'} = $line;
		}
		
		if ($$samplefiles{$sample}{'allIndel'}{'filehandle'} -> eof()) { # test for EOF # was $line
			$$samplefiles{$sample}{'allIndel'}{'GATK'} = 'NA';
			$$samplefiles{$sample}{'allIndel'}{'MPILEUP'} = 'NA';
			$$samplefiles{$sample}{'allIndel'}{'SHORE'} = 'NA';
			$$samplefiles{$sample}{'allIndel'}{'foundposition'} = 'NA';
			
			$flag = 0;
			last;
		}
		
		unless ($$samplefiles{$sample}{'allIndel'}{'line'} || $$samplefiles{$sample}{'allIndel'}{'line'} == '') {
			print "$sample , leads to missing line $$samplefiles{$sample}{'allIndel'}{'line'} \n";
		}
		
		# check if the line is a header line, if so read next line and move on
		if ($line =~ /^#/) {
			$line = $$samplefiles{$sample}{'allIndel'}{'filehandle'} -> getline();
			$$samplefiles{$sample}{'allIndel'}{'line'} = $line;
			next;
		}
		
		else {
			my @data = split("\t", $line);
			unless (scalar (@data) == 12) { #check for a correct input line
				print "merged.all.vcf of $sample does not contain correct number of elements in line\n", $line, "\n";
				next;
			}
			
			#ANNOVAR changes start positions "1       2938407 del_1_2938407   CAGA    C" gets transformed to "1,2938408,2938410,AGA,-"
			if (length($data[3]) > length($data[4])) { # ref is longer than observation -> deletion
				$data[1] += length($data[4]); # increment by length of observed nucleotide stretch, which would be done by ANNOVAR 
			}
			
			
			my $position = encryptposition($data[0], $data[1]);
			
			#locate the correct position
			if ($position > $currentposition) {
				#print 'oops missed the correct position (', $currentposition, ') in merged.all.vcf line of ', $sample, "\n", $line, "\n"; # DEBUG
				$$samplefiles{$sample}{'allIndel'}{'GATK'} = 'NA';
				$$samplefiles{$sample}{'allIndel'}{'MPILEUP'} = 'NA';
				$$samplefiles{$sample}{'allIndel'}{'SHORE'} = 'NA';
				$$samplefiles{$sample}{'allIndel'}{'foundposition'} = $position;
				$flag = 0;
				last;
			}
			elsif ($position == $currentposition) {
				#unless ($$samplefiles{'allSNP'}{'line'}) {print "position in non defined line is $position \n";}
				my @allsnpdata = split(/\t/, $line); # $$samplefiles{'allSNP'}{'line'}
				unless (missingvalue($data[9], $sample, 'merged.all.vcf')) {
					$data[9] = './.'; # GATK
				};
				unless (missingvalue($data[10], $sample, 'merged.all.vcf')) {
					$data[10] = './.'; # MPILEUP
				};
				unless (missingvalue($data[11], $sample, 'merged.all.vcf')) {
					$data[11] = './.'; # Shore
				};
				$$samplefiles{$sample}{'allIndel'}{'GATK'} = $data[9];
				$$samplefiles{$sample}{'allIndel'}{'MPILEUP'} = $data[10];
				$$samplefiles{$sample}{'allIndel'}{'SHORE'} = $data[11];
				$$samplefiles{$sample}{'allIndel'}{'foundposition'} = $position;
				$flag = 0;
				last;
			}
		}
		# get next line
		$$samplefiles{$sample}{'allIndel'}{'line'} = $$samplefiles{$sample}{'allIndel'}{'filehandle'} -> getline();
	}
}

sub checkReference {
	my $samplefiles = shift;
	my $sample = shift;
	my $currentposition = shift;

	my $flag = 1;
	while ($flag) {
		my $line;
		if (exists $$samplefiles{$sample}{'reference'}{'line'} ) {
			$line = $$samplefiles{$sample}{'reference'}{'line'};
		}
		else {
			$line = $$samplefiles{$sample}{'reference'}{'filehandle'} -> getline();
			$$samplefiles{$sample}{'reference'}{'line'} = $line;
		}
		
		if ($$samplefiles{$sample}{'reference'}{'filehandle'} -> eof()) { # test for EOF # was $line
			$$samplefiles{$sample}{'reference'}{'nucleotide'} = 'NA';
			$$samplefiles{$sample}{'reference'}{'foundposition'} = 'NA';
			$flag = 0;
			last;
		}
		
		my @data = split("\t", $line);
		unless (scalar (@data) == 10) {
			print "reference.shore.gz of $sample does not contain correct number of elements in line\n", $line, "\n";
			next;
		}
		
		my $position = encryptposition($data[1], $data[2]);
		
		if ($position < $currentposition) {
			$line = $$samplefiles{$sample}{'reference'}{'filehandle'} -> getline();
			$$samplefiles{$sample}{'reference'}{'line'} = $line;
			next;
		}
		elsif ($position > $currentposition) {
			#print 'oops missed the correct position (', $currentposition ,') in reference.shore.gz of ', $sample, ' in line', "\n", $line, "\n"; # DEBUG
			$$samplefiles{$sample}{'reference'}{'nucleotide'} = 'NA';
			$$samplefiles{$sample}{'reference'}{'foundposition'} = $position;
			$$samplefiles{$sample}{'reference'}{'state'} = 'NA';
			$flag = 0;
			last;
		}
		elsif ($position == $currentposition){
			if (($data[3] eq $data[4]) && ($data[7] >= REFERENCE_HOM)) { # $data[7] > 0.9 checks for homozygous reference call, if not it's not entered - decreased to 0.7
				$$samplefiles{$sample}{'reference'}{'nucleotide'} = $data[3];
				$$samplefiles{$sample}{'reference'}{'foundposition'} = $position;
				$$samplefiles{$sample}{'reference'}{'state'} = 'hom';
			}
			elsif (($data[3] eq $data[4]) && (REFERENCE_HET <= $data[7] && $data[7] < REFERENCE_HOM)) {
				$$samplefiles{$sample}{'reference'}{'nucleotide'} = $data[3];
                                $$samplefiles{$sample}{'reference'}{'foundposition'} = $position;
                                $$samplefiles{$sample}{'reference'}{'state'} = 'het';
			}
			else {
				$$samplefiles{$sample}{'reference'}{'nucleotide'} = 'N';
				$$samplefiles{$sample}{'reference'}{'foundposition'} = $position;
			}
			$flag = 0;
			
			last;
		}
		
		
	}
}

sub checkCoverage {
	my $samplefiles = shift;
	my $sample = shift;
	my $currentposition = shift;
	
	my $flag = 1;
	while ($flag) {
		my $line;
		
		if (exists $$samplefiles{$sample}{'coverage'}{'line'} ) {
			$line = $$samplefiles{$sample}{'coverage'}{'line'};
		}
		else {
			$line = $$samplefiles{$sample}{'coverage'}{'filehandle'} -> getline();
			$$samplefiles{$sample}{'coverage'}{'line'} = $line;
			
		}
		
		if ($$samplefiles{$sample}{'reference'}{'filehandle'} -> eof()) { # test for EOF # was $line
			$$samplefiles{$sample}{'coverage'}{'value'} = 'NA';
			$$samplefiles{$sample}{'coverage'}{'location'} = 'NA';
			$$samplefiles{$sample}{'coverage'}{'start'} = 'NA';
			$flag = 0;
			last;
		}
		
		if ($line =~ /^#/) {
			$$samplefiles{$sample}{'coverage'}{'line'} = $$samplefiles{$sample}{'coverage'}{'filehandle'} -> getline();
			next;
		}
		else {
			my @data = split("\t", $line);
			unless (scalar (@data) == 9) {
				print "coverage.0.gff.gz of $sample does not contain correct number of elements in line\n", $line, "\n";
				next;
			}
			
			my $start_position = encryptposition($data[0], $data[3]);
			my $end_position =  encryptposition($data[0], $data[4]);
			
			if ($end_position < $currentposition) {
				$$samplefiles{$sample}{'coverage'}{'line'} = $$samplefiles{$sample}{'coverage'}{'filehandle'} -> getline();
				next;
			}
			elsif ($start_position > $currentposition) {
				#print 'oops missed the correct position (', $currentposition, ') in coverage.0.gff.gz of ', $sample, ' in line', "\n", $line, "\n"; # DEBUG
				$$samplefiles{$sample}{'coverage'}{'value'} = 0;
				$$samplefiles{$sample}{'coverage'}{'location'} = $currentposition; 
				$$samplefiles{$sample}{'coverage'}{'start'} = $start_position; 
				$flag = 0;
				last;
			}
			elsif ($start_position <= $currentposition && $end_position >= $currentposition){
				$$samplefiles{$sample}{'coverage'}{'value'} = sprintf("%d", $data[5]); # integers are sufficiently precise
				$$samplefiles{$sample}{'coverage'}{'location'} = $currentposition; 
				$$samplefiles{$sample}{'coverage'}{'start'} = $start_position; 
				$flag = 0;
				last;
			}
		}
	}
}

sub clean {
	my $toClean = shift;
	$toClean =~ s/"//g;
	return($toClean);
}

sub closeAllFiles {
	my $samplefiles = shift;
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		$$samplefiles{$sample}{'coverage'}{'filehandle'} -> close();
		$$samplefiles{$sample}{'reference'}{'filehandle'} -> close();
		#$$samplefiles{$sample}{'SNP'}{'filehandle'} -> close();
		#$$samplefiles{$sample}{'allSNP'}{'filehandle'} -> close();
		$$samplefiles{$sample}{'Indel'}{'filehandle'} -> close();
		$$samplefiles{$sample}{'allIndel'}{'filehandle'} -> close();
	}
}


sub decifferposition {
	my $position = shift;
	my ($chr_position, $chr_number) = (0) x 2;
	$chr_position = $position % 1000000000;
	$chr_number = int ($position / 1000000000);
	return($chr_number, $chr_position);
}

sub deleteCondel{
	my $samplefiles = shift;
	#unlink ($$samplefiles{'allsamples'}{'condelin'}); #DEBUG
	#unlink ($$samplefiles{'allsamples'}{'condelout'}); #DEBUG
}
sub deleteSIFT {
	my $samplefiles = shift;
	
	unlink ($$samplefiles{'allsamples'}{'siftconf'});
	unlink ($$samplefiles{'allsamples'}{'siftinput'});
	unlink ($$samplefiles{'allsamples'}{'siftoutput'});
	
}
sub deleteUnzipped {
	my $samplefiles = shift;
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		unlink ($$samplefiles{$sample}{'coverage'}{'unzipped'});
		unlink ($$samplefiles{$sample}{'reference'}{'unzipped'});
	}
	unlink ($$samplefiles{'allsamples'}{'condelin'});
	unlink ($$samplefiles{'allsamples'}{'condelout'});
}

sub determineAllelFrequencies {
	my $samplefiles = shift;
	my $funcfilter = shift;
	my $exonicfuncfilter = shift;
		
	$$samplefiles{'allsamples'}{'AllelFreq'} = 0;
	$$samplefiles{'allsamples'}{'Freq'} = 0;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		next unless (defined($$samplefiles{$sample}{'allIndel'}{'line'}));
		my @data = split (/\t/, $$samplefiles{$sample}{'allIndel'}{'line'});
		my %values = ('GATK' => $data[9], 'MPILEUP' => $data[10], 'SHORE' => $data[11]);
		
		my $position = encryptposition($data[0], $data[1]);
		
		#ANNOVAR changes start positions "1       2938407 del_1_2938407   CAGA    C" gets transformed to "1,2938408,2938410,AGA,-"
		if (length($data[3]) > length($data[4])) { # ref is longer than observation -> deletion
			$data[1] += length($data[4]); # increment by length of observed nucleotide stretch, which would be done by ANNOVAR 
		}
		
		if ($position == $$samplefiles{'allsamples'}{'currentposition'}) {

			$$samplefiles{'allsamples'}{'Freq'} = 0;

			if ($$samplefiles{$sample}{'allIndel'}{'GATK'} ne 'NA') {
				if ($$samplefiles{$sample}{'allIndel'}{'GATK'} =~ /0\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 1;
				}
				elsif ($$samplefiles{$sample}{'allIndel'}{'GATK'} =~ /1\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 2;
				}
				else {
					
				}
				next; # use only GATK
			}
			
			if ($$samplefiles{$sample}{'allIndel'}{'MPILEUP'} ne 'NA') {
				if ($$samplefiles{$sample}{'allIndel'}{'MPILEUP'} =~ /0\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 1;
				}
				elsif ($$samplefiles{$sample}{'allIndel'}{'MPILEUP'} =~ /1\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 2;
				}
				else {
					
				}
				next; # second check MPILEUP
			}
			
			if ($$samplefiles{$sample}{'allIndel'}{'SHORE'} ne 'NA') {
				if ($$samplefiles{$sample}{'allIndel'}{'SHORE'} =~ /0\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 1;
				}
				elsif ($$samplefiles{$sample}{'allIndel'}{'SHORE'} =~ /1\/1/) {
					$$samplefiles{'allsamples'}{'AllelFreq'} += 2;
				}
				else {
					
				}
				next; # lastly use SHORE (if GATK and MPILEUP failed)
			}
		}
	}
}

sub determineFreq{
	my $positiveSamples = shift;
	my $SNPsamples = shift;
	
	my @positives = split(/\s/, $positiveSamples);
	my @negatives = split(/\s/, $SNPsamples);
	
	my $freq = 0;
	foreach my $sample (@positives) {
		unless ($sample eq 'NA') {
			$freq++;
		}
	}
	foreach my $sample (@negatives) {
		unless ($sample eq 'NA') {
			$freq++;
		}
	}
	
	#my $freq = scalar(@positives) + scalar(@negatives);
	
	return($freq);
	
}

sub determineGeneFrequencies {
	my $positiveSamples = shift;
	my $aachange = shift;
	my $allGENEdata = shift;
	my $funcfilter = shift;
	my $exonicfuncfilter = shift;
	my $exonic = shift;
	my $func = shift;
	
	if ($funcfilter) { 
		if ($func =~ /^.?$funcfilter.?$/i && ($exonic =~ /^.?$exonicfuncfilter.?$/i || $exonic  =~ /^.?$exonicfuncfilter SNV.?$/i)) {
	#	"skipped entry because of Funcfilter and Exonicfuncfilter"
			return;
		}
	}
	else {
		if ($exonic =~ /^.?$exonicfuncfilter.?$/i || $exonic  =~ /^.?$exonicfuncfilter SNV.?$/i) {
	#	"skipped entry because of Exonicfuncfilter"
			return;
		}
	}
	
	if ($aachange ne 'NA') {
		my @samples = split(/\s/, $positiveSamples);
		#how often is a gene hit and does it exceed the cutoff?
		my @aachanges = split(/:/, $aachange); # NM_152486:c.G503A:p.R168Q
		$$allGENEdata{$aachanges[0]} += scalar(@samples);
	}
	
	
};

sub doCondel{
	my $samplefiles = shift;
	my $condelcall = shift;
	
	system("$condelcall $$samplefiles{'allsamples'}{'condelin'} > $$samplefiles{'allsamples'}{'condelout'}");	
}

sub doSIFT{
	my $samplefiles = shift;
	my $siftcall = shift;
	my $siftconflocation = File::Spec->rel2abs( $$samplefiles{'allsamples'}{siftconf} );
	my $logfilelocation = File::Spec->rel2abs( $$samplefiles{'allsamples'}{siftconf} ).'log' ;
	#print "$siftcall $$samplefiles{'allsamples'}{siftconf} > $$samplefiles{'allsamples'}{siftconf}.log \n";
	#system("$siftcall $$samplefiles{'allsamples'}{siftconf} > $$samplefiles{'allsamples'}{siftconf}.log"); # DEBUG - REAL
	
	print "$siftcall $siftconflocation > $logfilelocation \n";
	system("$siftcall $siftconflocation > $logfilelocation");
	
	open(LOGFH, "<$logfilelocation"); # DEBUG - REAL
	my $file_content = do { local $/; <LOGFH> }; # DEBUG - REAL
	$file_content =~ /Job_ID:(\d+)/; #Job_ID:19620 # DEBUG - REAL
	$$samplefiles{'allsamples'}{'JobID'} = $1; # DEBUG - REAL
	
}

sub encryptposition {
	my $chr = shift;
	my $position = shift;	
	$chr = translate($chr); #X&Y to 23&24 
	my $encryptedPosition = $chr * 1000000000 + $position;
	
	return ($encryptedPosition);
}

sub enterDataForNegativeSamples {
	my $samplefiles = shift;
	my $negativeSamples = shift;
	my $delimiter = shift;
	my $currentposition = $$samplefiles{'allsamples'}{'currentposition'};
	
	foreach my $negativeSample (@{$negativeSamples}) {
		$$samplefiles{$negativeSample}{'data_pos'} = 0;
		$$samplefiles{$negativeSample}{'data_neg'} = 1;
		identifySNPposition($samplefiles, $negativeSample, $delimiter);
		
		#print "#####checkLowQualSNPs\n";
		checkLowQualSNPs($samplefiles, $negativeSample, $currentposition);
		#print "#####checkReference\n";
		checkReference($samplefiles, $negativeSample, $currentposition);
		#print "#####checkCoverage\n";
		checkCoverage($samplefiles, $negativeSample, $currentposition);
	}
}

sub enterDataForPositveSamples{
	my $samplefiles = shift;
	my $positiveSamples = shift;
	my $allGENEdata = shift;
	my $delimiter = shift;
	my $currentposition = $$samplefiles{'allsamples'}{'currentposition'};
	
	foreach my $positiveSample (@{$positiveSamples}) {
		$$samplefiles{$positiveSample}{'data_pos'} = 1;
		$$samplefiles{$positiveSample}{'data_neg'} = 0;
		
		checkLowQualSNPs($samplefiles, $positiveSample, $currentposition); # used for Allelfrequency count
		checkReference($samplefiles, $positiveSample, $currentposition); # used for sample details
		checkCoverage($samplefiles, $positiveSample, $currentposition);
		identifySNPposition($samplefiles, $positiveSample, $delimiter);# update position for next round
	}
}

sub enterEVSdata {
	my $samplefiles = shift;
	my $evsfilefound = shift;
	
	unless ($evsfilefound) {
		fillEVSvalues($samplefiles);
	}
	
	else {
		if ($$samplefiles{'allsamples'}{'evsfilehandle'} -> eof()) {
			fillEVSvalues($samplefiles);
			return;
		}
		
		
		my $currentposition = $$samplefiles{'allsamples'}{'currentposition'};
		my $flag = 1;
		
		
		
		while($flag) {
			my $line = $$samplefiles{'allsamples'}{'evsfilehandle'} -> getline();
			if (!defined($line) || $line eq '') { #
				fillEVSvalues($samplefiles);
				print "read error in evs file\n"; #
				$flag = 0;
				last;
			}
			
			chomp($line);
			next if ($line =~ /^#/);
			
			my @splitline = split(/\t/, $line); # at this position ^I and \t is not recognized.
			my $chr = translate($splitline[0]);
			my $position = encryptposition($chr, $splitline[1]);
			
			if ($position < $currentposition) {
				next;
			}
			elsif ($position > $currentposition) {
				$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'} = 'NA';
				$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'} = 'NA';
				$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'} = 'NA';
				$samplefiles{'allsamples'}{'evsdata'}{'FG'} = 'NA'; #functionGVS
				$samplefiles{'allsamples'}{'evsdata'}{'PH'} = 'NA'; #PolyPhen
				
				$flag = 0;
				last;
			}
			elsif ($position == $currentposition) {
				$splitline[7] =~ /MAF=(.+?);/;
				my $MAF = $1;
				$MAF =~ s/,/|/g;
				$splitline[7] =~ /FG=(.+?);/;
				my $FG = $1;
				$FG =~ s/,/|/g;
				$splitline[7] =~ /PH=(.+)/;
				my $PH = $1;
				$PH =~ s/,/|/g;
				
				
				my @evsmaf = split(/\|/, $MAF);
				
				# EA_AC=0,248;AA_AC=1,363;TAC=1,611
				# allel counts
				# european americans
				$splitline[7] =~ /EA_AC=(.+?);/;
				my $EA_AC = $1;
				my @EA_ACs = split(/\,/, $EA_AC);
				#african americans
				$splitline[7] =~ /AA_AC=(.+?);/;
				my $AA_AC = $1;
				my @AA_ACs = split(/\,/, $AA_AC);
				#total allel count
				$splitline[7] =~ /TAC=(.+?);/;
				my $TAC = $1;
				my @TACs = split(/\,/, $TAC);
				
				#if ($EA_ACs[0] < $EA_ACs[1]) { # if the reference is never seen
				if (($EA_ACs[0] + $EA_ACs[1]) != 0) {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'}  = sprintf("%.3f", $EA_ACs[0] / ($EA_ACs[0] + $EA_ACs[1]));
				}
				else {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'}  = 0;
				}
				if (($AA_ACs[0] + $AA_ACs[1]) != 0) {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'}  = sprintf("%.3f", $AA_ACs[0] / ($AA_ACs[0] + $AA_ACs[1]));
				}
				else {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'}  = 0;
				}
				if (($TACs[0] + $TACs[1]) != 0) {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'}  = sprintf("%.3f", $TACs[0] / ($TACs[0] + $TACs[1]));
				}
				else {
					$$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'}  = 0;
				}
				# if the reference is never seen (now with ternary factors)
				#$$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'}  = ($AA_ACs[0] < $AA_ACs[1]) ? sprintf("%.3f", $AA_ACs[0] / ($AA_ACs[0] + $AA_ACs[1])) : sprintf("%.3f", 1 - $AA_ACs[0] / $AA_ACs[1]);
				#$$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'} = ($TACs[0] < $TACs[1]) ? sprintf("%.3f", $TACs[0] / $TACs[1]) : sprintf("%.3f", 1 - $TACs[0] / $TACs[1]);
				$samplefiles{'allsamples'}{'evsdata'}{'FG'} = $FG; #functionGVS
				$samplefiles{'allsamples'}{'evsdata'}{'PH'} = $PH; #PolyPhen
				
				#if ($TACs[0] < $TACs[1]) { #EVS always provides minor allel frequency even if minor allel is reference
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'} = $evsmaf[0] / 100; # adapt 100% values of EVS to 0.01 values of 1000Genomes
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'} = $evsmaf[1] / 100;
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'} = $evsmaf[2] / 100;
				#}
				#else {
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'} = 1 - ($evsmaf[0] / 100); # adapt 100% values of EVS to 0.01 values of 1000Genomes
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'} = 1 - ($evsmaf[1] / 100);
				#	$$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'} = 1 - ($evsmaf[2] / 100);
				#}
				
				$flag = 0;
				last;
			}
		}
	}

}

sub enterMutationAssessor {
	my $samplefiles = shift;
	my $mutassfilefound = shift;
	my $delimiter = shift;
	my @annovardata;
	$annovardata[15] = 0; # initialize only position 15 and 16, because they are needed later on. Rest will be filled at SNP read
	$annovardata[16] = 0; # this is dirty.
	
	unless ($mutassfilefound == 1) {
		fillMutAssvalues($samplefiles);
	}
	
	else {
		if ($$samplefiles{'allsamples'}{'mutassfilehandle'} -> eof()) {
			fillMutAssvalues($samplefiles);
			return;
		}
		
		my $currentposition = $$samplefiles{'allsamples'}{'currentposition'};
		my $flag = 1;
		
		while($flag) {
			# read line and check if it worked
			my $line = $$samplefiles{'allsamples'}{'mutassfilehandle'} -> getline();
			if (!defined($line) || $line eq '') { #
				fillMutAssvalues($samplefiles);
				print "read error in mutation assessor file\n"; #
				$flag = 0;
				last;
			}
			
			chomp($line);
			
			#calculate position read from mutationassessor file
			my @splitline = split(/\,/, $line); 
			my $position = encryptposition($splitline[1], $splitline[2]);
			
			# fetch one SNP line to identify which nucleotide change has been seen and if it's the same in mutationassessor
			if ($position != encryptposition(translate($annovardata[15]), $annovardata[16])) {
				foreach my $sample (keys %{$samplefiles}) {
					next if $sample eq 'allsamples';
					@annovardata = split($delimiter, $$samplefiles{$sample}{'SNP'}{'lastline'});
					next if ($position != encryptposition(translate($annovardata[15]), $annovardata[16]));
					last;
				}
			}
			
			#
			if ($position < $currentposition) {
				next;
			}
			elsif ($position > $currentposition) {
				$$samplefiles{'allsamples'}{'mutass'} = 'NA';
				
				$flag = 0;
				last;
			}
			elsif ($position == $currentposition && ($annovardata[18] ne $splitline[3] || $annovardata[19] ne $splitline[4])) { #$annovardata[18] = ref, $annovardata[19] = observed nucleotide
				next;
			}
			elsif ($position == $currentposition && $annovardata[18] eq $splitline[3] && $annovardata[19] eq $splitline[4]) {
				if (defined($splitline[6])) {# mutationassessor value
					$$samplefiles{'allsamples'}{'mutass'} = $splitline[6];
					$flag = 0;
				}
				else {
					$$samplefiles{'allsamples'}{'mutass'} = 'NA';
				}
			}
			else {
				fillMutAssvalues($samplefiles);
			}
		}
	}
}

sub filenotfoundexception {
	my $snpfile = shift;
	
	print "\n\n", 'W A R N I N G: could not find sum.exome_summary.csv neither in AnnovarIntersection nor Annovar folder.', "\n";
	print $snpfile, "\n";
	exit;
}

sub filetest {
	my $testfile = shift;
	my $errormessages_ref = shift;
	my @errormessages = @$errormessages_ref;
	unless (-e $testfile) {
		print "$testfile does not exist \n", "proceeding to next sample.\n";
		reporterrors(\@errormessages, "$testfile does not exist \n");
		return 0;
	}
	else {
		return 1;
	}
	

}

sub fillEmptyValues {
	my $printout = shift;
	my @musthave = ('chromosome', 'nucleotide', 'positivesamples', 'samplereference', 'negativelowqualsnp', 'samplecoverage', 'func', 'gene', 'exonicfunc', 'aachange', 'conserved', 'segdup', 'maf_1000g', 'evsmaf_ea', 'evsmaf_aa', 'evsmaf_tac', 'evsfg', 'evsph', 'dbsnp', 'sift', 'polyphen', 'ljb_phylop', 'ljb_mutationtaster', 'ljb_lrt', 'ref', 'obs', 'qualitymedian');
	foreach my $obligary (@musthave) {
		unless(defined($$printout{$obligary})) {$$printout{$obligary} = 'NA'};
	}
	
	foreach my $key ( keys %{$printout}) {
			unless(defined($$printout{$key})) {$$printout{$key} = 'NA'}; # undefined values get compensated by 'NA'
			
			if (($$printout{$key} eq 'NA' || $$printout{$key} eq '') && ($key eq 'sift' || $key eq 'ljb_mutationtaster' || $key eq 'ljb_lrt')) { # fill empty "must be numeric" cells
				$$printout{$key} = 2;
			}
			elsif (($$printout{$key} eq 'NA' || $$printout{$key} eq '') && ($key eq 'exonicfunc' || $key eq 'aachange' || $key eq 'dbsnp' || $key eq 'ref' || $key eq 'obs' || $key eq 'other' || $key eq 'evsph' || $key eq 'evsfg')) { # fill in NA into cells, which should contain text
				$$printout{$key} = 'NA';
			}
			elsif (($$printout{$key} eq 'NA' || $$printout{$key} eq '') && ($key eq 'maf_1000g' || $key eq 'evsmaf_ea' || $key eq 'evsmaf_aa' || $key eq 'evsmaf_tac' || $key eq 'segdup')) { # fill in 0 into MAF cells, hasn't been seen before so MAF = 0 and EVS values. no segduppp value is good so it gets a 0. as well # TODO check if evs ph and evs fg are 0 based
				$$printout{$key} = 0;
			}
			elsif (($$printout{$key} eq 'NA' || $$printout{$key} eq '') && ($key eq 'polyphen' || $key eq 'conserved' || $key eq 'ljb_phylop')) { # fill in -1 into Polyphen2 cells, it's megabenign
				$$printout{$key} = -1;
			}
			#elsif (($$printout{$key} eq 'NA' || $$printout{$key} eq '') && ()) { # fill in -1 into Polyphen2 cells, it's megabenign
			#	$$printout{$key} = 1000;
			#}
		}
	if ($$printout{'sift'} == 2 && $$printout{'exonicfunc'} =~ /stoploss|stopgain|splicing/ ) {
		$$printout{'sift'} = 0;
	}
	if ($$printout{'polyphen'} == 2 && $$printout{'exonicfunc'}  =~ /stoploss|stopgain|splicing/ ) {
		$$printout{'polyphen'} = 1;
	}
}

sub fillEVSvalues{
	my $samplefiles = shift;
	$samplefiles{'allsamples'}{'evsdata'}{'MAF_ea'} = 'NA';
	$samplefiles{'allsamples'}{'evsdata'}{'MAF_aa'} = 'NA';
	$samplefiles{'allsamples'}{'evsdata'}{'MAF_tac'} = 'NA';
	$samplefiles{'allsamples'}{'evsdata'}{'FG'} = 'NA'; #functionGVS
	$samplefiles{'allsamples'}{'evsdata'}{'PH'} = 'NA'; #PolyPhen
};

sub fillLastSNPline {
	my $samplefiles = shift;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		$$samplefiles{$sample}{'Indel'}{'lastline'} = $$samplefiles{$sample}{'Indel'}{'filehandle'} -> getline();
		while ($$samplefiles{$sample}{'Indel'}{'lastline'} =~ /^Func.*/) {
			$$samplefiles{$sample}{'Indel'}{'lastline'} = $$samplefiles{$sample}{'Indel'}{'filehandle'} -> getline();
		}
		
		
		### replace wrong , letters in lastline
		my $line = $$samplefiles{$sample}{'Indel'}{'lastline'};
		chomp($line);
		my $newline;
		my $doubleqoteflag = 0;
		for (my $key = 0; $key < length($line); $key++) {
			my $letter = substr($line, $key, 1);
			$doubleqoteflag++ if ($letter eq '"');
			$letter = ';' if ($letter eq ',' && $doubleqoteflag % 2);
			$newline .= $letter;
		}
		$$samplefiles{$sample}{'Indel'}{'lastline'} = $newline;
		
		
	}
}

sub fillMutAssvalues{
	my $samplefiles = shift;
	$samplefiles{'allsamples'}{'mutass'} = 'NA';
};

sub fillSampleFields {
	my $printout = shift;
	
	unless ($$printout{'positivesamples'}) {$$printout{'positivesamples'} = 'NA'};
	unless ($$printout{'samplereference'}) {$$printout{'samplereference'} = 'NA'};
	unless ($$printout{'samplereference_lowqual'}) {$$printout{'samplereference_lowqual'} = 'NA'};
	unless ($$printout{'negativelowqualsnp'}) {$$printout{'negativelowqualsnp'} = 'NA'};
}

sub filterGeneFrequency {
	my $aachange = shift;
	my $allGENEdata = shift;
	my $genecutoff = shift;
	
	
	my @aachangedata = split(/:/, $aachange);
	if ($$allGENEdata{$aachangedata[0]}) {
		if ($$allGENEdata{$aachangedata[0]} >= $genecutoff) {
			return $$allGENEdata{$aachangedata[0]};
		}
		else {
			return 0;
		}
	}
	else {
		return 'NA';
	}
}

sub filterValues{
	my $data = shift;
	my $MAFcutoff = shift;
	my $funcfilter = shift;
	my $exonicfuncfilter = shift;

	#====== apply filter for 1000 genome MAF =========
# 	if not found in 1000 Genomes Project it's interesting
	
	if ($$data[16] > $MAFcutoff) {  #ignore SNPs that are more frequent in 1000Genomes than --MAFcutoff
		return 0;
	}
	#====== apply filter for Func and ExonicFunc =========
	if ($funcfilter) { # {'Func'} = $data[9]; {'ExonicFunc'} = $data[11];
		if ($$data[9] =~ /^.?$funcfilter.?$/i && ($$data[11] =~ /^.?$exonicfuncfilter.?$/i || $$data[11]  =~ /^.?$exonicfuncfilter SNV.?$/i)) {
#			"skipped entry because of Funcfilter and Exonicfuncfilter"
			return 0;
		}
	}
	else {
		if ($$data[11] =~ /^.?$exonicfuncfilter.?$/i || $$data[11]  =~ /^.?$exonicfuncfilter SNV.?$/i) {
#			"skipped entry because of Funcfilter"
			return 0;
		}
	}
	
	return 1;
};

sub groupvalues{
	my $alldata = shift;
		
	foreach my $position (keys %{$alldata}) {
		#my $maxMAF = ;
		my $MAF = sprintf("%.3f", $$alldata{$position}{'MAF_1000G'});
		$$alldata{$position}{'MAF_bin'} = int( $MAF / 0.001); # find the biggest value and form bins. (good 1-1000 bad)
		my $SegDup = sprintf("%.2f", $$alldata{$position}{'SegDup'} );
		$$alldata{$position}{'SegDup_bin'} = (int( $SegDup / 0.01) + 1) * 10; # form bins for segmental duplications ( good 10-1000 bad)
		#my $ConDel = sprintf("%.3f", $$alldata{$position}{'CondelValue'} );
		#$$alldata{$position}{'Condel_bin'} = int( $ConDel  / 0.001); # binning of condel values (bad 1-100 good)
		$$alldata{$position}{'Conserved_bin'} =  $$alldata{$position}{'Conserved'};#int( ( 1001 - ( $$alldata{$position}{'Conserved'} + 1 ) ) / 10) + 1; # (1-100 means that this position is conserved) [missing values had -1, but are put into lowest rank] no binning since 1000 stages are good
		my $PhyloP = sprintf("%.3f", $$alldata{$position}{'LJB_PhyloP'} );
		if ($PhyloP != -1.000) {
			$$alldata{$position}{'LJB_PhyloP_bin'} = int(  $PhyloP  / 0.001) ; # a value of 1 means, that the position is conserved [bad 1 - 1000 good]
		}
		else {
			$$alldata{$position}{'LJB_PhyloP_bin'} = 0;
		}
	}
}

sub identifyCurrentPosition {
	my $samplefiles = shift;
	my $delimiter = shift;
	my @positions;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		my @data = split($delimiter , $$samplefiles{$sample}{'Indel'}{'lastline'});
		
		print $sample, "\n";
		#print Dumper(@data); # DEBUG
		
		my $location = encryptposition($data[15], $data[16]);
		push (@positions, $location);
	}
	
	# calculate min position
	my $min_value = List::Util::min(@positions);
	# return min position
	return ($min_value);
	
}

sub identifyNextPosition{
	my $samplefiles = shift; # is % ref
	my $delimiter = shift; 
	my $flag = shift; # is $ ref
	my $positivesamples = shift; # is @ ref
	
	foreach my $sample (@{$positivesamples}) {
		next if ($sample eq 'allsamples');
		my $line = $$samplefiles{$sample}{'Indel'}{'filehandle'} -> getline();
		while (!$line) {
			if ($$samplefiles{$sample}{'Indel'}{'filehandle'} -> eof()) {
				$$samplefiles{$sample}{'Indel'}{'foundposition'} = 'NA';
				$line = 'eof';
				last;
			}
			$line = $$samplefiles{$sample}{'Indel'}{'filehandle'} -> getline(); # avoid $line to be empty
		}
		$$samplefiles{$sample}{'Indel'}{'lastline'} = $line;
		
		### replace wrong , letters in lastline
		$line = $$samplefiles{$sample}{'Indel'}{'lastline'};
		chomp($line);
		my $newline;
		my $doubleqoteflag = 0;
		for (my $key = 0; $key < length($line); $key++) {
			my $letter = substr($line, $key, 1);
			$doubleqoteflag++ if ($letter eq '"');
			$letter = ';' if ($letter eq ',' && $doubleqoteflag % 2);
			$newline .= $letter;
		}
		$$samplefiles{$sample}{'Indel'}{'lastline'} = $newline;
		
		
		identifySNPposition($samplefiles, $sample, $delimiter);
	}
	
	my @allpositions;
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		next if ($$samplefiles{$sample}{'Indel'}{'foundposition'} eq 'NA');# skip "end of file" entries
		push(@allpositions, $$samplefiles{$sample}{'Indel'}{'foundposition'});
	}
	
	### skip "end of file" entries
	if (scalar(@allpositions) == 0) {
		$$flag = 0;
		return;
	}
	### done - skip "end of file" entries
	
	my $min_value = List::Util::min(@allpositions);
	return ($min_value);
}

sub identifyQualityMedian {
	my $samplefiles = shift;
	my $delimiter = shift;
	my @qualities;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		if ($$samplefiles{$sample}{'data_pos'}) {
			my @data = split($delimiter, $$samplefiles{$sample}{'Indel'}{'lastline'});
			my $value = $data[21];
			$value =~ s/\"//g;
			push(@qualities, $value);
		}
		
		my @qualities_sorted = sort { $a <=> $b } @qualities;
		$$samplefiles{'allsamples'}{'quality'} = $qualities_sorted[int(scalar(@qualities_sorted) / 2)]; # calculate median is the middle element in a sorted array
	}
}

sub identifySamples{
	my $samplefiles = shift;
	my $currentPosition = $$samplefiles{'allsamples'}{'currentposition'};
	my $delimiter = shift;
	my $positiveSamples = shift;
	my $negativeSamples = shift;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		if ($$samplefiles{$sample}{'Indel'}{'lastline'} eq 'eof') {
			push (@{$negativeSamples}, $sample);
			next;
		};
		my @data = split($delimiter , $$samplefiles{$sample}{'Indel'}{'lastline'});
		my $location = encryptposition($data[15], $data[16]);
		if ($location == $currentPosition) {
			push (@{$positiveSamples}, $sample);
		}
		else {
			push (@{$negativeSamples}, $sample);
		};
	}	
}

sub identifySNPposition {
	my $samplefiles = shift;
	my $sample = shift;
	my $delimiter = shift;
	
	unless ($$samplefiles{$sample}{'Indel'}{'lastline'} eq 'eof') {
		# annovar separates gene names by , - have to get rid of them before splitting
		#replace , in ""
		#sometimes lines look like this "splicing","NCAPH2(NM_001185011:exon14:c.1233+2C>T,NM_152299:exon14:c.1233+2C>T)",,,
		#the , between C>T and NM_ has to be replaced
		if ($delimiter eq ',') {
			my $line = $$samplefiles{$sample}{'Indel'}{'lastline'};
			my $newline;
			my $doubleqoteflag = 0;
			for (my $key = 0; $key < length($line); $key++) {
				my $letter = substr($line, $key, 1);
				$doubleqoteflag++ if ($letter eq '"');
				$letter = ';' if ($letter eq ',' && $doubleqoteflag % 2);
				$newline .= $letter;
			}
			$$samplefiles{$sample}{'Indel'}{'lastline'} = $newline;
		}
		
		my @data = split($delimiter, $$samplefiles{$sample}{'Indel'}{'lastline'});
		my $position = encryptposition($data[15], $data[16]);
		$$samplefiles{$sample}{'Indel'}{'foundposition'} = $position;
	}
	else {
		$$samplefiles{$sample}{'Indel'}{'foundposition'} = 'NA';
	}
}

sub isExcludedGene{
	my $gene = shift;
	my $excludedgenes = shift;
	
	# sometimes Annovar puts multiple gene names into gene field
	my @geneNames = split(/\;/, $gene);
	my $judgement = 0;
	foreach my $singleGene (@geneNames) {
		if (defined($$excludedgenes{$singleGene}) && $$excludedgenes{$singleGene} == 1) {
			$judgement = 1;
		}
	}
	
	if ($judgement == 1) {
		return(1);
	}
	else {
		return(0);
	}
}

sub missingvalue {
	my $value = shift;
	my $sample = shift;
	my $filetype = shift;
	
	unless ($value) {
		print STDERR 'missing value in sample ', $sample, ' in file for ', $filetype, ' entering default value [no call]. Please check your data.', "\n";
		return 0;
	}
	else {
		return 1;
	}
}

sub nothingToDo {
	my $samplefiles = shift;
	my @non_eofsamples;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		if (defined ($$samplefiles{$sample}{'Indel'}{'lastline'}) && $$samplefiles{$sample}{'Indel'}{'lastline'} ne 'eof') {
			push @non_eofsamples, $sample;
		}
	}
	
	if (scalar @non_eofsamples == 0) {
		return 0;
	}
	else {
		return 1;
	}
}

sub openAllFiles {
	my $samplefiles = shift;
	my $outfile = shift;
	my $EVSfile = shift;
	my $evsfilefound = shift;
	my $mutationAssessor = shift;
	my $mutassfilefound = shift;
	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples');
		my $coverageFH = new IO::File($$samplefiles{$sample}{'coverage'}{'unzipped'}, O_RDONLY) or die "$$samplefiles{$sample}{'coverage'}{'unzipped'} $!";
		my $referenceFH = new IO::File($$samplefiles{$sample}{'reference'}{'unzipped'}, O_RDONLY) or die "$$samplefiles{$sample}{'reference'}{'unzipped'} $!";
		#my $snpFH = new IO::File($$samplefiles{$sample}{'SNP'}{'name'}, O_RDONLY) or die "$$samplefiles{$sample}{'SNP'}{'name'} $!";
		#my $allsnpFH = new IO::File($$samplefiles{$sample}{'allSNP'}{'name'}, O_RDONLY) or die "$$samplefiles{$sample}{'allSNP'}{'name'} $!";
		$$samplefiles{$sample}{'coverage'}{'filehandle'} = $coverageFH;
		$$samplefiles{$sample}{'reference'}{'filehandle'} = $referenceFH;
		#$$samplefiles{$sample}{'SNP'}{'filehandle'} = $snpFH;
		#$$samplefiles{$sample}{'allSNP'}{'filehandle'} = $allsnpFH;
		
		my $indelFH = new IO::File($$samplefiles{$sample}{'Indel'}{'name'}, O_RDONLY) or die "$$samplefiles{$sample}{'Indel'}{'name'} $!";
		my $allindelFH = new IO::File($$samplefiles{$sample}{'allIndel'}{'name'}, O_RDONLY) or die "$$samplefiles{$sample}{'allIndel'}{'name'} $!";
		$$samplefiles{$sample}{'Indel'}{'filehandle'} = $indelFH;
		$$samplefiles{$sample}{'allIndel'}{'filehandle'} = $allindelFH;
	}
	
	my $outfile_tmp = $outfile.".collection";
	my $collectionFH = new IO::File($outfile_tmp, O_EXCL|O_CREAT|O_WRONLY) or die "\n\n $outfile_tmp $! \n\n"; #REAL
	$$samplefiles{'allsamples'}{'collectionfilehandle'} = $collectionFH; #REAL
	
	my $EVSfh;
	if ($EVSfh = new IO::File($EVSfile, O_RDONLY)) {
		$$samplefiles{'allsamples'}{'evsfilehandle'} = $EVSfh;
		$$evsfilefound = 1;
	}
	else {
		print "Could not find $EVSfile will enter \"NA\" in EVS related fields. \n";
		$$evsfilefound = 0;
	}
	
	my $mutassfh;
	if ($mutassfh = new IO::File($mutationAssessor, O_RDONLY)) {
		$$samplefiles{'allsamples'}{'mutassfilehandle'} = $mutassfh;
		$$mutassfilefound = 1;
	}
	else {
		print "Could not find $mutationAssessor. Will enter \"NA\" in Mutation Assessor related fields. \n";
		$$mutassfilefound = 0;
	}
	
	### possible file permissions ###
	#O_RDWR		Read and Write
	#O_RDONLY	Read Only
	#O_WRONLY	Write Only
	#O_CREAT	Create the file
	#O_APPEND	Append the file
	#O_TRUNC	Truncate the file
	#O_EXCL		Stops if file already exists
	#O_NONBLOCK	Non-Blocking usability
	
};

sub openOutFile {
	my $samplefiles = shift;
	my $outfile = shift;
	my $exclusionList = shift;
	my $useexclusionlist = shift;
	my $random_number = int(rand(1000000));
	$$samplefiles{'allsamples'}{random} = $random_number;
	
	my $outfile_tmp = $outfile.".collection";
	my $collectionFH = new IO::File($outfile_tmp, O_RDWR) or die "\n\n $outfile_tmp $! \n\n";
	$$samplefiles{'allsamples'}{'collectionfilehandle'} = $collectionFH;
	
	my $outfile_final = $outfile.".final";
	my $outFH = new IO::File($outfile_final, O_EXCL|O_CREAT|O_WRONLY) or die "\n\n $outfile_final $! \n\n";
	$$samplefiles{'allsamples'}{'outfilehandle'} = $outFH;
	
	my $outfile_sift = $outfile.$random_number.".sift";
	
	my $siftFH = new IO::File($outfile_sift, O_EXCL|O_CREAT|O_WRONLY) or die "\n\n $outfile_sift $! \n\n";
	$$samplefiles{'allsamples'}{'siftfilehandle'} = $siftFH;
	$$samplefiles{'allsamples'}{'siftfilelocation'} = File::Spec->rel2abs( $outfile_sift );
	
	my $outfile_siftconf = $outfile.$random_number.".siftconf";
	my $siftconfFH = new IO::File($outfile_siftconf, O_EXCL|O_CREAT|O_WRONLY) or die "\n\n $outfile_siftconf $! \n\n";
	$$samplefiles{'allsamples'}{'siftconffilehandle'} = $siftconfFH;
	
	$$samplefiles{'allsamples'}{'siftconf'} = $outfile_siftconf;
	$$samplefiles{'allsamples'}{'siftinput'} = $outfile_sift;
	$$samplefiles{'allsamples'}{'siftoutput'} = $outfile.$random_number.".siftout";
	
	if ($useexclusionlist == 1) {
		my $exclusionFH = new IO::File($exclusionList, O_RDONLY) or die "\n\n $exclusionList $! \n\n";
		$$samplefiles{'allsamples'}{'exclusionfilehandle'} = $exclusionFH;
	}
	else {
		$$samplefiles{'allsamples'}{'exclusionfilehandle'} = 'none';
	}
}

sub preparenegativename {
	my $negativesample = shift;
	my $gatk = shift;
	my $mpileup = shift;
	my $shore = shift;
	
	if ($gatk !~ /^\.\/\./ && $gatk ne 'NA') {
		$negativesample = $negativesample.'_GATK';
	}
	if ($mpileup !~ /^\.\/\./ && $mpileup ne 'NA') {
		$negativesample = $negativesample.'_MPILEUP';
	}
	if ($shore !~ /^\.\/\./ && $shore ne 'NA') {
		$negativesample = $negativesample.'_SHORE';
	}
	
	return $negativesample;
}

sub rank {
	my $parameter = shift;
	my $order = shift;
	my $rankname = shift;
	my $alldata_ref = shift;
	my %alldata = %$alldata_ref;
	my $rank = 0;
	my $old_pos;
	
	if ($order eq 'ascending') {
		foreach my $pos ( sort { $alldata{$a}{$parameter} <=> $alldata{$b}{$parameter} } (keys %alldata)) {
			# $grades{$a} <=> $grades{$b}; would be the easy way to sort
			#====== initialize in first round ====
			if ($rank == 0) {
				$old_pos = $pos;
				$rank = 1;
			}
			#====== do ranking and give same ranks to equal values ====
			unless ($alldata{$pos}{$parameter} == $alldata{$old_pos}{$parameter}) {
				$rank++;
			}
			$alldata{$pos}{$rankname} = $rank; # enter rank value to %alldata
		
			$old_pos = $pos; # prepare next position
		}
	}
	elsif ($order eq 'descending') {
		foreach my $pos ( sort { $alldata{$b}{$parameter} <=> $alldata{$a}{$parameter} } (keys %alldata)) {
			#====== initialize in first round ====
			if ($rank == 0) {
				$old_pos = $pos;
				$rank = 1;
			}
			#====== do ranking and give same ranks to equal values ====
			unless ($alldata{$pos}{$parameter} == $alldata{$old_pos}{$parameter}) {
				$rank++;
			}
			$alldata{$pos}{$rankname} = $rank; # enter rank value to %alldata
		
			$old_pos = $pos; # prepare next position
		}
	}
	else {
		die 'error in sorting parameter. check parameters for rank subroutine! Its neither ascending nor descending';
	}
}

sub rankproduct{
	my $alldata = shift;
	
	foreach my $position (keys %{$alldata}) {
		my $toolcount = 1; # will contain the maximum amount of ranks that can be reached
		# take care on conservation values - if everything is not set and a segmental duplication is present put in the worst rank.
		my $conservation_rank;
		my $tool;
		if ($$alldata{$position}{'LJB_PhyloP'} != -1) {
			$conservation_rank = $$alldata{$position}{'LJB_PhyloP_Rank'};
			$toolcount *= 1000;
			$tool = 'LJB_PhyloP';
		}
		elsif ($$alldata{$position}{'Conserved'} != -1) {
			$conservation_rank = $$alldata{$position}{'Conserved'};
			$toolcount *= 1000;
			$tool = 'Conserved';
		}
		elsif ($$alldata{$position}{'Conserved'} == -1 && $$alldata{$position}{'SegDup'} >= 0.8 ) {
			$conservation_rank = 1000;
			$toolcount *= 1000;
			$tool = 'noToolHighSegDup';
		}
		else {
			$conservation_rank = 1000;
			$toolcount *= 1000;
			$tool = 'noTool';
		}
		$$alldata{$position}{'Conservation_Rank_Used'} = $conservation_rank;
		$$alldata{$position}{'Conservation_Tool_Used'} = $tool;
		
		if ($$alldata{$position}{'SegDup_Rank'} =~ /\d*/) {
			$toolcount *= 1000;
		}
		
		if ($$alldata{$position}{'MAF_Rank'} =~ /\d*/) {
			$toolcount *= 1000;
		}
		
		#my $condel_rank;
		#if ($$alldata{$position}{'Condel_Rank'} =~ /\d*/ && $$alldata{$position}{'Condel_Rank'} ne 'not_computable') {
		#	$toolcount *= 1000;
		#	$condel_rank = $$alldata{$position}{'Condel_Rank'};
		#}
		#else { # would mean 'not_computable'
		#	$condel_rank = 1; # actually removing condel from product, since it gets neutral
		#	# and do not increase $toolcount
		#}
		
		if ($$alldata{$position}{'SegDup_Rank'} =~ /\d*/ && $$alldata{$position}{'MAF_Rank'} =~ /\d*/  &&  $conservation_rank =~ /\d*/) {# &&  $$alldata{$position}{'Condel_Rank'} =~ /\d*/
			$$alldata{$position}{'RankProduct'} = ( $$alldata{$position}{'SegDup_Rank'} * $$alldata{$position}{'MAF_Rank'}  * $conservation_rank ) / $toolcount;# * $condel_rank
		}
	}
};

sub readCondelOutput{
	my $samplefiles = shift;
	my $alldata = shift;
	open (CONDEL, "<$$samplefiles{'allsamples'}{'condelout'}");
	while(<CONDEL>) {
		chomp;
		my @data = split(/\t/);
		if ($data[1] eq 'not_computable_was') {
			$$alldata{$data[0]}{'CondelValue'} = 0;
			$$alldata{$data[0]}{'CondelJudgement'} = 'not_computable';
		}
		else {
			$$alldata{$data[0]}{'CondelValue'} = $data[1];
			$$alldata{$data[0]}{'CondelJudgement'} = $data[2];
		}
	}
}

sub readExclusionList{
	my $samplefiles = shift;
	my $excludedgenes = shift;
	
	if ($$samplefiles{'allsamples'}{'exclusionfilehandle'} ne 'none') {
		while (my $line = $$samplefiles{'allsamples'}{'exclusionfilehandle'} -> getline()) {
			if (!defined($line) || $line eq '') {
				if ($$samplefiles{'allsamples'}{'exclusionfilehandle'} -> eof) {
					last;
				}
				next;
			}
			else {
				chomp($line);
				$$excludedgenes{$line} = 1;
			}
			if ($$samplefiles{'allsamples'}{'exclusionfilehandle'} -> eof) {
				last;
			}
			
		}
	}
	else {
		return;
	}
	
	
	
}

sub readInitialOutput {
	my $samplefiles = shift;
	my $alldata = shift;
	my $delimiter = shift;
	my $allGENEdata = shift;
	my $genecutoff = shift;
	my $MAFcutoff = shift;
	my $funcfilter = shift;
	my $exonicfuncfilter = shift;
	
	while (my $line = $$samplefiles{'allsamples'}{'collectionfilehandle'} -> getline()) {
		next if ($line =~ /PositiveSamples/); #skip header
		last unless ($line);
		my @data = split ($delimiter, $line);
		
		next unless (filterValues(\@data, $MAFcutoff, $funcfilter, $exonicfuncfilter) );
		
		
		my $position = encryptposition($data[0], $data[1]);
		
		$$alldata{$position}{'endposition'} = $data[2];
		$$alldata{$position}{'PositiveSamples'} = $data[3];
		$$alldata{$position}{'NegativeSamplesRefCall'} = $data[4];
                $$alldata{$position}{'NegativeSamplesRefCall_lowqual'} = $data[5];
		$$alldata{$position}{'NegativeSamplesSNPCall'} = $data[6];
		$$alldata{$position}{'AverageCoverage'} = $data[7];
		$$alldata{$position}{'SamplesCoverage'} = $data[8];
		$$alldata{$position}{'AllelFreq'} = $data[9];
		$$alldata{$position}{'Freq'} = $data[10];
		$$alldata{$position}{'Func'} = $data[11];
		$$alldata{$position}{'Gene'} = $data[12];
		$$alldata{$position}{'ExonicFunc'} = $data[13];
		$$alldata{$position}{'AAChange'} = $data[14];
		$$alldata{$position}{'Conserved'} = $data[15];
		$$alldata{$position}{'LJB_PhyloP'} = $data[16];
		#$$alldata{$position}{'LJB_PhyloP_Rank'} = 'NA';
		$$alldata{$position}{'SegDup'} = $data[17];
		$$alldata{$position}{'MAF_1000G'} = $data[18];
		#$$alldata{$position}{'MAF_Rank'} = 'NA';
		#$$alldata{$position}{'EVS_MAF_EA'} = $data[18];
		#$$alldata{$position}{'EVS_MAF_AA'} = $data[19];
		#$$alldata{$position}{'EVS_MAF_TAC'} = $data[20];
		#$$alldata{$position}{'EVS_mutation'} = $data[21];
		#$$alldata{$position}{'EVS_PolyPhen'} = $data[22];
		$$alldata{$position}{'dbSNP132'} = $data[19];
		$$alldata{$position}{'SIFT'} = $data[20];
		#$$alldata{$position}{'SIFT_Rank'} = 'NA';
		$$alldata{$position}{'PolyPhen2'} = $data[21];
		#$$alldata{$position}{'PolyPhen2_Rank'} = 'NA';
		$$alldata{$position}{'LJB_MutationTaster'} = $data[22];
		#$$alldata{$position}{'LJB_MutationTaster_Rank'} = 'NA';
		#$$alldata{$position}{'MutationAssessor'} = $data[27];
		$$alldata{$position}{'LJB_LRT'} = $data[23];
		#$$alldata{$position}{'LJB_LRT_Rank'} = 'NA';
		$$alldata{$position}{'CondelValue'} = 'NA';
		$$alldata{$position}{'CondelJudgement'} = 'NA',
		#$$alldata{$position}{'Condel_bin'} = 'NA',
		#$$alldata{$position}{'RankProduct'} = 'NA';
		$$alldata{$position}{'Ref'} = $data[24];
		$$alldata{$position}{'Obs'} = $data[25];
		$$alldata{$position}{'Qual_median'} = $data[26];
		$$alldata{$position}{'SampleDetails'} = $data[27];
		$$alldata{$position}{'Conservation_Rank_Used'} = 1000;
		
		
		$$alldata{$position}{'Freq'} = determineFreq($$alldata{$position}{'PositiveSamples'}, $$alldata{$position}{'NegativeSamplesSNPCall'});
		determineGeneFrequencies($$alldata{$position}{'PositiveSamples'}, $$alldata{$position}{'AAChange'}, $allGENEdata, $funcfilter, $exonicfuncfilter, $$alldata{$position}{'ExonicFunc'}, $$alldata{$position}{'Func'}); # count gene occurences
		
		#print Dumper($alldata); #DEBUG
		#exit; #DEBUG
	}
};


sub readSIFT{
	my $samplefiles = shift;
	my $alldata = shift;
	
	my $siftfolder = File::Spec->rel2abs($$samplefiles{'allsamples'}{'siftoutput'});
	my $predictions = $$samplefiles{'allsamples'}{'JobID'}."_predictions.tsv";
	
	#print "$siftfolder/$predictions", "\n"; #DEBUG
	
	open(siftFH, "$siftfolder/$predictions") or die $!;
	
	while (<siftFH>) {
		
		#print $_, "\n"; #DEBUG
		
		if ($_ =~ /^Coordinates/) {
			next;
		}
		my @splitline = split(/\t/, $_);
		my @coordinates = split(/,/, $splitline[0]);
		my $position = encryptposition($coordinates[0], $coordinates[1]);
		
		if ($splitline[0] =~ /Warning/) {
			$$alldata{$position}{'ref_error'} = 'yes'; # Warning: NCBI reference miscall!
		}
		else {
			$$alldata{$position}{'ref_error'} = 'no';
		}
		
		if ($splitline[3] eq '') {
			$$alldata{$position}{'type'} = $splitline[4]; # intron
		}
		else {
			$$alldata{$position}{'type'} = $splitline[3]; # aa_insertion, frameshift
		}
		
		# effect
		if ($splitline[6] eq '') {
			$$alldata{$position}{'effect'} = 'NA';
		}
		else {
			$$alldata{$position}{'effect'} = $splitline[6];
		}
		
		# aa change
		if ($splitline[10] eq '') {
			$$alldata{$position}{'siftAAchange'} = 'NA';
		}
		else {
			$$alldata{$position}{'siftAAchange'} = $splitline[10];
		}
		
		# location in transcript
		if ($splitline[11] eq '') {
			$$alldata{$position}{'location'} = 'NA';
		}
		else {
			$$alldata{$position}{'location'} = $splitline[11];
		}
		
		# causes NMD
		if ($splitline[12] eq '') {
			$$alldata{$position}{'nmd'} = 'NA';
		}
		else {
			$$alldata{$position}{'nmd'} = $splitline[12];
		}
		
		# repeats?
		if ($splitline[13] eq '') {
			$$alldata{$position}{'repeat'} = 'NA';
		}
		else {
			$$alldata{$position}{'repeat'} = $splitline[13];
		}
		
		# OMIM
		if ($splitline[23] eq '') {
			$$alldata{$position}{'OMIM'} = 'NA';
		}
		else {
			$$alldata{$position}{'OMIM'} = $splitline[23];
			$$alldata{$position}{'OMIM'} =~ s/,/;/g;
		}
		
		#print "position: $position", "\n";
		#print "ref_error: $$alldata{$position}{'ref_error'}", "\n";
		#print "type: $$alldata{$position}{'type'}", "\n";
		#print "effect: $$alldata{$position}{'effect'}", "\n";
		#print "siftAAchange: $$alldata{$position}{'siftAAchange'}", "\n";
		#print "location: $$alldata{$position}{'location'}", "\n";
		#print "nmd: $$alldata{$position}{'nmd'}", "\n";
		#print "repeat: $$alldata{$position}{'repeat'}", "\n";
		#print "OMIM: $$alldata{$position}{'OMIM'}", "\n";
		
		
		## 1000 genomes
		#if ($splitline[24] eq '') {
		#	$$alldata{$position}{'1000genomes'} = 'NA';
		#}
		#else {
		#	$$alldata{$position}{'1000genomes'} = $splitline[24];
		#}
	}
};

sub reporterrors {
	my $errormessages = shift;
	my $message = shift;
	if ($message) {
		push(@$errormessages, $message);
	}
	else {
		if (scalar(@$errormessages) >= 1) {
			print STDERR "The script did not run without errors: \n";
			foreach my $line(@$errormessages) {
				print STDERR $line, "\n";
			}
		}
	}
}

sub translate {
	my $Chr = shift;
	return 23 if ($Chr =~ /x/i);
	return 24 if ($Chr =~ /y/i);
	return 25 if ($Chr =~ /mt/i);
	return $Chr;
}

sub unzipfile {
	my $sample = shift;
	my $datatype = shift;
	my $samplefiles = shift;
	my $yesall = shift;
	my $unzipped = $$samplefiles{$sample}{$datatype}{'name'}."unzipped";
	
	if (-e $unzipped && $allyes == 1) { # $allyes does not work with threads
		if (usefile($unzipped, $yesall)) {
			print "using $unzipped \n"; #MESSAGE
		}
		else {
			print "inflating $sample $datatype data. this takes normally about 1,5min.\n"; #MESSAGE
			system("gunzip -c $$samplefiles{$sample}{$datatype}{'name'} > $unzipped"); #REAL
			print "inflating $sample $datatype data. Done.\n"; #MESSAGE
		}
	}
	else {
		print "inflating $sample $datatype data. this takes normally about 1,5min.\n"; #MESSAGE
		system("gunzip -c $$samplefiles{$sample}{$datatype}{'name'} > $unzipped"); #REAL
		print "inflating $sample $datatype data. Done.\n"; #MESSAGE
	}
	return($unzipped);
}

sub usage {
	print   "usage SNP_comparison.pl --dirlistinfile [file] --outfile \n",
		" or SNP_comparison.pl --collectionfile [file] \n",
		"invoke with required arguments: \n",
		"--dirlistinfile provides a list of directories where to find the files to investigate INSTEAD of --infile and --listinfile\n\tEach directory should be given to depth where to find \'Indel_Intersection\', \'SNP_Intersection\' and \'shore\'\n",
		"--collectionfile use a file that was previously built by SNP_comparison\n",
		"--useexclusionlist use a filter list that contains genes that frequently show up in ExomeSeq experiments \n\t[Analysis of DNA sequence variants detected by high-throughput sequencing; DOI: 10.1002/humu.22035]\n",
		"--evsfile location of the ExomeVariantServer data. [default = /users/GD/projects/genome_indices/human/hg19/EVS/ESP5400.chrall.snps.vcf]\n",
		"--outfile specifies the name of the result file\n", 
		"--mailaddress is required for SIFT usage\n", 
		"optional arguments: \n",
		#"--delimiter may be set if not tab, must be a unique letter in your file. [default ',']\n", 
		#"--frequencycutoff report only SNPs that are found in at least as many files as cutoff [default = 0].\n", 
		#"--genecutoff report only GENEs that are found in at least as many files as cutoff [default = 0]. The positions of SNPs doesn\'t affect this cutoff.\n", 
		#"--MAFcutoff take only SNPs into account that occur less than cutoff [default = 1.00].\n", 
		#"--FUNCfilter excludes all entries of the \'Func\' field, e.g. exonic or splicing. Is applied in combination (AND/&& linkage) with --EXONICFUNCfilter.\n", 
		#"--EXONICFUNCfilter excludes all entries of the \'ExonicFunc\' field, e.g. synonymous SNV. Is applied in combination (AND/&& linkage) with --FUNCfilter.\n",
		"--debug print more details on each SNP. Like ranking for each parameter in rank product.",
		"--nodescription print no description of each SNP once.\n", 
		"--allyes reply all questions regarding of reusing unzipped files with 'yes' [default = 'no']\n", 
		"get help: \n",
		"--help show this usage information\n",
		"Or: Call Oliver for direct help.\n\n";
	exit;
}

sub usefile {
	my $file = shift;
	my $yesall = shift;
	
	return(1) if ($yesall);
	
	my $yesno = 'x'; #put back to x
	while ($yesno ne 'y' || $yesno ne 'n') {
		print "use already unzipped file: $file? (y/n)";
		$yesno = <STDIN>;#REAL
		chomp($yesno);
		if ($yesno eq 'y') {
			return(1);
			last;
		}
		elsif ($yesno eq 'n') {
			return(0);
			last;
		}
	}
}

sub writeCondelInput{
	my $samplefiles = shift;
	my $position = shift;
	my $sift = shift;
	my $logre = shift;
	my $polyphen = shift;
	my $mapp = shift;
	my $mutass = shift;
	
	my $fh = $$samplefiles{'allsamples'}{'condelfilehandle'};
	print $fh join ("\t", $position, $sift, $logre, $polyphen, $mapp, $mutass, "\n");
}

sub writeFinalResults {
	my $samplefiles = shift;
	my $alldata = shift;
	my $delimiter = shift;
	my $debug = shift;
	my $allGENEdata = shift;
	my $genecutoff = shift;
	my $excludedgenes = shift;
	my $usecollectionfile = shift;
	
	my $outFH = $$samplefiles{'allsamples'}{'outfilehandle'};
	
	writeFinalHeader($samplefiles, $delimiter, $debug, 'outfilehandle');
	
	foreach my $position (keys %{$alldata}) {
		my $genefrequency;
		next unless ($genefrequency = filterGeneFrequency($alldata{$position}{'AAChange'}, $allGENEdata, $genecutoff));
		my ($chromosome, $nucleotide) = decifferposition($position);
		
		if (!defined($$alldata{$position}{'ref_error'})){
			$$alldata{$position}{'ref_error'} = 'NA'
		}
		if (!defined($$alldata{$position}{'type'})){
			$$alldata{$position}{'type'} = 'NA'
		}
		if (!defined($$alldata{$position}{'effect'})){
			$$alldata{$position}{'effect'} = 'NA'
		}
		if (!defined($$alldata{$position}{'siftAAchange'})){
			$$alldata{$position}{'siftAAchange'} = 'NA'
		}
		if (!defined($$alldata{$position}{'nmd'})){
			$$alldata{$position}{'nmd'} = 'NA'
		}
		if (!defined($$alldata{$position}{'OMIM'})){
			$$alldata{$position}{'OMIM'} = 'NA'
		}
		
		my $complete = join($delimiter,
			       $chromosome,
			       $nucleotide,
			       $$alldata{$position}{'Ref'},
			       $$alldata{$position}{'Obs'},
			       "0",#$$alldata{$position}{'RankProduct_Sorted'}, #DEBUG - REAL
			       $$alldata{$position}{'PositiveSamples'},
			       $$alldata{$position}{'NegativeSamplesSNPCall'},
			       $$alldata{$position}{'NegativeSamplesRefCall'},
                               $$alldata{$position}{'NegativeSamplesRefCall_lowqual'},
			       $$alldata{$position}{'AverageCoverage'},
			       $$alldata{$position}{'SamplesCoverage'},
			       $$alldata{$position}{'AllelFreq'},
			       $$alldata{$position}{'Freq'},
			       $$alldata{$position}{'Gene'},
			       $$alldata{$position}{'Conserved'},
			       $$alldata{$position}{'SegDup'},
			       $$alldata{$position}{'dbSNP132'},
			       $$alldata{$position}{'Func'},
			       $$alldata{$position}{'ExonicFunc'},
			       $$alldata{$position}{'AAChange'},
			       $$alldata{$position}{'ref_error'},
			       $$alldata{$position}{'type'},
			       $$alldata{$position}{'effect'},
			       $$alldata{$position}{'siftAAchange'},
			       $$alldata{$position}{'nmd'},
			       $$alldata{$position}{'OMIM'},
			       $$alldata{$position}{'Qual_median'},
			       $$alldata{$position}{'SampleDetails'},
			       "\n");
		
		my $shortened = join($delimiter,
			       $chromosome,
			       $nucleotide,
			       $$alldata{$position}{'Ref'},
			       $$alldata{$position}{'Obs'},
			       "0",#$$alldata{$position}{'RankProduct_Sorted'}, #DEBUG - REAL
			       $$alldata{$position}{'PositiveSamples'},
			       $$alldata{$position}{'NegativeSamplesSNPCall'},
			       $$alldata{$position}{'NegativeSamplesRefCall'},
                               $$alldata{$position}{'NegativeSamplesRefCall_lowqual'},
			       $$alldata{$position}{'AverageCoverage'},
			       $$alldata{$position}{'SamplesCoverage'},
			       $$alldata{$position}{'AllelFreq'},
			       $$alldata{$position}{'Freq'},
			       $$alldata{$position}{'Gene'},
			       $$alldata{$position}{'Conserved'},
			       $$alldata{$position}{'SegDup'},
			       $$alldata{$position}{'dbSNP132'},
			       $$alldata{$position}{'Func'},
			       $$alldata{$position}{'ExonicFunc'},
			       $$alldata{$position}{'AAChange'},
			       $$alldata{$position}{'ref_error'},
			       $$alldata{$position}{'type'},
			       $$alldata{$position}{'effect'},
			       $$alldata{$position}{'siftAAchange'},
			       $$alldata{$position}{'nmd'},
			       $$alldata{$position}{'OMIM'},
			       $$alldata{$position}{'Qual_median'},
			       $$alldata{$position}{'SampleDetails'},
			       "\n"); # 'Chr', 'Pos', 'Ref', 'Obs', 'Rank', 'SNPcall:highQual', 'SNPcall:lowQual', 'Refcall:hom', 'Refcall:het', 'AverageCoverage', 'SamplesCoverage', 'AllelFreq', 'Freq', 'Gene', 'Conserved', 'PhyloP', 'SegDup', 'MAF_1000G', 'EVS_MAF_EA', 'EVS_MAF_AA', 'EVS_MAF_TAC', 'dbSNP132', 'SIFT', 'PolyPhen2', 'MutationAssessor', 'MutationTaster', 'LRT', 'Condel', 'CondelJudgement', 'EVS_PolyPhen', 'Function', 'ExonicFunc', 'EVS_exonic_function', 'AAChange',  'Qual_median', 'SampleDetails(SNP[0|1|2]Ref[0|1|2])'
		
		
		if ($useexclusionlist == 0 && $debug == 1) { # if exclusion list is not to be used print out directly
			print $outFH $complete;
		}
		elsif ($useexclusionlist == 0 && $debug == 0) {
			print $outFH $shortened;
		}
		elsif ($useexclusionlist == 1 && $debug == 1) { # if exclusion list is going to be used, test if 'Gene' is element of excludedgenes
			if (isExcludedGene($$alldata{$position}{'Gene'}, $excludedgenes) == 0) {
				print $outFH $complete;
			}
		}
		elsif ($useexclusionlist == 1 && $debug == 0) {
			if (isExcludedGene($$alldata{$position}{'Gene'}, $excludedgenes) == 0) {
				print $outFH $shortened;
			}
		}
		else {
			#don't do anything
		}
	}
	

}

sub writeFinalHeader {
	my $samplefiles = shift;
	my $delimiter = shift;
	my $debug = shift;
	my $filehandle = shift;
	my $fh = $$samplefiles{'allsamples'}{$filehandle}; # accessing the FileHandle directly doesn't work
	#print $fh join($delimiter, 'Chr', 'Pos', 'PositiveSamples', 'NegativeSamplesRefCall_hom', 'NegativeSamplesRefCall_het', 'NegativeSamplesSNPCall', 'SamplesCoverage', 'AllelFreq', 'Freq', 'Func', 'Gene', 'ExonicFunc', 'AAChange', 'Conserved', 'LJB_PhyloP', 'ConservationValue_Rank', 'ConservationTool', 'SegDup', 'SegDup_Rank', 'MAF_1000G', 'MAF_Rank', 'EVS_MAF_EA', 'EVS_MAF_AA', 'EVS_MAF_TAC', 'EVS_mutation', 'EVS_PolyPhen', 'dbSNP132', 'SIFT', 'PolyPhen2', 'LJB_MutationTaster', 'LJB_LRT', 'CondelValue', 'CondelJudgement', 'Condel_Rank', 'RankProduct', 'Ref', 'Obs',  'Qual_median', 'SampleDetails(SNP[0|1|2]Ref[0|1|2])', "\n"); # , 'SIFT_Rank', 'PolyPhen2_Rank', 'LJB_MutationTaster_Rank', 'LJB_LRT_Rank' 'GenFreq',
	if ($debug == 1) {
		print $fh join($delimiter, 'Chr', 'Pos', 'Ref', 'Obs', 'Rank', 'SNPcall:highQual', 'SNPcall:lowQual', 'Refcall:hom', 'Refcall:het', 'AverageCoverage', 'SamplesCoverage', 'AllelFreq', 'Freq', 'Gene', 'Conserved', 'SegDup', 'dbSNP132',  'Function', 'ExonicFunc', 'AAChange', 'SIFT_ReferenceError', 'SIFT_Function', 'SIFT_impact', 'SIFT_AAChange', 'SIFT_NMD', 'SIFT_OMIM', 'Qual_median', 'SampleDetails(SNP[0|1|2]Ref[0|1|2])', "\n"); #, 'EVS_MAF_EA', 'EVS_MAF_AA', 'EVS_MAF_TAC', 'MutationAssessor', 'EVS_exonic_function', 'EVS_PolyPhen''ConservationValue_Rank',, 'PhyloP' 'ConservationTool', , 'MAF_1000G', 'MAF_Rank''SIFT', 'PolyPhen2', 'MutationTaster', 'LRT', 'Condel', 'CondelJudgement', 'Condel_Rank',, 'SegDup_Rank'
	}
	elsif ($debug == 0) {
		print $fh join($delimiter, 'Chr', 'Pos', 'Ref', 'Obs', 'Rank', 'SNPcall:highQual', 'SNPcall:lowQual', 'Refcall:hom', 'Refcall:het', 'AverageCoverage', 'SamplesCoverage', 'AllelFreq', 'Freq', 'Gene', 'Conserved', 'SegDup', 'dbSNP132', 'Function', 'ExonicFunc', 'AAChange', 'SIFT_ReferenceError', 'SIFT_Function', 'SIFT_impact', 'SIFT_AAChange', 'SIFT_NMD', 'SIFT_OMIM',  'Qual_median', 'SampleDetails(SNP[0|1|2]Ref[0|1|2])', "\n"); #'MAF_1000G',, 'SIFT', 'PolyPhen2', 'MutationAssessor', 'MutationTaster', 'LRT', 'Condel', 'CondelJudgement', 'EVS_MAF_EA', 'EVS_MAF_AA', 'EVS_MAF_TAC', 'EVS_exonic_function', 'EVS_PolyPhen', 'PhyloP'
	}
}

sub writeHeader {
	my $samplefiles = shift;
	my $delimiter = shift;
	my $filehandle = shift;
	my $fh = $$samplefiles{'allsamples'}{$filehandle}; # accessing the FileHandle directly doesn't work
	print $fh join($delimiter, 'Chr', 'Pos', 'EndPos', 'PositiveSamples', 'NegativeSamplesRefCall_hom', 'NegativeSamplesRefcall_het', 'NegativeSamplesSNPCall', 'AverageCoverage', 'SamplesCoverage', 'AllelFreq', 'Freq', 'Func', 'Gene', 'ExonicFunc', 'AAChange', 'Conserved', 'LJB_PhyloP', 'SegDup', 'MAF_1000G', 'dbSNP132', 'SIFT', 'PolyPhen2', 'LJB_MutationTaster', 'LJB_LRT', 'Ref', 'Obs', 'Qual_median', 'SampleDetails(SNP[0|1|2]Ref[0|1|2])', "\n"); #'EVS_MAF_EA', 'EVS_MAF_AA', 'EVS_MAF_TAC', 'EVS_mutation', 'EVS_PolyPhen', ... , 'MutationAssessor'
}
sub writeResults {
	my $samplefiles = shift;
	my $delimiter = shift;
	my $frequencycutoff = shift;
	my $genecutoff = shift;
	my $MAFcutoff = shift;
	my $funcfilter = shift;
	my $exonicfuncfilter = shift;
	
	my @positiveSamples;
	my @negativeSamples;
	my %printout; 
	
	# identify positive and negative samples	
	foreach my $sample (keys %{$samplefiles}) {
		next if ($sample eq 'allsamples'); # skip the general parameter set
		if ($$samplefiles{$sample}{'data_pos'} == 1) {
			push (@positiveSamples, $sample);
		}
		elsif ($$samplefiles{$sample}{'data_neg'} == 1) {
			push (@negativeSamples, $sample);
		}
	}

	# enter EVS data to printout array
	my ($chrposition, $nucleotideposition) = decifferposition($$samplefiles{'allsamples'}{'currentposition'});
	$printout{'chromosome'} = $chrposition;
	$printout{'nucleotide'} = $nucleotideposition;
	$printout{'qualitymedian'} = $$samplefiles{'allsamples'}{'quality'};
	$printout{'coverage_average'} = $$samplefiles{'allsamples'}{'coverage_average'};
	
	# enter allele frequencies
	$printout{'AllelFreq'} = $$samplefiles{'allsamples'}{'AllelFreq'};
	$printout{'Freq'} = $$samplefiles{'allsamples'}{'Freq'};
	## enter mutationAssessor
	#$printout{'mutass'} = $$samplefiles{'allsamples'}{'mutass'};
	
	
	foreach my $positivesample (@positiveSamples) {
		next if ($$samplefiles{$positivesample}{'Indel'}{'lastline'} eq 'eof');
		# compile all positive samples
		if ($printout{'positivesamples'}) {
			$printout{'positivesamples'} = $printout{'positivesamples'}.' '.$positivesample;
		}
		else {
			$printout{'positivesamples'} = $positivesample;
		}
		# add coverage to sample name
		if ($printout{'samplecoverage'}) {
			$printout{'samplecoverage'} = $printout{'samplecoverage'}.' '.$positivesample.'_'.$$samplefiles{$positivesample}{'coverage'}{'value'};
		}
		else {
			$printout{'samplecoverage'} = $positivesample.'_'.$$samplefiles{$positivesample}{'coverage'}{'value'};
		}
		
		# add sample details
		# fill in hom/het info for ref and SNP position
		my $ref = 0;
		my $snp = 0;
		
		# first look at ref position
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$positivesample}{'reference'}{'foundposition'} && $$samplefiles{$positivesample}{'reference'}{'nucleotide'} ne 'N') {
			if ($$samplefiles{$positivesample}{'reference'}{'state'} eq 'hom') {
				$ref = 2;
			}
			elsif ($$samplefiles{$positivesample}{'reference'}{'state'} eq 'het') {
				$ref = 1;
			}
			else {
				$ref = 0;
			}
		}
		# now check out the SNP position
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$positivesample}{'allIndel'}{'foundposition'}) {
			my @SNPtools = ('GATK', 'MPILEUP', 'SHORE');
			foreach my $tool (@SNPtools) { 
				if ($$samplefiles{$positivesample}{'allIndel'}{$tool} =~ /^\.\/\./ || $$samplefiles{$positivesample}{'allIndel'}{$tool} eq 'NA') {
					$snp = 0;
					next;
				}
				elsif ($$samplefiles{$positivesample}{'allIndel'}{$tool} =~ /^0\/1/) {
					$snp = 1;
					
					#adjust hom ref, if a heterozygous indel is seen
					if ($$samplefiles{$positivesample}{'reference'}{'state'} eq 'hom') {
						# if ref is alredy 0 do not decrease
						if ($ref != 0) {$ref -= 1;};
					}
					
					last; # use the tools in the order GATK, MPILEUP, SHORE - first positive hit found will be used
				}
				elsif ($$samplefiles{$positivesample}{'allIndel'}{$tool} =~ /^1\/1/) {
					$snp = 2;
					#adjust hom ref and indel, if a both indel and ref are seen - obviously the algorithms confuse hom/het with indels...
					if ($$samplefiles{$positivesample}{'reference'}{'state'} eq 'hom') {
						$ref -= 1;
						$snp -= 1;
					}
					
					last; # use the tools in the order GATK, MPILEUP, SHORE - first positive hit found will be used
				}
			}
		}
		# enter the name 
		if ($printout{'sampledetails'}) { 
			$printout{'sampledetails'} = $printout{'sampledetails'}.' '.$positivesample.'_'.$snp.$ref;
			}
		else {
			$printout{'sampledetails'} = $positivesample.'_'.$snp.$ref;
		}
		
		unless ($printout{'func'}) { #test if the stuff was seen before
			my @annovardata = split($delimiter, $$samplefiles{$positivesample}{'Indel'}{'lastline'});
			unless($annovardata[23]) {
				$annovardata[23] = 'NA';
			}
			$printout{'func'} = clean($annovardata[0]);
			$printout{'gene'} = clean($annovardata[1]);
			$printout{'exonicfunc'} = clean($annovardata[2]);
			$printout{'aachange'} = clean($annovardata[3]);
			$annovardata[4] = clean($annovardata[4]);
			my @conserveddata = split(/\;/, $annovardata[4]); # see below
			$printout{'conserved'} = $conserveddata[0];
			$printout{'segdup'} = clean($annovardata[5]);
			$printout{'maf_1000g'} = clean($annovardata[6]);
			$printout{'dbsnp'} = clean($annovardata[9]);
			$printout{'sift'} = clean($annovardata[10]);
			$printout{'polyphen'} = clean($annovardata[11]);
			$printout{'ljb_phylop'} = clean($annovardata[12]);
			$printout{'ljb_mutationtaster'} = clean($annovardata[13]);
			$printout{'ljb_lrt'} = clean($annovardata[14]);
			$printout{'endposition'} = clean($annovardata[17]);
			$printout{'ref'} = clean($annovardata[18]);
			$printout{'obs'} = clean($annovardata[19]);
			$printout{'other'} = clean($annovardata[20]);
			$printout{'qual'} = clean($annovardata[21]);
			$printout{'depth'} = clean($annovardata[22]);
			$printout{'qualitybydepth'} = clean($annovardata[23]);
			chomp($printout{'qualitybydepth'});
			# conserved data in Annovar:
			# http://www.openbioinformatics.org/annovar/annovar_region.html
			# The output is saved in the varlist.hg18_phastConsElements44way file. The first column in the output is "mce44way" indicating the type of annotation. The second column is the normalized score assigned by UCSC Genome Browser, and this score range from 0 to 1000. (Note that the --score_threshold or --normscore_threshold can also be used to filter out specific variants with low conservation scores.) The second column also shows "Name=lod=x" which is used to tell the user the LOD score for the region. All other columns are identical as those in the input file. Only variants that actually are located within a conserved region will be printed in the output file. As a result, only 5 variants are in the output file. 
		}
	}
	
	
	
	foreach my $negativesample (@negativeSamples) {
		# fill in the coverage for the    position if there's some at this position
		if ($$samplefiles{$negativesample}{'coverage'}{'value'} && $$samplefiles{$negativesample}{'coverage'}{'value'} ne 'NA') { 
                        if ($printout{'samplecoverage'}) {
                                $printout{'samplecoverage'} = $printout{'samplecoverage'}.' '.$negativesample.'_'.$$samplefiles{$negativesample}{'coverage'}{'value'};
                                }
                        else {
                                $printout{'samplecoverage'} = $negativesample.'_'.$$samplefiles{$negativesample}{'coverage'}{'value'};
                        }
                }		

		# fill in a low quality SNP if there's one for this position
		my $negativeforprintout = $negativesample;		
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$negativesample}{'allIndel'}{'foundposition'}) { 
			$negativeforprintout = preparenegativename($negativesample,
								   $$samplefiles{$negativesample}{'allIndel'}{'GATK'},
								   $$samplefiles{$negativesample}{'allIndel'}{'MPILEUP'},
								   $$samplefiles{$negativesample}{'allIndel'}{'SHORE'});
			if ($printout{'negativelowqualsnp'}) {
				$printout{'negativelowqualsnp'} = $printout{'negativelowqualsnp'}.' '.$negativeforprintout;
			}
			else {
				$printout{'negativelowqualsnp'} = $negativeforprintout;
			}
		}
		
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$negativesample}{'reference'}{'foundposition'} && $$samplefiles{$negativesample}{'reference'}{'nucleotide'} ne 'N') {
			if ($$samplefiles{$negativesample}{'reference'}{'state'} eq 'hom' && $$samplefiles{'allsamples'}{'currentposition'} != $$samplefiles{$negativesample}{'allIndel'}{'foundposition'})  { # check if we have a reference call or if we have a low quality SNP, which would have been covered in the block above
				if ($printout{'samplereference'}) { 
					$printout{'samplereference'} = $printout{'samplereference'}.' '.$negativesample;
					}
				else {
					$printout{'samplereference'} = $negativesample;
				}
			}
			elsif ($$samplefiles{$negativesample}{'reference'}{'state'} eq 'het' && $$samplefiles{'allsamples'}{'currentposition'} != $$samplefiles{$negativesample}{'allIndel'}{'foundposition'}) { # check if we have a low support reference call or if we have a low quality SNP, which would have been covered in the block above
				if ($printout{'samplereference_lowqual'}) {
                                        $printout{'samplereference_lowqual'} = $printout{'samplereference_lowqual'}.' '.$negativesample;
                                        }
                                else {
                                        $printout{'samplereference_lowqual'} = $negativesample;
                                }
			}
		}
		
		# fill in hom/het info for ref and SNP position
		my $ref = 0;
		my $snp = 0;
		
		# first look at ref position
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$negativesample}{'reference'}{'foundposition'} && $$samplefiles{$negativesample}{'reference'}{'nucleotide'} ne 'N') {
			if ($$samplefiles{$negativesample}{'reference'}{'state'} eq 'hom') {
				$ref = 2;
			}
			elsif ($$samplefiles{$negativesample}{'reference'}{'state'} eq 'het') {
				$ref = 1;
			}
			else {
				$ref = 0;
			}
		}
		# now check out the SNP position
		if ($$samplefiles{'allsamples'}{'currentposition'} == $$samplefiles{$negativesample}{'allIndel'}{'foundposition'}) {
			my @SNPtools = ('GATK', 'MPILEUP', 'SHORE');
			foreach my $tool (@SNPtools) { 
				if ($$samplefiles{$negativesample}{'allIndel'}{$tool} =~ /^\.\/\./ || $$samplefiles{$negativesample}{'allIndel'}{$tool} eq 'NA') {
					$snp = 0;
					next;
				}
				elsif ($$samplefiles{$negativesample}{'allIndel'}{$tool} =~ /^0\/1/) {
					$snp = 1;
					last; # use the tools in the order GATK, MPILEUP, SHORE - first positive hit found will be used
				}
				elsif ($$samplefiles{$negativesample}{'allIndel'}{$tool} =~ /^1\/1/) {
					$snp = 2;
					last; # use the tools in the order GATK, MPILEUP, SHORE - first positive hit found will be used
				}
			}
		}
		# enter the name 
		if ($printout{'sampledetails'}) { 
			$printout{'sampledetails'} = $printout{'sampledetails'}.' '.$negativesample.'_'.$snp.$ref;
			}
		else {
			$printout{'sampledetails'} = $negativesample.'_'.$snp.$ref;
		} 
	}
	
	# fill fields with standard values if empty (including pos/neg/ref samples)
	fillEmptyValues(\%printout);
	fillSampleFields(\%printout);
	
	# finally print to file
	my $fh = $$samplefiles{'allsamples'}{'collectionfilehandle'}; # accessing the FileHandle driectly doesn't work
	print $fh join($delimiter,
	       $printout{'chromosome'},
	       $printout{'nucleotide'},
	       $printout{'endposition'},
	       $printout{'positivesamples'},
	       $printout{'samplereference'},
	       $printout{'samplereference_lowqual'},
	       $printout{'negativelowqualsnp'},
	       $printout{'coverage_average'},
	       $printout{'samplecoverage'},
	       $printout{'AllelFreq'},
	       $printout{'Freq'},
	       $printout{'func'},
	       $printout{'gene'},
	       $printout{'exonicfunc'},
	       $printout{'aachange'},
	       $printout{'conserved'},
	       $printout{'ljb_phylop'},
	       $printout{'segdup'},
	       $printout{'maf_1000g'},
	       $printout{'dbsnp'},
	       $printout{'sift'},
	       $printout{'polyphen'},
	       $printout{'ljb_mutationtaster'},
	       $printout{'ljb_lrt'},
	       $printout{'ref'},
	       $printout{'obs'},
	       $printout{'qualitymedian'},
	       $printout{'sampledetails'},
	       "\n");
}

sub writeSIFTconfig{
	my $samplefiles = shift;
	my $siftconfFH = $$samplefiles{'allsamples'}{'siftconffilehandle'};
	my $infile = shift;
	#my $infile_sift = $$samplefiles{'allsamples'}{'siftinput'};
	my $infile_sift = $$samplefiles{'allsamples'}{'siftfilelocation'};
	my $outfile = $$samplefiles{'allsamples'}{'siftoutput'};
	my $fulloutfile = File::Spec->rel2abs( $outfile );
	my $mail = $$samplefiles{'allsamples'}{'mailaddress'};
	
	mkdir($outfile, 0700);
	
	print $siftconfFH join("\n", 	"ORGANISM='hs37'",	
					"FILE1='$infile_sift'	#path to input file ",
					"EMAIL='$mail'	#fill in your email address",
					"FINALOUTDIR='$fulloutfile'	#location to download the output file",
					
					"#Variables below take value of 0 or 1; 0 = OFF, 1 = ON",
					"ALL_TRANSCRIPTS=1",
					#"EBL_GENE_ID=1",
					"GENE_NAME=1",
					"GENE_DESC=1",
					"EBL_TRANS_STATUS=1",
					"EBL_PROTEIN_FAMILY_ID=1",
					"EBL_PROTEIN_FAMILY_ID_DESC=1",
					"PROTEIN_FAMILY_SIZE=1",
					"KA_KS_HUMAN_MOUSE=1",
					"KA_KS_HUMAN_MACAQUE=1",
					"OMIM_DISEASE=1",
					"OVERLAP_WITH_1000_GENOME_INDEL=1",
					"\n");
					#"ALL_WEIGHTED_AVERAGE=1",
					#"CEU_POPULATION=1",
					#"HCB_POPULATION=1",
					#"JPT_POPULATION=1",
					#"YRI_POPULATION=1",
					#"ALL_AVAIL_POPULATIONS=1",
					#"EUROPEAN_POPULATION=1",
					#"EAST_ASIAN_POPULATION=1",
					#"AFRICAN_POPULATION=1",
					#"SOUTH_ASIAN_POPULATION=1",
					#"AMERICAN_POPULATION=1");
}

sub writeSIFTinput{
	my $samplefiles = shift;
	my $alldata = shift;
	my $siftFH = $$samplefiles{'allsamples'}{'siftfilehandle'};
	my $counter = 0;
	
	#insertion or deletion?
	foreach my $position (sort(keys %{$alldata})) { # sorted, because then it would later on be easier to fuse results
		my ($chromosome, $nucleotide) = decifferposition($position);
		if (length($$alldata{$position}{'Ref'}) > length($$alldata{$position}{'Obs'})) {
			# that's a deletion
			my $endposition = $nucleotide + length($$alldata{$position}{'Ref'});
			print $siftFH join(',', $chromosome, $nucleotide, $endposition, '1', "$$alldata{$position}{'Ref'}\/-");
			print $siftFH "\n";
		}
		else {
			#that's an insertion
			print $siftFH join(',', $chromosome, $nucleotide, $nucleotide, '1', "-/$$alldata{$position}{'Obs'}");
			print $siftFH "\n";
		}
	}
}
__END__
