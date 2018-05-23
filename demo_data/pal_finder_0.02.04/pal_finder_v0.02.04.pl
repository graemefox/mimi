#!/usr/bin/perl

##########################################################################
#                                                                        #
# File: pal_finder_v0.02.04                                              #
# Authors: Alex Poole                                                    #
# Version: 0.02.04                                                       #
#                                                                        #
##########################################################################

=license

	Copyright (c) 2010, 2011, 2012 Alex Poole
	
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;

######################################### DECLARING GLOBAL VARIABLES ###############################################

# Reading control information
my $globals = {};
my $programName = "pal_finder_v0.02.04.pl";

############################################### MAIN PROGRAM #######################################################

# Read input parameters from configuration file
my $configFile = shift(@ARGV);
die "$programName: No control file given\nUsage: perl $programName <config_file>\n" unless $configFile;

ReadConfigFile($configFile, $globals);
createPr3Hash($globals);


my $findPrimers = $globals->{findPrimers};
print "Creating Microsat Hash...\n";
my $masterMicrosatHash = createMicrosatHashCount($globals);
my $microsatTally; # = createMicrosatTallyHash($globals);
my $palHash = {};

if ($globals->{platform} eq "454") {
	if( ($globals->{inputFormat} eq "fastq") or ($globals->{inputFormat} eq "scarf") ) { 
		die "\nCurrently fastq and scarf input formats work only with Illumina paired end data.\n";
	}
	print "Finding primers for 454 ";
	if ($globals->{pairedEnd} == 0) {
		print "single reads...\n";
		print "Scanning reads for microsatellites...\n";
		$microsatTally = createMicrosatTallyHash($globals);			
		if($findPrimers) { 
			open(BOULDER_IO, ">$globals->{primer3input}") or die "could not open file: $!\n";
			$microsatTally->{totalReads} = scan454reads($globals->{input454reads},
							    									   \&createBoulderIOfile,
							    									   $globals,
							    									   $palHash,
							    									   $masterMicrosatHash,
							      									   $microsatTally,
							    									   \*BOULDER_IO);
			close(BOULDER_IO);
			my $commandLine = "$globals->{primer3executable} < $globals->{primer3input} > $globals->{primer3output}";
			print "Finding primer pairs with primer3...\n";
			print `$commandLine`;
		}
		else { $microsatTally->{totalReads} = scan454reads($globals->{input454reads},
							    									   \&createBoulderIOfile,
							    									   $globals,
							    									   $palHash,
							    									   $masterMicrosatHash,
							      									   $microsatTally);
		}
		print "Writing results...\n";
		writeResults($globals, $microsatTally, $palHash);
		if( ($findPrimers) and (!($globals->{keepPrimer3files})) ) {
			my $delCommand = "rm $globals->{primer3input} $globals->{primer3output}";
			print `$delCommand`;
		}
	}
	else {
		print "paired end reads\nProgram does not currently support 454 paired end data\n";
		exit(1);
	}
}
elsif ($globals->{platform} eq "Illumina") {
	if( $globals->{inputFormat} eq "fasta" ) { 
		die "\nCurrently fasta input format works only with 454 single end read data.\n";
	}
	print "Finding primers for Illumina ";
	if($globals->{pairedEnd} == 1) {
		print "paired end reads...\n";
		print "Scanning reads for microsatellites...\n";
		$microsatTally = createPEmicrosatTallyHash($globals);
		if($findPrimers) {
			open(BOULDER_IO, ">$globals->{primer3input}") or die "could not open file: $!\n";
			if($globals->{inputFormat} eq "scarf") {			
				$microsatTally->{totalReads} = scanScarfPEreads($globals->{inputReadFile},
									   $globals->{pairedReadFile},
									   \&createPairEndBoulderIOfile,
									   $globals,
									   $masterMicrosatHash,
									   $microsatTally,
									   $palHash,
									   \*BOULDER_IO);
			}
			elsif($globals->{inputFormat} eq "fastq") {
				$microsatTally->{totalReads} = scanFastqPEreads($globals->{inputReadFile},
									   $globals->{pairedReadFile},
									   \&createPairEndBoulderIOfile,
									   $globals,
									   $masterMicrosatHash,
									   $microsatTally,
									   $palHash,
									   \*BOULDER_IO);
			}
			else { die "\nInvalid input file format $globals->{inputFormat}\n"; }
			close(BOULDER_IO);
			if ($microsatTally->{readsWithMicrosat}) {
				my $commandLine = "$globals->{primer3executable} < $globals->{primer3input} > $globals->{primer3output}";
				print "Finding primer pairs with primer3...\n";
				print `$commandLine`;
			}
			else {
				print "No microsatellites found in any reads. Ending script.\n";
				exit(0);
			}
		}
		else {
			if($globals->{inputFormat} eq "scarf") {			
				$microsatTally->{totalReads} = scanScarfPEreads($globals->{inputReadFile},
									   $globals->{pairedReadFile},
									   \&createPairEndBoulderIOfile,
									   $globals,
									   $masterMicrosatHash,
									   $microsatTally,
									   $palHash);
			}
			elsif($globals->{inputFormat} eq "fastq") {
				$microsatTally->{totalReads} = scanFastqPEreads($globals->{inputReadFile},
									   $globals->{pairedReadFile},
									   \&createPairEndBoulderIOfile,
									   $globals,
									   $masterMicrosatHash,
									   $microsatTally,
									   $palHash);
			}
			else { die "\nInvalid input file format $globals->{inputFormat}\n"; }
		}
		print "Writing results...\n";
		writePEresults($globals, $microsatTally, $palHash);
		if( ($findPrimers) and (!($globals->{keepPrimer3files})) ) {
			my $delCommand = "rm $globals->{primer3input} $globals->{primer3output}";
			print `$delCommand`;
		}
	}
	else {
		print "single end reads.\nprogram does not currently support Illumina single end data\n";
		exit(1);
	}
}
else {
	print "Invalid platform request $globals->{platform} \n";
	exit(1);
}

print "Done!!\n";

################################################# END MAIN #########################################################

sub getConfigFile {
	foreach my $arg (@_) {
		print $arg;
	}
}

sub createPr3Hash {
	my $globals = shift;
	
	my $rangeMin = "60";
	my $rangeMax = "800";
	my $pr3hash = {};
	foreach my $variable (keys(%$globals)) {
		my $wordsize = length($variable);
		if($variable =~ m/[A-Z_]{$wordsize}/) {
			$pr3hash->{$variable} = $globals->{$variable};
			delete $globals->{$variable};
		}
		elsif($variable eq "pr3ProductSizeRangeMinVal") { 
			$rangeMin = $globals->{$variable};
			delete $globals->{$variable};
		}
		elsif($variable eq "pr3ProductSizeRangeMaxVal") {
			$rangeMax = $globals->{$variable};
			delete $globals->{$variable};
		}
	}
	$pr3hash->{PRIMER_PRODUCT_SIZE_RANGE} = "$rangeMin-$rangeMax";
	$globals->{pr3hash} = $pr3hash;
}

sub createHeaderHash {
	my %headerHash = ();
	my $columnNum = scalar(@_);
	for (my $i=0; $i < $columnNum; $i++) { 
		my $title = shift(@_);
		$headerHash{$title} = $i;
	}
	return \%headerHash;
}

sub createBoulderIOfile {
	my $accno = shift;
	my $seq = shift;
	my $globs = shift;
	my $palHash = shift;
	my $microsatHash = shift;
	my $allRepeatTally = shift;
	my $param_ref = $globs->{pr3hash};
	my $findPrimers = $globs->{findPrimers};
	my $boulderIO_FH = "";
	if($findPrimers) { $boulderIO_FH = shift; }

	my @primer3List = ();
	my @monoLengths = keys(%$microsatHash);

	my $seqID = "";
	my $seqInfo = scanSeqAgainstHash($seq, $microsatHash, $globs, @monoLengths);
	$allRepeatTally->{totalBases}+= length($seq);
	updateTally($allRepeatTally, $seqInfo, $palHash, $accno);
	if (!($seqInfo)) { $seqID = $accno; }
	else {
		$allRepeatTally->{readsWithMicrosat}++;
		$seqID = "$accno" . "_" . $seqInfo;
		if($findPrimers) {
			$param_ref->{SEQUENCE_ID} = $seqID;
			updatePr3Info($seqInfo, $seq, $param_ref, $allRepeatTally);
			foreach my $parameter (keys(%$param_ref)) {
				print $boulderIO_FH "$parameter=$$param_ref{$parameter}\n";
			}
			print $boulderIO_FH "=\n";
		}
	}
	return \$allRepeatTally;
}

sub createPairEndBoulderIOfile {
	my $accno = shift; shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $qual1 = shift;
	my $qual2 = shift;
	my $globs = shift;
	my $microsatHash = shift;
	my $allRepeatTally = shift;
	my $PALhash = shift;
	my $param_ref = $globs->{pr3hash};
	my $findPrimers = $globs->{findPrimers};
	my $boulderIO_FH = "";
	$seq2 = reverseComplement($seq2);
	$qual2 = reverse($qual2);
	if($findPrimers) { $boulderIO_FH = shift; }

	my @primer3List = ();
	my @monoLengths = keys(%$microsatHash);

	my $seqID = ""; my $seqInfo = "";
	my $seq1Info = scanSeqAgainstHash($seq1, $microsatHash, $globs, @monoLengths);
	my $seq2Info = scanSeqAgainstHash($seq2, $microsatHash, $globs, @monoLengths);
	$allRepeatTally->{totalBases}+= length($seq1);
	$allRepeatTally->{totalBases}+= length($seq2);
	$seqInfo = "($seq1Info)($seq2Info)";
	updatePEtally($globs, $allRepeatTally, $PALhash, $accno, $seq1Info, $seq2Info, length($seq1), length($seq2));
	#$updateTally($allRepeatTally, $seq2Info);
	
	if ( !($seq1Info) and !($seq2Info) ) { $seqID = $accno; }
	else {
		$allRepeatTally->{readsWithMicrosat}++ unless ( not $seq1Info );
		$allRepeatTally->{readsWithMicrosat}++ unless ( not $seq2Info );
		$seqID = $seqInfo . "$accno";
		if($findPrimers) {
			$param_ref->{SEQUENCE_ID} = $seqID;
			updatePr3InfoForPairEnds($seq1Info,
						 $seq2Info,
						 $seq1, 
						 $seq2, 
						 $param_ref,
						 $allRepeatTally);
			foreach my $parameter (keys(%$param_ref)) {
				print $boulderIO_FH "$parameter=$$param_ref{$parameter}\n";
			}
			print $boulderIO_FH "=\n";
		}
	}
	return \$allRepeatTally;
}

sub scan454reads {
	my $readfile = shift;
	my $subref = shift;
	#my $globs = shift;
	my @args = @_;
	#my $readfile = $globs->{input454reads};

	my $readTitle = '';
	my $accno = '';
	my $seq = '';
	my $totalreads = 0;

	open(READS, "<$readfile") or die "could not open file $readfile\n";
	while(my $currLine = <READS>) {
		chomp($currLine);
		if ((substr($currLine, 0, 1)) eq '>') {
			if($totalreads > 0) {
				if($subref) { $subref->($accno, $seq, @args); }
			}
			$readTitle = substr($currLine, 1);
			if($readTitle =~ m/uaccno=/i) { $accno = (split(/uaccno=/, $readTitle))[1]; }
			elsif($readTitle =~ m/rank=/i) { $accno = (split(/\s/, $readTitle))[0]; }
			elsif($readTitle =~ m/length=\d+ xy=.+ region=\d{1,2} run=.+/i) {
				$accno = (split(/\s/, $readTitle))[0];
			}
			else { $accno = $readTitle; }
			$seq = '';
			$totalreads++;
		}
		else { $seq .= $currLine; }
	}
	# Act on last sequence
	if($totalreads > 0) {
		if($subref) { $subref->($accno, $seq, @args); }
	}
	close(READS);
	return $totalreads;
}

sub scanScarfPEreads {
	my $pair1file = shift;
	my $pair2file = shift;
	my $subRef = shift;
	my @args = @_;
	#my $globs = shift;
	#my $msatHash = shift;
	#my $tallyHash = shift;
	#my $pr3_FH = shift;
	
	my $totalReads = 0;

	open(PE1, "<$pair1file") or die "cannot open file $pair1file\n";
	open(PE2, "<$pair2file") or die "cannot open file $pair2file\n";
	
	my($line1, $line2, $title1, $title2, $seq1, $seq2, $qual1, $qual2);
	while( defined($line1 = <PE1>) and defined($line2 = <PE2>) ) {
		if ( ($line1 =~ m/^(.*:\d+:.+:[1-8]:\d+:\d+:\d+) [12]:[YN]:\d+:[ATCGN0]+:([ATCGN]+):(.*)/) ) { 
			$title1 = $1;
			$seq1 = $2;
			$qual1 = $3;
		}
		elsif ( ($line1 =~ m/(^.*-\d+_\d+(?::\d+){4}#\d\/)1:([ATCGN]+):(.+)/) ) {
			$title1 = $1;
			$seq1 = $2;
			$qual1 = $3;
		}
		elsif( ($line1 =~ m/^.*-(.+)#\d\/\d:([ATCGN]+):(.+)/) ) {
			$title1 = $1;
			$seq1 = $2;
			$qual1 = $3;
		}
		else {
			die "invalid sequence title in scarf file line:\n$line1\n";
		}
		if ( ($line2 =~ m/^(.*:\d+:.+:[1-8]:\d+:\d+:\d+) [12]:[YN]:\d+:[ATCGN0]+:([ATCGN]+):(.*)/) ) { 
			$title2 = $1;
			$seq2 = $2;
			$qual2 = $3;
		}
		elsif ( ($line2 =~ m/(^.*-\d+_\d+(?::\d+){4}#\d\/)2:([ATCGN]+):(.+)/) ) {
			$title2 = $1;
			$seq2 = $2;
			$qual2 = $3;
		}
		elsif( ($line2 =~ m/^.*-(.+)#\d\/\d:([ATCGN]+):(.+)/) ) {
			$title2 = $1;
			$seq2 = $2;
			$qual2 = $3;
		}
		else {
			die "invalid sequence title in scarf file line:\n$line2\n";
		}
		if (not validPEread ($title1, $title2, $seq1, $seq2, $qual1, $qual2) ) {
			print "Non-valid paired end read:\n$title1/1:$seq1:$qual1\n$title2/2:$seq2:$qual2\n";
			exit(1);
		}
		else { $totalReads++; }
		if ($subRef) { $subRef->($title1, $title2, $seq1, $seq2, $qual1, $qual2, @args); }
	}

	close(PE1);
	close(PE2);
	return $totalReads;
}

sub scanFastqPEreads {
	my $pair1file = shift;
	my $pair2file = shift;
	my $subRef = shift;
	my @args = @_;
	
	my $totalReads = 0;
	my $lineNum = 0;
	my $lineType = 0;

	open(PE1, "<$pair1file") or die "cannot open file $pair1file\n";
	open(PE2, "<$pair2file") or die "cannot open file $pair2file\n";
	
	my($line1, $line2, $title1, $title2, $seq1, $seq2, $qual1, $qual2);
	while( defined($line1 = <PE1>) and defined($line2 = <PE2>) ) {
		$lineType = $lineNum % 4;
		$line1 = stripStr($line1);
		$line2 = stripStr($line2);
		if ($line1 eq "") {
			next;	
		}
		elsif ($lineType == 0) {
			if($lineNum) {
				if (not validPEread ($title1, $title2, $seq1, $seq2, $qual1, $qual2) ) {
					print "Non-valid paired end read:\n$title1/1:$seq1:$qual1\n$title2/2:$seq2:$qual2\n";
					exit(1);
				}
				if ($subRef) { $subRef->($title1, $title2, $seq1, $seq2, $qual1, $qual2, @args); }
			}				
			if ( ($line1 =~ m/^@(.*:\d+:.+:[1-8]:\d+:\d+:\d+) [12]:[YN]:\d+:/) ) { $title1 = $1; }		
			elsif ( ($line1 =~ m/^@(.+-.+_.+(?::\d+){4}#[01ATCG]+)\/\d/) ) { $title1 = $1; }
			elsif( ($line1 =~ m/^@(.*)/ ) ) { $title1 = $1; }
			else { die "invalid sequence title in fastq file line:\n$line1\n"; }
			
			if ( ($line2 =~ m/^@(.*:\d+:.+:[1-8]:\d+:\d+:\d+) [12]:[YN]:\d+:/) ) { $title2 = $1; }		
			elsif ( ($line2 =~ m/^@(.+-.+_.+(?::\d+){4}#[01ATCG]+)\/\d/) ) { $title2 = $1; }
			elsif( ($line2 =~ m/^@(.*)/ ) ) { $title2 = $1; }
			else { die "invalid sequence title in fastq file line:\n$line2\n"; }
			$totalReads++;
		}
		elsif ($lineType == 1) {
			$seq1 = $line1;
			$seq2 = $line2;
		}
		elsif ($lineType == 2) {
			die "\nInvalid line: $line1 where +... expected\n" unless $line1 =~ m/^\+.*/;
			die "\nInvalid line: $line1 where +... expected\n" unless $line2 =~ m/^\+.*/;
		}
		elsif ($lineType == 3) {
			$qual1 = $line1;
			$qual2 = $line2;
		}
		$lineNum += 1;
	}
	
	if($totalReads) {
		if (not validPEread ($title1, $title2, $seq1, $seq2, $qual1, $qual2) ) {
			print "Non-valid paired end read:\n$title1/1:$seq1:$qual1\n$title2/2:$seq2:$qual2\n";
			exit(1);
		}
		if ($subRef) { $subRef->($title1, $title2, $seq1, $seq2, $qual1, $qual2, @args); }
	}	
	
	close(PE1);
	close(PE2);
	return $totalReads;
}


sub validPEread {
	my $title1 = shift;
	my $title2 = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $qual1 =shift;
	my $qual2 = shift;
	
	return 0 unless ($title1 eq $title2);
	return 0 unless ( length($seq1) == length($qual1) );
	return 0 unless ( length($seq2) == length($qual2) );
	
	return 1;
}

sub updatePr3InfoForPairEnds {
	my $seq1Info = shift;
	my $seq2Info = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $pr3params = shift;
	my $tallyHash = shift;
	my (@allRepeats1, @allRepeats2, $leftStart, $leftEnd, $total_microsats, $targetRegions, $shortest,
	   $rightStart, $rightEnd, $targetLength, $index, $read1cutoff, $read2cutoff, $tmpEnd, $tmpSize);
	my $read1len = length($seq1);
	my $read2len = length($seq2);
	
	my $leftDistToEnd = -1;
	my $rightDistToBeg = -1;
	$targetRegions = "";

	if ($seq1Info) {
		@allRepeats1 = split(/_/, $seq1Info);
		$index = (scalar(@allRepeats1)/3) - 1;
		$leftStart = $allRepeats1[$index * 3 + 1];
		$leftEnd = $allRepeats1[$index * 3 + 2];
		$leftDistToEnd = $read1len - $leftEnd;
	}
	if ($seq2Info) {
		@allRepeats2 = split(/_/, $seq2Info);
		$rightStart = $allRepeats2[1];
		$rightEnd = $allRepeats2[2];
		$rightDistToBeg = $rightStart - 1;		
	}
	if(($leftDistToEnd == -1) and ($rightDistToBeg == -1)) {
		die "feeding non-PAL to primer3: ($seq1Info)($seq2Info)";
	}
	elsif( ($leftDistToEnd >= 0) and ($rightDistToBeg == -1) ) {
		$targetLength = ($read1len - $leftStart) + 2;
		$pr3params->{SEQUENCE_TARGET} = "$leftStart,$targetLength";
	}
	elsif( ($leftDistToEnd == -1) and ($rightDistToBeg >= 0) ) {
		$targetLength = $rightEnd + 1;
		$rightStart = $read1len + 1;
		$pr3params->{SEQUENCE_TARGET} = "$rightStart,$targetLength";
	}
	elsif( $leftDistToEnd >= $rightDistToBeg ) {
		$targetLength = $rightEnd + 1;
		$rightStart = $read1len + 1;
		$pr3params->{SEQUENCE_TARGET} = "$rightStart,$targetLength";
	}
	else {
		$targetLength = ($read1len - $leftStart) + 2;
		$pr3params->{SEQUENCE_TARGET} = "$leftStart,$targetLength";
	}
	$pr3params->{SEQUENCE_TEMPLATE} = $seq1 . "n" . $seq2;
	my $maxSize = length($pr3params->{SEQUENCE_TEMPLATE});
	my $minSize = $pr3params->{PRIMER_MIN_SIZE} * 2 + $targetLength;
	my $mid = int(($maxSize + $minSize) / 2);
	$pr3params->{PRIMER_PRODUCT_SIZE_RANGE} = "$mid-$maxSize $minSize-$mid";
}

sub updatePr3Info {
	my $seqInfo = shift;
	my $seq = shift; 
	my $pr3params = shift;
	my $tallyHash = shift;

	my @allRepeats = split(/_/, $seqInfo);
	#shift(@allRepeats);
	my $start = $allRepeats[1];
	my $length = ($allRepeats[2] - $start);
	my $total_microsats = scalar(@allRepeats)/3;
	$pr3params->{SEQUENCE_TEMPLATE} = $seq;
	$pr3params->{SEQUENCE_INTERNAL_EXCLUDED_REGION} = "$start,$length";
	my $type = $allRepeats[0];
	my $tmpStart = 0; my $tmpLength = 0;
	for(my $index = 1; $index < $total_microsats; $index++) {
		$tmpStart = $allRepeats[$index * 3 + 1];
		$tmpLength = ($allRepeats[$index * 3 + 2] - $tmpStart);
		if($tmpLength > $length) {
			$type = $allRepeats[$index * 3];
			$start = $tmpStart;
			$length = $tmpLength;
		}
		$pr3params->{SEQUENCE_INTERNAL_EXCLUDED_REGION}.= " $tmpStart,$tmpLength";
	}
	$pr3params->{SEQUENCE_TARGET} = "$start,$length";
	$tallyHash->{$type}->{longest}++;
}

sub getTargetRegionFromID {
    my @inList = @_;
    my $seqID = $inList[0];
    my @allValues = split(/_/, $seqID);
    shift(@allValues);
    my $start = $allValues[1];
    my $length = ($allValues[2] - $start);
    my $total_microsats = (scalar(@allValues) / 3);
    my $excluded_regions = "$start,$length";
    my $type = $allValues[0];
    my $target_region = "";
    my $tmpStart = 0;
    my $tmpLength = 0;
    my $index = 1;
    while ($index < $total_microsats) {
	$tmpStart = $allValues[$index * 3 + 1];
	$tmpLength = ($allValues[$index * 3 + 2] - $tmpStart);
	if($tmpLength > $length) {
	    $type = $allValues[$index * 3];
	    $start = $tmpStart;
	    $length = $tmpLength;
	}
	$excluded_regions = $excluded_regions . " $tmpStart,$tmpLength";
	$index++;
    }
    $target_region = "$start,$length";
    return ($target_region, $excluded_regions, $type);
}

sub updateSumHash {
	my $sumHashRef = shift;
	my $accno = shift;
	my $monomerLen = shift;
	my $monomerType = shift;
	my $numRepeats = shift;
	my $fPrimer = shift;
	my $lPrimer = shift;
	my $repsAmplified = shift;

	my %infoHash  = ();
	$infoHash{monoLength} = $monomerLen;
	$infoHash{monoType} = $monomerType;
	$infoHash{numRepeats} = $numRepeats;
	$infoHash{fprimer} = $fPrimer;
	$infoHash{rprimer} = $lPrimer;
	$infoHash{repsAmplified} = $repsAmplified;

	$sumHashRef->{$accno} = \%infoHash;
}

sub updatePEsumHash {
	my $sumPEhashRef = shift;
	my $accno = shift;
	my $monomerLen = shift;
	my $monomerType = shift;
	my $numRepeats = shift;
	my $fPrimer = shift;
	my $lPrimer = shift;
	my $repsAmplified = shift;

	my %infoHash  = ();
	$infoHash{monoLength} = $monomerLen;
	$infoHash{monoType} = $monomerType;
	$infoHash{numRepeats} = $numRepeats;
	$infoHash{fprimer} = $fPrimer;
	$infoHash{rprimer} = $lPrimer;
	$infoHash{repsAmplified} = $repsAmplified;

	$sumPEhashRef->{$accno} = \%infoHash;
}

sub updatePrimerHash {
	my $hashRef = shift;
	my $primer = shift;
	
	if($hashRef->{$primer}) { $hashRef->{$primer}++; }
	else { $hashRef->{$primer} = 1; }
}

sub updatePrimerSizeRange {
	my $short = shift;
	my $long = shift;
	my $lprimer = shift;
	my $rprimer = shift;
	
	if(length($lprimer) < $short) { $short = length($lprimer); }
	if(length($rprimer) < $short) { $short = length($rprimer); }
	if(length($lprimer) > $long) { $long = length($lprimer); }
	if(length($rprimer) > $long) { $long = length($rprimer); }

	return ($short, $long);
}

sub printMicrosatInfo {
	my $globs = shift;
	my $microData = shift;
	my $outputFile = $globs->{MicrosatSumOut};
	my $findPrimers = $globs->{findPrimers};
	my $pe = $globs->{pairedEnd};
	
	my $satTypes = {};
	my $maxSize = $globs->{maxMonoSize};
	open(OUTPUT, ">$outputFile") or die "could not open file: $!\n";
	foreach my $var (keys(%$microData)) {
		if($var =~ m/\A[ATCGatcg]{2,$maxSize}\z/) {
			$satTypes->{$var} = $microData->{$var};
			#delete $microData->{$var};
		}
		else {
			my $value = $microData->{$var};
			if ($var eq "totalReads") {
				my $singles = $microData->{$var};
				if($pe) {
					$singles = $singles * 2;
					print OUTPUT "$var:\t$singles\t(2 x $value)\n";
				}
				else { print OUTPUT "$var:\t$singles\n"; }			
			}
			elsif ( ($var eq "readsWithPrimers") && (!($findPrimers)) ) { next; }
			else { print OUTPUT "$var:\t$value\n"; }
		}
	}
	if (!($pe)) {
		if($findPrimers) { print OUTPUT "\n\nMicrosat Type\tmonomer length\treads with loci\ttotal loci\tloci amplified\ttotal loci bases\tloci bases amplified\n"; }
		else { print OUTPUT "\n\nMicrosat Type\tmonomer length\ttotal loci\treads with loci\ttotal bases\n"; }
		foreach my $type (keys(%$satTypes)) {
			my $monoLen = length($type);
			my $all = $satTypes->{$type}->{totalLoci};
			my $longest = $satTypes->{$type}->{longest};
			my $primers = $satTypes->{$type}->{longestWithPrimers};
			my $bases = $satTypes->{$type}->{totalBases};
			my $ampedBases = $satTypes->{$type}->{basesAmplified};
			my $lociReads = $satTypes->{$type}->{readsWithLoci};
			my $ampedLoci = $satTypes->{$type}->{lociAmplified};
			if($findPrimers) { print OUTPUT "$type\t$monoLen\t$lociReads\t$all\t$ampedLoci\t$bases\t$ampedBases\n"; }
			else { print OUTPUT "$type\t$monoLen\t$all\t$lociReads\t$bases\n" }
		}
	}
	else {
		if($findPrimers) { print OUTPUT "\n\nMicrosat Type\tmonomer length\ttotal loci\tloci w/ primers\treads with loci\ttotal bases\textended\textended w/ primers\tspanning\tspanning w/ primers\n"; }
		else { print OUTPUT "\n\nMicrosat Type\tmonomer length\ttotal loci\treads with loci\ttotal bases\textended\tspanning\n"; }
		foreach my $type (keys(%$satTypes)) {
			my $monoLen = length($type);
			my $all = $satTypes->{$type}->{totalLoci};
			my $primers = $satTypes->{$type}->{lociWithPrimers};
			my $bases = $satTypes->{$type}->{totalBases};
			my $lociReads = $satTypes->{$type}->{readsWithLoci};
			my $expanded = $satTypes->{$type}->{potentialExtended};
			my $extPrimers = $satTypes->{$type}->{extendWithPrimers};
			my $span = $satTypes->{$type}->{potentialSpan};
			my $spanPrimer = $satTypes->{$type}->{spanWithPrimers};
			
			if($findPrimers) { print OUTPUT "$type\t$monoLen\t$all\t$primers\t$lociReads\t$bases\t$expanded\t$extPrimers\t$span\t$spanPrimer\n"; }
			else { print OUTPUT "$type\t$monoLen\t$all\t$lociReads\t$bases\t$expanded\t$span\n"; }
		}
	}
	close(OUTPUT);
}

sub writeResults {
	my $globs = shift;
	my $allMicrosatRef = shift;
	my $palHash = shift;
	my $dataOutputFile = $globs->{PALsummaryOut};
	my $pr3output = $globs->{primer3output};
	my $allReadsFile = $globs->{input454reads};
	my $findPrimers = $globs->{findPrimers};
	my $platform = $globs->{platform};
	my $pairEnd = $globs->{pairedEnd};

	if(!$findPrimers) {
		printMicrosatInfo($globs, $allMicrosatRef);
		printReadInfoNoPrimers($dataOutputFile, $palHash);
		return 0;
	}
	my %microsatSummary = ();
	my %prSingleHash = ();
	my %prPairHash = ();
	my $longestPrimerLen = 0;
	my $shortestPrimerLen = 99999999;
	my $total_multiple = 0;
	my $total_seqs = 0;
	my $total_nt = 0;
	my $curr_seq = "";
    
	#open(OUTPUT, ">$dataOutputFile");
	open(PR3IN, "<$pr3output") or die "could not open file $pr3output\n";
	my $type = "";
	my $mono_size = 0;
	my $primersReturned = 0;
	my $lprimerSeq = "";
	my $rprimerSeq = "";
	my $numRepeats = 0;
	my $totalCompound = 0;
	my $totalBroken = 0;
	my $seqID = "";
	my $accno = "";
	my $repsBetweenPrimers = "";
	while (my $curr_line = <PR3IN>){
		chomp($curr_line);
		if($curr_line eq '=') {
			$total_seqs++;
			if ($primersReturned > 0) {
				$repsBetweenPrimers = getNumberRepeatsBetweenPrimers($curr_seq, $seqID, $lprimerSeq, $rprimerSeq, $allMicrosatRef);
				($shortestPrimerLen, $longestPrimerLen) = updatePrimerSizeRange($shortestPrimerLen, $longestPrimerLen, $lprimerSeq, $rprimerSeq);				
				updatePrimerHash(\%prSingleHash, $lprimerSeq);
				updatePrimerHash(\%prSingleHash, $rprimerSeq);
				updatePrPairHash4pals(\%prPairHash, $lprimerSeq, $rprimerSeq);
				$allMicrosatRef->{readsWithPrimers} += 1;
			}
			#$allMicrosats{$type}[0]++;
			#$allMicrosats{$type}[1]+= $primersReturned;
			updateSumHash(\%microsatSummary, $accno, $mono_size, $type, $numRepeats, $lprimerSeq, $rprimerSeq, $repsBetweenPrimers);
			
			$lprimerSeq = "";
			$rprimerSeq = "";
			$primersReturned = 0;
			$repsBetweenPrimers = "";
		}
		my @values = split(/=/, $curr_line);
		if((scalar(@values)) == 2) {
			if($values[0] eq "SEQUENCE_TEMPLATE") {
				$curr_seq = $values[1];
			}
			elsif ($values[0] eq "PRIMER_LEFT_0_SEQUENCE") {
				$lprimerSeq = $values[1];
			}
			elsif ($values[0] eq "PRIMER_RIGHT_0_SEQUENCE") {
				$rprimerSeq = $values[1];
			}
			elsif ($values[0] eq "SEQUENCE_ID") {
				$seqID = $values[1];
				$accno = (split(/_/, $seqID))[0];
				my @seqTypeList = getLongestType($seqID);
				$type = $seqTypeList[0];
				$total_multiple += ($seqTypeList[1] - 1);
				$numRepeats = $seqTypeList[2];
				$mono_size = length($type);
				$totalCompound += $seqTypeList[3];
				$totalBroken += $seqTypeList[4];
			}
			elsif ($values[0] eq "PRIMER_PAIR_NUM_RETURNED") {
				$primersReturned = $values[1];
			}
		}
	}
	my $prNameHash = createPrimerNames(\%prSingleHash, $globs->{prNamePrefix});
	my $singlePrAllReadHash = createSinglePrAllReadHash(\%prSingleHash);
	#my $pairPrAllReadHash = createPairPrAllReadHash(\%prPairHash);
	scan454reads($allReadsFile, \&scanPrimersAlongSeq, $shortestPrimerLen, $longestPrimerLen, $singlePrAllReadHash, \%prPairHash);
	printMicrosatReadFile($dataOutputFile, \%microsatSummary, $prNameHash, $singlePrAllReadHash, \%prPairHash);
	printMicrosatInfo($globs, $allMicrosatRef);

	close(PR3IN);
	return 0;
}

sub writePEresults {
	my $globs = shift;
	my $allMicrosatRef = shift;
	my $palHash = shift;
	my $dataOutputFile = $globs->{PALsummaryOut};
	my $pr3output = $globs->{primer3output};
	#my $allReadsFile = $globs->{input454reads}; # Original fasta read file
	my $readFile1 = $globs->{inputReadFile};
	my $readFile2 = $globs->{pairedReadFile};
	my $findPrimers = $globs->{findPrimers};
	my $platform = $globs->{platform};
	my $pairEnd = $globs->{pairedEnd};

	if(!$findPrimers) {
		printMicrosatInfo($globs, $allMicrosatRef);
		printMicrosatPEreadFile($dataOutputFile, $palHash);
		return 0;
	}
	my %microsatSummary = ();
	my %allMicrosats = %$allMicrosatRef;
	my %prSingleHash = ();
	my %prPairHash = ();
	my $longestPrimerLen = 0;
	my $shortestPrimerLen = 99999999;
	my $total_multiple = 0;
	my $total_seqs = 0;
	my $total_nt = 0;
	my $curr_seq = "";
    
	#open(OUTPUT, ">$dataOutputFile");
	open(PR3IN, "<$pr3output") or die "could not open file $pr3output\n";
	my $type = "";
	my $mono_size = 0;
	my $primersReturned = 0;
	my $lprimerSeq = "";
	my $rprimerSeq = "";
	my $numRepeats = 0;
	my $totalCompound = 0;
	my $totalBroken = 0;
	my $seqID = "";
	my $accno = "";
	my $repsBetweenPrimers = "";
	my $repBasesBetweenPrimers = 0;
	while (my $curr_line = <PR3IN>){
		chomp($curr_line);
		if($curr_line eq '=') {
			$total_seqs++;
			if ($primersReturned > 0) {
				($repBasesBetweenPrimers, $repsBetweenPrimers) = getNumberRepeatsBetweenPrimers4PE($curr_seq,
																								   $seqID,
																								   $accno,
																								   $lprimerSeq,
																								   $rprimerSeq,
																								   $palHash,
																								   \%allMicrosats);
				($shortestPrimerLen, $longestPrimerLen) = updatePrimerSizeRange($shortestPrimerLen,
																				$longestPrimerLen,
																				$lprimerSeq,
																				$rprimerSeq);				
				updatePrimerHash(\%prSingleHash, $lprimerSeq);
				updatePrimerHash(\%prSingleHash, $rprimerSeq);
				updatePrPairHash4pals(\%prPairHash, $lprimerSeq, $rprimerSeq);
			}
			#$allMicrosats{$type}[0]++;
			#$allMicrosats{$type}[1]+= $primersReturned;
			updatePEsumHash(\%microsatSummary, $accno, $mono_size, $type, $numRepeats, $lprimerSeq, $rprimerSeq, $repsBetweenPrimers);
			
			$lprimerSeq = "";
			$rprimerSeq = "";
			$primersReturned = 0;
			$repsBetweenPrimers = "";
		}
		my @values = split(/=/, $curr_line);
		if((scalar(@values)) == 2) {
			if($values[0] eq "SEQUENCE_TEMPLATE") {
				$curr_seq = $values[1];
			}
			elsif ($values[0] eq "PRIMER_LEFT_0_SEQUENCE") {
				$lprimerSeq = $values[1];
			}
			elsif ($values[0] eq "PRIMER_RIGHT_0_SEQUENCE") {
				$rprimerSeq = $values[1];
			}
			elsif ($values[0] eq "SEQUENCE_ID") {
				$seqID = $values[1];
				#$seqID =~ /^(\((?:[ATCG]*_*\d*_*\d*_*\))\((?:[ATCG]*_*\d*_*\d*_*)\))/;
				if ( ! ($seqID =~ /^\((.*)\)\((.*)\)(.+)/ ) ) {
					die "$seqID is not in the correct format\n";
				}
				my $seq1Info = $1;
				my $seq2Info = $2;
				$accno = $3;
				#$seq1Info = substr($seq1Info, 1);
				#$seq2Info = substr($seq2Info, 0, length($seq2Info) - 1);
				my @seqInfo = ($seq1Info, $seq2Info);
				my $seqCombinedInfo = combineSeqInfo($seq1Info, $seq2Info);
				#$accno = $';
				#$accno = (split(/_/, $seqID))[0];
				foreach my $seqInfo (@seqInfo) {
					if ($seqInfo eq ""){
						next;
					}
					my @seqTypeList = getLongestType("X_" . $seqInfo);
					$type = $seqTypeList[0];
					$total_multiple += ($seqTypeList[1] - 1);
					$numRepeats = $seqTypeList[2];
					$mono_size = length($type);
					$totalCompound += $seqTypeList[3];
					$totalBroken += $seqTypeList[4];
				}
			}
			elsif ($values[0] eq "PRIMER_PAIR_NUM_RETURNED") {
				$primersReturned = $values[1];
			}
		}
	}
	my $prNameHash = createPrimerNames(\%prSingleHash, $globs->{prNamePrefix});
	my $singlePrAllReadHash = createSinglePrAllReadHash(\%prSingleHash);
	#my $pairPrAllReadHash = createPairPrAllReadHash(\%prPairHash);
	#scan454reads($allReadsFile, \&scanPrimersAlongSeq, $shortestPrimerLen, $longestPrimerLen, $singlePrAllReadHash, \%prPairHash);
	if ($globals->{inputFormat} eq "scarf") {
		scanScarfPEreads($globals->{inputReadFile},
				    $globals->{pairedReadFile},
				    \&scanPrimersAlongPEseqs,
				    $shortestPrimerLen,
				    $longestPrimerLen,
				    $singlePrAllReadHash,
				    \%prPairHash);
	}
	elsif ($globals->{inputFormat} eq "fastq") {
		scanFastqPEreads($globals->{inputReadFile},
				    $globals->{pairedReadFile},
				    \&scanPrimersAlongPEseqs,
				    $shortestPrimerLen,
				    $longestPrimerLen,
				    $singlePrAllReadHash,
				    \%prPairHash);
	}
	#printMicrosatReadFile($dataOutputFile, \%microsatSummary, $prNameHash, $singlePrAllReadHash, \%prPairHash);
	printMicrosatPEreadFile($dataOutputFile, $palHash, $prNameHash, $singlePrAllReadHash, \%prPairHash);
	printMicrosatInfo($globs, $allMicrosatRef);

	close(PR3IN);
	return 0;


}

sub combineSeqInfo {
	my $seq1Info = shift;
	my $seq2Info = shift;
	if( $seq1Info and $seq2Info ) { return $seq1Info . "_" . $seq2Info; }
	elsif( $seq1Info and ( !($seq2Info) ) ) { return $seq1Info; }
	elsif( (!($seq1Info)) and ($seq2Info) ) { return $seq2Info; }
	else { return ""; }
}

sub scanPrimersAlongSeq {
	my $accno = shift;
	my $seq = shift;
	my $shortest = shift;
	my $longest = shift;
	my $singleAllReadHash = shift;
	my $primerPairHash = shift;

	my @lengths2scan = ($shortest .. $longest);
	my @primersFound = ();
	my $seqLen = length($seq);

	for(my $offset=0; $offset < ($seqLen - $shortest) + 1; $offset++) {
		foreach my $len (@lengths2scan) {
			if($offset + $len < $seqLen + 1) {
				my $frag = substr($seq, $offset, $len);
				if(defined $singleAllReadHash->{$frag}) { 
					$singleAllReadHash->{$frag}++;
					my $revFrag = reverseComplement($frag);
					$singleAllReadHash->{$revFrag}++;
					push(@primersFound, $frag);
				}
			}
		}
	}
	if(scalar(@primersFound) >= 2) {
		updatePrPairHash4all($primerPairHash, @primersFound);
	}
}

sub scanPrimersAlongPEseqs {
	my $title1 = shift;
	my $title2 = shift;
	my $seq1 = shift;
	my $seq2 = shift; shift; shift;
	my $shortest = shift;
	my $longest = shift;
	my $singleAllReadHash = shift;
	my $primerPairHash = shift;

	my $seq = $seq1 . "n" . reverseComplement($seq2);
	my @lengths2scan = ($shortest .. $longest);
	my @primersFound = ();
	my $seqLen = length($seq);

	for(my $offset=0; $offset < ($seqLen - $shortest) + 1; $offset++) {
		foreach my $len (@lengths2scan) {
			if($offset + $len < $seqLen + 1) {
				my $frag = substr($seq, $offset, $len);
				if(defined $singleAllReadHash->{$frag}) { 
					$singleAllReadHash->{$frag}++;
					my $revFrag = reverseComplement($frag);
					$singleAllReadHash->{$revFrag}++;
					push(@primersFound, $frag);
				}
			}
		}
	}
	if(scalar(@primersFound) >= 2) {
		updatePrPairHash4all($primerPairHash, @primersFound);
	}
}

sub updatePrPairHash4all {
	my $pairHash = shift;
	my @primerList = @_;
	
	for(my $i=0; $i < scalar(@primerList) - 1; $i++) {
		my $currPrimer = $primerList[$i];
		my $pair = hasPair($pairHash, $currPrimer, $i, @primerList);
		if($pair) {
			my $rcPrimer = reverseComplement($currPrimer);
			my $rcPair = reverseComplement($pair);
			$pairHash->{$currPrimer}->{$pair}->{all}++;
			$pairHash->{$rcPair}->{$rcPrimer}->{all}++;
		}
	}
}

sub hasPair {
	my $pairHash = shift;
	my $primer = shift;
	my $index = shift;
	my @prList = @_;
	my $pair = "";

	if($pairHash->{$primer}->{$prList[$index]}) {
		return $prList[$index];
	}
	elsif($index == scalar(@prList) - 1) {
		return "";
	}
	else {
		$pair = hasPair($pairHash, $primer, $index + 1, @prList);
		return $pair;
	}
}

sub printReadInfoNoPrimers {
	my $dataOutputFile = shift;
	my $palHash = shift;
	
	open(OUT, ">$dataOutputFile") or die "could not find $dataOutputFile\n";
	
	print OUT "SequenceID\tLargest motif size\tLargest Motif\tMotifs(bases)\tTotal motif bases\n";
	
	foreach my $accno (keys(%$palHash)) {
		my $long = $palHash->{$accno}->{longest};
		my $size = length($long);
		my $loci = $palHash->{$accno}->{allLoci};
		my $bases = $palHash->{$accno}->{allLociBases};
		print OUT "$accno\t$size\t$long\t$loci\t$bases\n";
	}
	close(OUT);
}

sub printMicrosatReadFile {
	my $outFile = shift;
	my $usatReadSummary = shift;
	my $nameHash = shift;
	my $primerAllReads = shift;
	my $primerPairHash = shift;

	open(OUT, ">$outFile") or die "could not open file $outFile\n";
	print OUT "SequenceID\tRepeat Motif Size\tRepeat Motif\tNumber Tandem Repeats\tPrimer Designed (1=y,0=n)\tF Primer Name\t";
	print OUT "Forward Primer\tR Primer Name\tReverse Primer\tTotal Repeats In Amplicon\tOccurances of Forward Primer in Reads\t";
	print OUT "Occurances of Reverse Primer in Reads\tOccurances of Amplifiable Primer Pair in Reads\t";
	print OUT "Occurances of Amplifiable Primer Pair in PALs\n";
	
	foreach my $readName (keys(%$usatReadSummary)) {
		my $motifSize = $usatReadSummary->{$readName}->{monoLength};
		my $motifType = $usatReadSummary->{$readName}->{monoType};
		my $motifReps = $usatReadSummary->{$readName}->{numRepeats};
		my $fPrimer = $usatReadSummary->{$readName}->{fprimer};
		my $rPrimer = $usatReadSummary->{$readName}->{rprimer};
		my $repsAmped = $usatReadSummary->{$readName}->{repsAmplified};
		my $ifRet = 0;
		my $fInReads = '';
		my $rInReads = '';
		my $fName = '';
		my $rName = '';
		my $PALpairs = '';
		my $allReadPairs = '';
		if($fPrimer) { 
			$ifRet = 1;
			$fInReads = $primerAllReads->{$fPrimer};
			$rInReads = $primerAllReads->{$rPrimer};
			$fName = $nameHash->{$fPrimer};
			$rName = $nameHash->{$rPrimer};
			my $rcRprimer = reverseComplement($rPrimer);
			$PALpairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{pal};
			$allReadPairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{all};
		}
		print OUT "$readName\t$motifSize\t$motifType\t$motifReps\t$ifRet\t$fName\t$fPrimer\t$rName\t$rPrimer\t$repsAmped\t$fInReads\t$rInReads\t$allReadPairs\t$PALpairs\n";
	}

	close(OUT);
}

sub printMicrosatPEreadFile {
	my $findPrimers = $globals->{findPrimers};
	my $outFile = shift;
	my $palHash = shift;
	my $nameHash = shift unless !($findPrimers); 
	my $primerAllReads = shift unless !($findPrimers);
	my $primerPairHash = shift unless !($findPrimers);
	

	my($fPrimer, $rPrimer, $fInReads, $rInReads, $fName, $rName, $ifRet, $PALpairs, $allReadPairs, $lociInAplicon,
	   $lociBasesInAmp, $lociBasesInAmpm, $primersSepReads, $extPrimer, $spanPrimer);
	
	open(OUT, ">$outFile") or die "could not open file $outFile\n";
	print OUT "readPairID\tMotifs(bases)\tBases in all Motifs\tPossible Extended\tPossible Spanning";
	if ($findPrimers) {
		print OUT "\tPrimers found (1=y,0=n)\tF Primer Name\tForward Primer\tR Primer Name\tReverse Primer\t";
		print OUT "Amplicon Motifs\tNumber motif bases in amplicon\tPrimers on sep reads\tExtend with primers\t";
		print OUT "Spand with primers\tOccurances of Forward Primer in Reads\t";
		print OUT "Occurances of Reverse Primer in Reads\tOccurances of Amplifiable Primer Pair in Reads\t";
		print OUT "Occurances of Amplifiable Primer Pair in PALs\n";
	}
	else { print OUT "\n"; }
	
	
	foreach my $readName (keys(%$palHash)) {
		my $motifs = $palHash->{$readName}->{allLoci};
		my $allBases = $palHash->{$readName}->{allLociBases};
		my $extended = $palHash->{$readName}->{extend}; $extended = "" unless (defined($extended));
		my $span = $palHash->{$readName}->{span}; $span = "" unless (defined($span));
		my $primersFound = 0;
		if ($findPrimers) {
			$fPrimer = $palHash->{$readName}->{lprimer};
			$rPrimer = $palHash->{$readName}->{rprimer};
			if($fPrimer) { $primersFound = 1; }
		}
		if ($primersFound) {
			$lociInAplicon = $palHash->{$readName}->{lociWithPrimers};
			$lociBasesInAmp = $palHash->{$readName}->{repeatBasesInAmplicon};
			$primersSepReads = $palHash->{$readName}->{primersSepReads};
			$extPrimer = $palHash->{$readName}->{extPrimer}; $extPrimer = "" unless defined($extPrimer);
			$spanPrimer = $palHash->{$readName}->{spanWithPrimers}; $spanPrimer = "" unless defined($spanPrimer);
			$fInReads = $primerAllReads->{$fPrimer};
			$rInReads = $primerAllReads->{$rPrimer};
			$fName = $nameHash->{$fPrimer};
			$rName = $nameHash->{$rPrimer};
			my $rcRprimer = reverseComplement($rPrimer);
			$PALpairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{pal};
			$allReadPairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{all};
		}
		else {
			$lociInAplicon = "";
			$lociBasesInAmp = "";
			$primersSepReads = "";
			$extPrimer = "";
			$spanPrimer = "";
			$fInReads = "";
			$rInReads = "";
			$fName = "";
			$fPrimer = "";
			$rPrimer = "";
			$rName = "";
			$PALpairs = "";
			$allReadPairs = "";
		}
		#my $ifRet = 0;
		#my $fInReads = '';
		#my $rInReads = '';
		#my $fName = '';
		#my $rName = '';
		#my $PALpairs = '';
		#my $allReadPairs = '';
		#if($fPrimer) { 
		#	$ifRet = 1;
		#	$fInReads = $primerAllReads->{$fPrimer};
		#	$rInReads = $primerAllReads->{$rPrimer};
		#	$fName = $nameHash->{$fPrimer};
		#	$rName = $nameHash->{$rPrimer};
		#	my $rcRprimer = reverseComplement($rPrimer);
		#	$PALpairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{pal};
		#	$allReadPairs = $primerPairHash->{$fPrimer}->{$rcRprimer}->{all};
		#}
		print OUT "$readName\t$motifs\t$allBases\t$extended\t$span";
		if($findPrimers){
		#print OUT "$readName\t$motifSize\t$motifType\t$motifReps\t$ifRet\t$fName\t$fPrimer\t$rName\t$rPrimer\t$repsAmped\t$fInReads\t$rInReads\t$allReadPairs\t$PALpairs\n";
			print OUT "\t$primersFound\t$fName\t$fPrimer\t$rName\t$rPrimer\t$lociInAplicon\t$lociBasesInAmp\t$primersSepReads\t";
			print OUT "$extPrimer\t$spanPrimer\t$fInReads\t$rInReads\t$allReadPairs\t$PALpairs\n";
		} else { print OUT "\n"; }
	}

	close(OUT);
}

sub createSinglePrAllReadHash {
	my $satPrimerHash = shift;
	my %newHash = ();

	foreach my $primer (keys(%$satPrimerHash)) {
		$newHash{$primer} = 0;
		my $rcPrimer = reverseComplement($primer);
		$newHash{$rcPrimer} = 0;
	}
	return \%newHash;
}

sub createPairPrAllReadHash {
	my $satPrimerPairHash = shift;
	my %newHash = ();

	foreach my $primerPair (keys(%$satPrimerPairHash)) {
		my($fPrimer, $rPrimer) = split(/-/, $primerPair);
		my $rcFprimer = reverseComplement($fPrimer);
		my $rcRprimer = reverseComplement($rPrimer);
		$newHash{$fPrimer}->{$rcRprimer}->{all} = 0;
		$newHash{$rPrimer}->{$rcFprimer}->{all} = 0;
	}
	return \%newHash;
}


sub updatePrPairHash4pals {
	my $pairHash = shift;
	my $fPrimer = shift;
	my $rPrimer = shift;

	my $rcFprimer = reverseComplement($fPrimer);
	my $rcRprimer = reverseComplement($rPrimer);
	
	if(defined($pairHash->{$fPrimer}->{$rcRprimer}->{pal})) {
		$pairHash->{$fPrimer}->{$rcRprimer}->{pal}++;
		$pairHash->{$rPrimer}->{$rcFprimer}->{pal}++;
	}
	else {
		$pairHash->{$fPrimer}->{$rcRprimer}->{pal} = 1;
		$pairHash->{$fPrimer}->{$rcRprimer}->{all} = 0;
		$pairHash->{$rPrimer}->{$rcFprimer}->{pal} = 1;
		$pairHash->{$rPrimer}->{$rcFprimer}->{all} = 0;
	}
	return 0;
}

sub createPrimerNames {
	my $primerHash = shift;
	my $prefix = shift;
	
	my %nameHash = ();
	my @primerList = keys(%$primerHash);
	my $totalPrimers = scalar(@primerList);
	if($totalPrimers > 0) {
		my $current = 1;
		my $digits = int(log10($totalPrimers)) + 1;
		foreach my $primer (@primerList) {
			my $zeros = $digits - (int(log10($current)) + 1);
			my $zeroStr = "";
			for (my $i=0; $i<$zeros; $i++) { $zeroStr .= "0"; }
			$nameHash{$primer} = "$prefix" . "$zeroStr" . "$current";
			$current++;
		}
	}
	return \%nameHash;
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub stripStr {
    my @input = @_;
    my $string = $input[0];
    $string =~ s/\A[\s\t\n\r\f]*//;
    $string =~ s/[\s\t\n\r\f]*$//;
    return $string;
}

sub getLongestType {
    my $seqID = shift;
    my @allValues = split(/_/, $seqID);
    shift(@allValues);
    if ( scalar(@allValues) % 3 != 0 )  {
    	die "getLongestType: invalid seqID '$seqID'\n";
    }
    my $numValues = scalar(@allValues)/3;
    my $type = $allValues[0];
    my $firstType = $allValues[0];
    my $mono_size = length($type);
    my $repeats = 0;
    my $start = $allValues[1];
    my $stop = $allValues[2];
    my $length = ($stop - $start) + 1;
    my $tmpLength = 0;
    my $tmpStart = 0;
    my $tmpStop = 0;
    my $index = 1;
    my $compound = 0;
    my $broken = 0;
    my $curr_type = "";
    while ($index < $numValues) {
	$tmpStart = $allValues[($index * 3) + 1];
	$tmpStop = $allValues[($index * 3) + 2];
	$tmpLength = ($tmpStop - $tmpStart) + 1;
	$curr_type = $allValues[($index * 3)];
	if ($curr_type ne $firstType) {
	    $compound = 1;
	}
	if ($length < $tmpLength) {
	    $length = $tmpLength;
	    $type = $allValues[$index * 3];
	    $mono_size = length($type);
	}
	$index++;
    }
    if ( ($numValues > 1) && ($compound == 0)) {
	$broken = 1;
    }
    $repeats = int($length/$mono_size);
    return ($type, $numValues, $repeats, $compound, $broken);
}

sub getNumberRepeatsBetweenPrimers {
    my $seq = shift;
    my $seqID = shift;
    my $lPrimer = shift;
    my $rPrimer = shift;
    my $allMicrosatRef = shift;

    my $leftBookEnd = index($seq, $lPrimer);
    $rPrimer = (reverseComplement($rPrimer))[0];
    my $rightBookEnd = rindex($seq, $rPrimer);

    my @allValues = split(/_/, $seqID);
    shift(@allValues);
    my $totalMicrosats = scalar(@allValues)/3;

    my $currSat = $allValues[0];
    my $currStart = $allValues[1];
    my $currStop = $allValues[2];

    my $index = 0;
    my $totalRepeats = 0;

    while ($index < $totalMicrosats) {
	$currSat = $allValues[($index * 3)];
	$currStart = $allValues[($index * 3) + 1];
	$currStop = $allValues[($index * 3) + 2];
	if (($currStart >= ($leftBookEnd - 1)) && ($currStop <= ($rightBookEnd + 1))) {
	    my $currLength = ($currStop - $currStart) + 1;
	    my $currMonoSize = length($currSat);
	    my $currRepeats = int($currLength/$currMonoSize);
	    $totalRepeats += $currRepeats;
	    $allMicrosatRef->{$currSat}->{basesAmplified} += $currLength;
	    $allMicrosatRef->{$currSat}->{lociAmplified} += 1;
	}
	$index++;
    }
    return $totalRepeats;
}

sub getNumberRepeatsBetweenPrimers4PE {
    my $seq = shift;
    my $seqInfo = shift;
    my $accno = shift;
    my $lPrimer = shift;
    my $rPrimer = shift;
    my $palHash = shift;
    my $allMicrosatHash = shift;
    
    my ($currSeqInfo, $seq1Info, $seq2Info, $currSeq, $leftBookEnd, $rightBookEnd, $currLen,
    	$seq1info, $seq2info, $lExtType, $lExtStart, $lExtEnd, $rExtType, $rExtStart, $rExtEnd, $spacer,
    	$currStart, $currEnd, $currType, $lseq, $rseq, @seqs, $lastBaseInL, $lastBaseInR, $linl, $rinr,
    	$pairLbookEnd, $pairRbookEnd);
    

    @seqs = split(/n+/, $seq);
    $lseq = $seqs[0];
    $rseq = $seqs[1];
    
    $lastBaseInL = length($lseq);
    $lastBaseInR = length($rseq);
    
    $leftBookEnd = index($lseq, $lPrimer) + 1;
    if ($leftBookEnd == 0) { 
    	$leftBookEnd = index($rseq, $lPrimer) + 1;
    	$linl = 0;
    }
    else { $linl = 1; }
    
    
    $palHash->{$accno}->{lprimer} = $lPrimer;
	$palHash->{$accno}->{rprimer} = $rPrimer;
    $rPrimer = reverseComplement($rPrimer);
    $rightBookEnd = rindex($rseq, $rPrimer) + 1;
    if ($rightBookEnd == 0) {
    	$rightBookEnd = index($lseq, $rPrimer) + 1;
    	$rinr = 0;
    }
    else { $rinr = 1; }
    
    if ($linl == -1 or $rinr == -1) {
    	die "could not find both primers $lPrimer and $rPrimer in seq $seq for read pair $seqInfo\n";
    }
	$spacer = $lastBaseInL + (length($seq) - ($lastBaseInL + $lastBaseInR));
	if($linl) { $pairLbookEnd = $leftBookEnd; }
	else { $pairLbookEnd = $leftBookEnd + $spacer; }
	if($rinr) { $pairRbookEnd = $rightBookEnd + $spacer; }
	else { $pairRbookEnd = $rightBookEnd; }
	
	if($linl and $rinr) { $palHash->{$accno}->{primersSepReads} = 1; }
	else { $palHash->{$accno}->{primersSepReads} = 0; }
	
    if ( $seqInfo =~ /^\((.*)\)\((.*)\).*/ ) {
    	$seq1info = $1;
    	$seq2info = $2;
    }
    else {
    	die "Sequence Tag $seqInfo not in correct format\n";
    }
    my @lfields = split(/_/, $seq1info);
    my @rfields = split(/_/, $seq2info);
    my $lsats = int(scalar(@lfields)/3);
    my $rsats = int(scalar(@rfields)/3);
	
	my $totalRepsBetweenPrimers = 0;
	my $totalRepBasesBetweenPrimers = 0;
	
	$lExtType = "";
	for (my $i= 0; $i < $lsats; $i++ ) {
		$currType = $lfields[$i * 3];
		$currStart  = $lfields[($i * 3) + 1];
		#if (!($linl)) { $currStart = $currStart + $spacer; }
		$currEnd = $lfields[($i * 3) + 2];
		#if (!($linl)) { $currEnd = $currEnd + $spacer; }
		$currLen = ($currEnd - $currStart) + 1;
		if ( ($currEnd < $pairRbookEnd) and ($currStart > $pairLbookEnd) ) {
			if ($currEnd == $lastBaseInL) {
				$lExtType = $currType;
				$microsatTally->{$currType}->{extendWithPrimers}++;
				$palHash->{$accno}->{extPrimer} .= "$currType ";
			}
			$microsatTally->{$currType}->{lociWithPrimers}++;
			$totalRepBasesBetweenPrimers += ($currEnd - $currStart) + 1;
			$totalRepsBetweenPrimers += int((($currEnd - $currStart) + 1)/length($currType));
			$palHash->{$accno}->{lociWithPrimers} .= "$currType($currLen) ";
		}
	}
	
	$rExtType = "";
	for (my $j= 0; $j < $rsats; $j++ ) {
		$currType = $rfields[$j * 3];
		$currStart  = $rfields[($j * 3) + 1];
		$currStart += $spacer;
		$currEnd = $rfields[($j * 3) + 2];
		$currEnd += $spacer;
		$currLen = ($currEnd - $currStart) + 1;
		if ( ($currStart < $pairRbookEnd) and ($currEnd > $pairLbookEnd) ) {
			if ($currStart == $spacer + 1) {
				$rExtType = $currType;
				$microsatTally->{$currType}->{extendWithPrimers}++;
				$palHash->{$accno}->{extPrimer} .= "$currType ";
			}
			$microsatTally->{$currType}->{lociWithPrimers}++;
			$totalRepBasesBetweenPrimers += ($currEnd - $currStart) + 1;
			$totalRepsBetweenPrimers += int((($currEnd - $currStart) + 1)/length($currType));
			$palHash->{$accno}->{lociWithPrimers} .= "$currType($currLen) ";
		}
	
	}
    if ($rExtType and ($rExtType eq $lExtType)) {
    	if( ($leftBookEnd < index($seq, "n")) and ( rindex($seq, $rPrimer) > rindex($seq, "n") ) ) {
    		$microsatTally->{$rExtType}->{spanWithPrimers}++;
    		$palHash->{$accno}->{spanWithPrimers} = $rExtType;
    	}
    }
    $palHash->{$accno}->{repeatBasesInAmplicon} = $totalRepBasesBetweenPrimers;
    return ($totalRepBasesBetweenPrimers, $totalRepsBetweenPrimers);
}    


sub maskRegions {
    my ($seq, $excluded) = ((@_)[0], (@_)[1]);
    my $maskedSeq = "";
    my @excludedList = split(/ /, $excluded);
    foreach my $tuple (@excludedList) {
	my $start = (split(/,/, $tuple))[0];
	my $length = (split(/,/, $tuple))[1];
	my $subst = "";
	for(1..$length) { $subst = $subst . "N"; }
	#my $first3rd = substr($seq, 0, $start - 1);
	#my $second3rd = substr($seq, $start - 1, $length, $subst);
	substr($seq, $start - 1, $length, $subst);
	#my $third3rd = substr($seq, $start + $length);
	#$seq = $first3rd . $second3rd . $third3rd;
    }
    return $seq;
}

sub createAllNmers {
    my $nmerLength = $_[0];
    my @bases = qw(A T C G);
    my @nmerList = ();
    my @posIndex = (0);
    for (my $i = 1; $i < $nmerLength; $i++) {
	push(@posIndex, 0);
    }
    while($posIndex[0] != -1 ) {
	my $newNmer = "";
	for (my $j = 0; $j < $nmerLength; $j++) {
	    $newNmer = $newNmer . $bases[$posIndex[$j]];
	}
	push(@nmerList, $newNmer);
	@posIndex = updatePosIndex(@posIndex);
    }
    return @nmerList;
}

sub updatePosIndex {
    my @posIndex = @_;
    my $nmerLength = scalar(@posIndex);
    for (my $k = ($nmerLength - 1); $k >= 0; $k--) {
	if ($posIndex[$k] < 3) {
	    $posIndex[$k]++;
	    last;
	}
	elsif (($posIndex[$k] == 3) && ($k != 0)) {
	    $posIndex[$k] = 0;
	}
	elsif (($posIndex[$k] == 3) && ($k == 0)) {
	    return -1;
	}
    }
    return @posIndex;
}

sub makeGroup {
  my @inList = @_;
  my $curr_seq = $inList[0];
  my $seqSize = length($curr_seq);
  my $new_seq = "";
  my @group = ($curr_seq);
  for(my $i=0;$i < ($seqSize - 1); $i++){
    $new_seq = substr($curr_seq, ($seqSize - 1)) . substr($curr_seq, 0, ($seqSize - 1));
    push(@group, $new_seq);
    $curr_seq = $new_seq;
  }
  return @group;
}

sub removeFromArray {
  my @inList = @_;
  my $array_ref = $inList[0];
  my @array = @$array_ref;
  my $element = $inList[1];
  for (my $i=0;$i < scalar(@array); $i++) {
    if ($element eq $array[$i]) {
      splice(@array, $i, 1);
      $i--;
    }
  }
  return @array
}
sub groupAllNmers {
  my @list4mers = @_;
  my %groupHash = ();
  while (scalar(@list4mers) > 0) {
    my $leftInArray = scalar(@list4mers);
    my $curr4mer = shift(@list4mers);
    my @currGroup = makeGroup($curr4mer);
    $groupHash{$curr4mer} = \@currGroup;
    foreach my $item (@currGroup) {
      @list4mers = removeFromArray(\@list4mers, $item);
    }
  }
  return \%groupHash;
}

sub reverseComplement {
    my @seqList = @_;
    my $seqin = $seqList[0];
    my $rev_seq = reverse($seqin);
    $rev_seq =~ tr/ATCGatcg/TAGCtagc/;
    return $rev_seq;
}

sub removeRevCompGroups {
  my @inList = @_;
  my $hash_ref = $inList[0];
  my %big_hash = %$hash_ref;
  my $delete_self =0;
  
  my @allKeys = keys(%big_hash);
  my $index = 0;
  while(scalar(@allKeys) > 0) {
    $delete_self = 0;
    my $key = shift(@allKeys);
    my $seqLen = length($key);
    my $revKey = reverseComplement($key);
    foreach my $dictKey (keys(%big_hash)) {
      my $currListRef = $big_hash{$dictKey};
      my @currList = @$currListRef;
      for(my $i=0;$i<$seqLen;$i++) {
        if($revKey eq $currList[$i]) {
          $delete_self = 0;
          for(my $j=0;$j<$seqLen;$j++){
            if($key eq $currList[$j]){
              $delete_self = 1;
              last;
            }
          }
          if ($delete_self == 0) {
            $index++;
            delete $big_hash{$dictKey};
            @allKeys = removeFromArray(\@allKeys, $dictKey);
            last;
          }
        }
      }
    }
  }
  return \%big_hash;
}


sub createAllUniqueRepeats {
    my $nmerLength = $_[0];
    my @others2remove = ();
    foreach my $i (1..($nmerLength-1)) {
	if (($nmerLength % $i) == 0) { 
	    push(@others2remove, $i); 
	}
    }
    my @nmers2remove = smallerGroups2remove(\@others2remove, $nmerLength);
    my @allNmers = createAllNmers($nmerLength);
    foreach my $currNmer (@nmers2remove) {
	@allNmers = removeFromArray(\@allNmers, $currNmer);
    }
    my $allGroups = groupAllNmers(@allNmers);
    my $finalGroupRef = removeRevCompGroups($allGroups);
    return $finalGroupRef;
}

sub smallerGroups2remove {
    my @inList = @_;
    my $toRemoveRef = $inList[0];
    my @nmers2remove = @$toRemoveRef;
    my $nmerLength = $inList[1];
    my $multiplier = 0;
    my @all2remove = ();
    foreach my $len (@nmers2remove) {
	$multiplier = $nmerLength / $len;
	my @allmers = createAllNmers($len);
	foreach my $monomer (@allmers) {
	    my $currNmer = "";
	    for (1..$multiplier) { $currNmer = $currNmer . $monomer; }
	    push(@all2remove, $currNmer);
	}
    }
    return @all2remove;
}

sub createCountHash {
	my @microLengthList = @_;
	my $countHash = {};
	my $tmpHashRef;
	for my $microLength (@microLengthList) {
		$tmpHashRef = createAllUniqueRepeats($microLength);
		for my $satType (keys(%$tmpHashRef)){
			my $satGroupRef = $tmpHashRef -> {$satType};
			my @satGroup = @$satGroupRef;
			my $numSatsInGroup = scalar(@satGroup);
			for(my $i=0; $i < $numSatsInGroup; $i++) {
				my $currSat = $satGroup[$i];
				my $revSat = reverseComplement($currSat);
				
				$countHash->{$currSat} = {
							  reps => 0,
							  startPos => 0, 
							  satType => $satType
							 };
				$countHash->{$revSat} =  {
							  reps => 0,
							  startPos => 0, 
							  satType => $satType
							 };
			}
		}
	}
	return $countHash;
}

sub scanSeqAgainstHash {
    my $seq = shift;
    my $lookupHashRef = shift;
    my $globs = shift;
    my @chunkSizes = @_;
    my $seqID = "";
    foreach my $chunk (@chunkSizes){
	my $partialID = scanSeqWithWindow($seq, $chunk, $globs, $lookupHashRef);
	if ($partialID) {
	    if ($seqID eq "") {
		$seqID = $partialID;
	    }
	    else {
		$seqID = $seqID . "_$partialID";
	    }
	}
    }
    return $seqID;
}

sub createMicrosatHashCount {
	my $globals = shift;
	my @monoLengths = (2..$globals->{maxMonoSize});
	my $hashOfMicrosatHashes = {};
	for my $mono (@monoLengths) {
		for my $offset (0 .. $mono - 1) {
			$hashOfMicrosatHashes->{$mono}->{$offset} = {
				       prevSeq => "",
				       currSeq => "",
				       countHash => createCountHash($mono)
				      };
		}
	}
	return $hashOfMicrosatHashes;
}


sub scanSeqWithWindow {
	my $seq = shift;
	my $windowSize = shift;
	my $globs = shift;
	my $masterHashRef = shift;
	my $offsetHashRef = $masterHashRef->{$windowSize};
	my $seqID = "";
	my $seqLength = length($seq);
	my $thresholdLength = 12;
	my $minRepeats = $globs->{"$windowSize"."merMinReps"};
	#my $minRepeats = $thresholdLength/$windowSize;
	
	if ($seqLength < ($windowSize * 2)) {
		return "";
	}
	else{
		my @offsetList = (0 .. ($windowSize - 1));
		my @runningRepeats = ();
		for (1 .. $windowSize ) { push(@runningRepeats, 0); }
		#my %offsetHash = %$offsetHashRef;
		#my $currSeq = "";
		#my $prevSeq = "";
		#my $except = "";
		my $cutOut = 0;
		my $offset = 0;
		my $pos = 0;
		while ($pos + $offset + $windowSize < $seqLength + 1) { 
			foreach $offset (@offsetList) {
				my $currOffset = $offsetHashRef->{$offset};
				if ($pos + $offset + $windowSize < $seqLength + 1) {
					
					$currOffset->{prevSeq} = $currOffset->{currSeq};
					$currOffset->{currSeq} = substr($seq, ($pos + $offset), $windowSize);
					if (defined $currOffset->{countHash}->{$currOffset->{currSeq}}->{satType}) {
						if ($currOffset->{countHash}->{$currOffset->{currSeq}}->{reps} > 0) {
							$currOffset->{countHash}->{$currOffset->{currSeq}}->{reps}++;
							$runningRepeats[$offset] = $currOffset->{countHash}->{$currOffset->{currSeq}}->{reps};
						}
						elsif ($currOffset->{countHash}->{$currOffset->{currSeq}}->{reps} == 0) {
							if (defined $currOffset->{countHash}->{$currOffset->{prevSeq}}->{satType}) {
								if ($currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps} >= $minRepeats) {
									$runningRepeats[$offset] = $currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps};
									$cutOut = 1;
								}
								$currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps} = 0;
							}
							$currOffset->{countHash}->{$currOffset->{currSeq}}->{startPos} = $pos + $offset + 1;
							$currOffset->{countHash}->{$currOffset->{currSeq}}->{reps} = 1;
						}
					}
					else {
						if (defined $currOffset->{countHash}->{$currOffset->{prevSeq}}->{satType}) {
							if ($currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps} >= $minRepeats) {
								$runningRepeats[$offset] = $currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps};
								$cutOut = 1;
							}
							$currOffset->{countHash}->{$currOffset->{prevSeq}}->{reps} = 0;
						}
					}
				}
			}
			if ($cutOut == 1) {
				$seqID = constructSeqID($offsetHashRef, $seqID, $windowSize, $minRepeats, @runningRepeats);
				$cutOut = 0;
				for my $q (0 .. $windowSize - 1) {
					$runningRepeats[$q] = 0;
					if(defined $offsetHashRef->{$q}->{countHash}->{$offsetHashRef->{$q}->{prevSeq}}) {
						$offsetHashRef->{$q}->{countHash}->{$offsetHashRef->{$q}->{prevSeq}}->{reps} = 0;
					}
					if(defined $offsetHashRef->{$q}->{countHash}->{$offsetHashRef->{$q}->{currSeq}}) {
						$offsetHashRef->{$q}->{countHash}->{$offsetHashRef->{$q}->{currSeq}}->{reps} = 0;
					}
				}
			}
			$offset = 0;
			$pos+=$windowSize;
		}
		$seqID = constructSeqID($offsetHashRef, $seqID, $windowSize, $minRepeats, @runningRepeats);
		for my $k (0 .. $windowSize - 1) {
			if(defined $offsetHashRef->{$k}->{countHash}->{$offsetHashRef->{$k}->{prevSeq}}) {
				$offsetHashRef->{$k}->{countHash}->{$offsetHashRef->{$k}->{prevSeq}}->{reps} = 0;
			}
			if(defined $offsetHashRef->{$k}->{countHash}->{$offsetHashRef->{$k}->{currSeq}} ) {
				$offsetHashRef->{$k}->{countHash}->{$offsetHashRef->{$k}->{currSeq}}->{reps} = 0;
			}
		}
	}
	#clearCountHash($offsetHashRef);
	return $seqID;
}

sub constructSeqID {
	my $offsetHash = shift;
	my $seqID = shift;
	my $windowSize = shift;
	my $minRepeats = shift;
	my @runningRepeats = @_;
	my $longestRepeat = 0;
	my $longestOffset = 0;
	#my %offsetHash = %$offsetHashRef;
	my $ind = 0;
	for (0 .. ($windowSize - 1)) {
		if ($runningRepeats[$ind] > $longestRepeat) {
			$longestRepeat = $runningRepeats[$ind];
			$longestOffset = $ind;
		}
		$ind++;
	}
	if($longestRepeat >= $minRepeats) {
		#my $currRepeat = $offsetHash{$longestOffset}[0];
		my $currRepeat = $offsetHash->{$longestOffset}->{prevSeq};
		my $currStart = $offsetHash->{$longestOffset}->{countHash}->{$currRepeat}->{startPos};
		my $currType = $offsetHash->{$longestOffset}->{countHash}->{$currRepeat}->{satType};
		my $repLength = $longestRepeat * $windowSize;
		my $stop = ($currStart + $repLength) - 1;
		if ($seqID eq ""){
			$seqID = "$currType\_$currStart\_$stop";
		}
		else {
			$seqID = $seqID . "_$currType\_$currStart\_$stop";
		}
	}
	return $seqID;
}

sub clearCountHash {
	my $offsetHashRef = shift;
	foreach my $offset (keys(%$offsetHashRef)) {
		my $currOffset = $offsetHashRef->{$offset};
	       
		foreach my $mono (keys(%{$currOffset->{countHash}})) {
			$currOffset->{countHash}->{$mono}->{reps} = 0;
		}
	}
}

sub updateTally {
	my $allRepeats = shift;
	my $seq_info = shift;
	my $palHash = shift;
	my $accno = shift;
	
	my $numRepeats = 0;
	my $tmpHash = {};
	
	if ($seq_info) {
		my $longType = "";
		my $longVal = 0;
		my @allInfo = split(/_/, $seq_info);
		$numRepeats = (scalar(@allInfo))/3;
		for(my $i=0; $i < $numRepeats; $i++) {
			my $currRepeat = $allInfo[$i * 3];
			my $start = $allInfo[($i * 3) + 1];
			my $stop = $allInfo[($i * 3) + 2];
			my $length = ($stop - $start) + 1;
			$allRepeats->{$currRepeat}->{totalLoci}++;
			$allRepeats->{$currRepeat}->{totalBases}+= $length;
			if(defined $tmpHash->{$currRepeat}) { $tmpHash->{$currRepeat}++; }
			else { $tmpHash->{$currRepeat} = 1; }
			if ($length > $longVal) { $longVal = $length; $longType = $currRepeat; }
			$palHash->{$accno}->{allLoci} .= "$currRepeat($length) ";
			$palHash->{$accno}->{allLociBases} += $length;
		}
		my @repTypes = keys(%$tmpHash);
		my $numTypes = scalar(@repTypes);
		if ($numTypes > 1) { $allRepeats->{compound}++; }
		elsif(($numTypes == 1) && ($numRepeats > 1)) { $allRepeats->{broken}++; }
		foreach my $type (@repTypes) {
			$allRepeats->{$type}->{readsWithLoci}++;
		}
		$palHash->{$accno}->{longest} = $longType;
	}
	
	return \$allRepeats;
}

sub updatePEtally {
	my $globs = shift;
	my $allRepeats = shift;
	my $PALhash = shift;
	my $accno = shift;
	my $seq1 = shift;
	my $seq2 = shift;
	my $seq1len = shift;
	my $seq2len = shift;
	
	my @seqs = ($seq1, $seq2);
	my $numRepeats = 0;
	my $tmpHash = {};
	my($seqNum, $seq_info);
	my $longBases = 0;
	my $longType = "";
	my $lspan = ""; my $rspan = "";
	
	#my($seq1, $seq2) = split(/[\(\)]/, $seq_info);
	for($seqNum = 0; $seqNum < 2; $seqNum++) {
		$seq_info = $seqs[$seqNum];
		if ($seq_info) {
			$tmpHash = {};
			my @allInfo = split(/_/, $seq_info);
			$numRepeats = (scalar(@allInfo))/3;
			for(my $i=0; $i < $numRepeats; $i++) {
				my $currRepeat = $allInfo[$i * 3];
				my $start = $allInfo[($i * 3) + 1];
				my $stop = $allInfo[($i * 3) + 2];
				my $length = ($stop - $start) + 1;
				if ($length > $longBases) {$longBases = $length; $longType = $currRepeat; }
				$allRepeats->{$currRepeat}->{totalLoci}++;
				$allRepeats->{$currRepeat}->{totalBases}+= $length;
				$PALhash->{$accno}->{allLoci} .= "$currRepeat($length) ";
				$PALhash->{$accno}->{allLociBases} += $length;
				if(defined $tmpHash->{$currRepeat}) { $tmpHash->{$currRepeat}++; }
				else { $tmpHash->{$currRepeat} = 1; }
				if( ($stop == $seq1len) and ($seqNum == 0) ){
					$lspan = $currRepeat;
					$allRepeats->{$currRepeat}->{potentialExtended}++;
					$allRepeats->{allExtended}++;
					$PALhash->{$accno}->{extend} .= "$currRepeat ";
				}
				if ( ($start == 1) and ($seqNum == 1) ){
					$rspan = $currRepeat;
					$allRepeats->{$currRepeat}->{potentialExtended}++;
					$allRepeats->{allExtended}++;
					$PALhash->{$accno}->{extend} .= "$currRepeat ";
				}
				if ($rspan and ($rspan eq $lspan)){
					$allRepeats->{$currRepeat}->{potentialSpan}++;
					$allRepeats->{allSpan}++;
					$PALhash->{$accno}->{span} = $currRepeat;
				}		
			}				
			my @repTypes = keys(%$tmpHash);
			my $numTypes = scalar(@repTypes);
			if ($numTypes > 1) { $allRepeats->{compound}++; }
			elsif(($numTypes == 1) && ($numRepeats > 1)) { $allRepeats->{broken}++; }
			foreach my $type (@repTypes) {
				$allRepeats->{$type}->{readsWithLoci}++;
			}
			$PALhash->{$accno}->{longest} = "$longType";
		}
	}
	return \$allRepeats;
}

sub createMicrosatTallyHash {
	my $globals = shift;
	my @allMonos = (2 .. $globals->{maxMonoSize});
	my @allTypes = ();
	my $tallyHash = {};
	foreach my $mono (@allMonos) {
		my $tmpRef = createAllUniqueRepeats($mono);
		my @tmpList = keys(%$tmpRef);
		push(@allTypes, @tmpList);
	}
	foreach my $type (@allTypes) {
		$tallyHash->{$type}->{totalLoci} = 0;
		$tallyHash->{$type}->{longest} = 0;
		$tallyHash->{$type}->{longestWithPrimers} = 0;
		$tallyHash->{$type}->{totalBases} = 0;
		$tallyHash->{$type}->{readsWithLoci} = 0;
		$tallyHash->{$type}->{basesAmplified} = 0;
		$tallyHash->{$type}->{lociAmplified} = 0;
	}
	$tallyHash->{totalReads} = 0;
	$tallyHash->{totalBases} = 0;
	$tallyHash->{readsWithMicrosat} = 0;
	$tallyHash->{readsWithPrimers} = 0;
	$tallyHash->{broken} = 0;
	$tallyHash->{compound} = 0;
	return $tallyHash;
}

sub createPEmicrosatTallyHash {
	my $globals = shift;
	my @allMonos = (2 .. $globals->{maxMonoSize});
	my @allTypes = ();
	my $tallyHash = {};
	foreach my $mono (@allMonos) {
		my $tmpRef = createAllUniqueRepeats($mono);
		my @tmpList = keys(%$tmpRef);
		push(@allTypes, @tmpList);
	}
	foreach my $type (@allTypes) {
		$tallyHash->{$type}->{totalLoci} = 0;
		$tallyHash->{$type}->{lociWithPrimers} = 0;
		$tallyHash->{$type}->{totalBases} = 0;
		$tallyHash->{$type}->{readsWithLoci} = 0;
		$tallyHash->{$type}->{potentialExtended} = 0;
		$tallyHash->{$type}->{potentialSpan} = 0;
		$tallyHash->{$type}->{extendWithPrimers} = 0;
		$tallyHash->{$type}->{spanWithPrimers} = 0;
	}
	$tallyHash->{totalReads} = 0;
	$tallyHash->{totalBases} = 0;
	$tallyHash->{readsWithMicrosat} = 0;
	$tallyHash->{broken} = 0;
	$tallyHash->{compound} = 0;
	$tallyHash->{allExtended} = 0;
	$tallyHash->{allSpan} = 0;
	return $tallyHash;
}

sub hash_ref {
	my %hash = ();
	return \%hash;
}

sub getMonoList {
	my $globs = shift;
	my $maxSize = shift;
	my @mList = ();
	for my $i (2..$maxSize) {
		push(@mList, $i) unless ( not exists ($globs->{"$i"."merMinReps"}) );
	}
	return @mList;
}

sub ReadConfigFile{
	my $infile = shift;				# input file name string
	my $globs = shift;				# global variables
	my ($var, $value, $maxMono);
	my @range;
	my %sethash; my $setptr; my $rangeptr;
	$maxMono = 2;
	
	my @allVars = ("findPrimers", "platform", "pairedEnd", "input454reads", "inputReadFile", "pairedReadFile", "MicrosatSumOut",
				   "PALsummaryOut", "2merMinReps", "3merMinReps", "4merMinReps", "5merMinReps", "6merMinReps",  "primer3input",
				   "primer3output", "keepPrimer3files", "primer3executable", "prNamePrefix", "PRIMER_TASK", "PRIMER_OPT_SIZE",
				   "PRIMER_MIN_SIZE", "PRIMER_MAX_SIZE", "PRIMER_MAX_NS_ACCEPTED", "pr3ProductSizeRangeMinVal", "pr3ProductSizeRangeMaxVal",
				    "PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_GC_CLAMP", "PRIMER_MAX_END_GC", "PRIMER_MIN_TM", "PRIMER_MAX_TM",
				    "PRIMER_OPT_TM", "PRIMER_PAIR_MAX_DIFF_TM", "PRIMER_TM_FORMULA", "PRIMER_MAX_SELF_ANY", "PRIMER_PAIR_MAX_COMPL_ANY",
				    "PRIMER_MAX_SELF_END", "PRIMER_PAIR_MAX_COMPL_END", "PRIMER_MAX_POLY_X", "PRIMER_LOWERCASE_MASKING",
				    "PRIMER_NUM_RETURN", "PRIMER_MISPRIMING_LIBRARY", "PRIMER_MAX_LIBRARY_MISPRIMING",
				    "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS", "inputFormat");

	open(CONTROL,"<$infile");
	while(<CONTROL>)
	{		
		chomp($_);
		my @lfields= split(/#/, $_);
		$_ = $lfields[0];
		next if(!($_));
		($var, $value) = split(/[\s:=]+/, $_);

		if ($var eq "findPrimers") {
			if(!($value =~ /^[01]$/)) { die "\ninvalid entry \"$value\"for $var. Must be 0 or 1."; }
			$globs->{$var} = $value;
		}
		elsif ($var eq "platform") {
			if(!($value =~ /^(?:(Illumina)|(454))$/)) {die "\ninvalid entry \"$value\"for $var. Must be Illumina or 454 (case sensitive).";}
			$globs->{$var} = $value;
		}
		elsif ($var eq "inputFormat") {
			if(!($value =~ /^(?:(fasta)|(fastq)|(scarf))$/)) {die "\ninvalid entry \"$value\"for $var. Must be fasta, fastq or scarf (case sensitive).";}
			$globs->{$var} = $value;
		}
		elsif ($var eq "pairedEnd") {
			if(!($value =~ /^[01]$/)) { die "\ninvalid entry \"$value\"for $var. Must be 0 or 1."; }
			$globs->{$var} = $value;
		}
		elsif ($var eq "input454reads") { $globs->{$var} = $value; }
		elsif ($var eq "inputReadFile") { $globs->{$var} = $value; }
		elsif ($var eq "pairedReadFile") { $globs->{$var} = $value; }
		elsif ($var eq "MicrosatSumOut") { $globs->{$var} = $value; }
		elsif ($var eq "PALsummaryOut") { $globs->{$var} = $value; }
		elsif ($var =~ /^(?:[2-6])merMinReps$/) {
			if ( !($value =~ /^\d+$/) ) { die "\ninvalid entry \"$value\"for $var. Must be a positive integer.";}
			$var =~ /^([2-6])merMinReps$/;
			my $monoSize = int($1);
			$globs->{$var} = $value;
			if ( ($monoSize > $maxMono) and ($value > 0) ) { $maxMono = $monoSize; }
		}
		elsif ($var eq "primer3input") { $globs->{$var} = $value; }
		elsif ($var eq "primer3output") { $globs->{$var} = $value; }
		elsif ($var eq "keepPrimer3files") {
			if(!($value =~ /^[01]$/)) { die "\ninvalid entry \"$value\"for $var. Must be 0 or 1."; }
			$globs->{$var} = $value;
		}
		elsif ($var eq "primer3executable") { $globs->{$var} = $value; }
		elsif ($var eq "prNamePrefix") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_TASK") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_OPT_SIZE") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MIN_SIZE") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_SIZE") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_NS_ACCEPTED") { $globs->{$var} = $value; }
		elsif ($var eq "pr3ProductSizeRangeMinVal") { $globs->{$var} = $value; }	
		elsif ($var eq "pr3ProductSizeRangeMaxVal") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MIN_GC") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_GC") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_GC_CLAMP") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_END_GC") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MIN_TM") { $globs->{$var} = $value; }		
		elsif ($var eq "PRIMER_MAX_TM") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_OPT_TM") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_PAIR_MAX_DIFF_TM") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_TM_FORMULA") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_SELF_ANY") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_PAIR_MAX_COMPL_ANY") { $globs->{$var} = $value; }	
		elsif ($var eq "PRIMER_MAX_SELF_END") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_PAIR_MAX_COMPL_END") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_POLY_X") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_LOWERCASE_MASKING") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_NUM_RETURN") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MISPRIMING_LIBRARY") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_MAX_LIBRARY_MISPRIMING") { $globs->{$var} = $value; }
		elsif ($var eq "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS") { $globs->{$var} = $value; }
		else { print "\nUnrecognized variable '$var'. Ignoring...\n"; }
	}
	$globs->{maxMonoSize} = $maxMono;
	for my $requiredVar (@allVars) {
		die "Terminated: required variable '$requiredVar' not in configuration file '$infile'\n" if(!(exists $globs->{$requiredVar}));
	}
	print "Configuration File Read.\n";
}

sub array_ref{ my @array; return \@array; }

