#!/usr/bin/perl

#Maria Beatriz Walter Costa

#This script processes a spliceMAP file and returns to the user a list of orthologus lncRNAs in BED format. It also takes as input a reference lncRNA database in BED format, and returns the orthologs only from this database. In essence, it curates a spliceMap file that also contains transcript starts (TSS) and transcript ends and returns to the user a list of BED-format transcripts that are likely to be transcribed in the cell. It takes into account the splicing mechanism, in which a Donor and an Acceptor sites are both required for the intron splicing. If either one is missing, the intron is retained in the final transcript. It takes into account the transcription initiation model as well, in which a TSS and an upstream promotor are required. In the special case that a TSS is lost in the orthologous species, the whole transcript is regarded as if lost in evolution. This requirement is important to avoid including false orthologs in the results, or orthologous that have a weak transcription signal, due to the TSS being too far from the promoter. If the transcript's end is also missing, the resulting transcript is also regarded as lost in evolution. 

#To get the usage run the script with the help option: perl buildOrthologs.pl --help

use strict;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use Bio::Search::Hit::GenericHit;
use Bio::SearchIO;
use Getopt::Long qw(GetOptions);

my $map_file; #MAP_FILE input.map
my $ref_file; #reference.bed
my $threshold; #SPLICE_map_THRESHOLD
my $reference_species; #REFERENCE_SPECIES
my $blast_file; #BLAST_ADDRESS
my $cutoff; #BLAST_THRESHOLD
my $blastThread; #BLAST thread
my $DB_reg_exp; #Regular expression from the Database you are using (example: ENS for human genes from Gencode)
my $lowerThresholdBlastSequenceInput;
my $upperThresholdBlastSequenceInput;
my $help;

GetOptions(
	'input=s' => \$map_file,
	'db=s' => \$ref_file,
	'map_threshold=s' => \$threshold,
	'reference=s' => \$reference_species,
	'blast_address=s' => \$blast_file,
	'blast_threshold=s' => \$cutoff,
	'blast_thread=s' => \$blastThread,
	'db_reg_exp=s' => \$DB_reg_exp,
	'lower_seq_thresh=s' => \$lowerThresholdBlastSequenceInput,
	'upper_seq_thresh=s' => \$upperThresholdBlastSequenceInput,
	'help' => \$help,
) or die "Usage: $0 --input MAP_FILE --db REF_BED --map_threshold N --reference SPECIES --blast_address FILE --blast_threshold I --blast_thread T --db_reg_exp EXP --lower_seq_thresh L --upper_seq_thresh U\n";

if ($help){
	print "Reconstruction of lncRNA orthologs\n";
	print "Maria Beatriz Walter Costa (bia\@bioinf.uni-leipzig.de)\n\n";
	print "Usage: $0 --input MAP_FILE --db REF_BED --map_threshold N --reference SPECIES --blast_address FILE --blast_threshold I --blast_thread T --db_reg_exp EXP --lower_seq_thresh L --upper_seq_thresh U\n\n";
	print "IMPORTANT: this script requires a high baseline memory. For a 2.6 Mb REF_BED and a 670 Mb MAP_FILE, it required 25 Gb of memory, running on a 32 Gb memory machine.\n";
	print "REQUIRED: you must have a working retrieve-fasta or path to it in your working directory as well as soft links to the genomes (for retrieve-fasta)!\n";
	print "This script uses BLASTN as well as retrieve-fasta internally. BLASTN is used to reconstruct eventual starts and end sites that failed to be reconstructed by SpliceMap. retrieve-fasta is used to retrieve temporary FASTA sequences from temporary BED files.\n\n";
	print "The output of this script is in BED12 format, with a 13th row indicating the species. A BED line starting with # indicates a BLASTN reconstruction. \n\n";
	print "All parameters are required for running the script:\n";
	print "--input: site MAP_FILE from the SpliceMap Pipeline\n";
	print "--db: REF_BED, or BED12 reference lncRNA database\n";
	print "--map_threshold: splice map threshold for excluding unlikely real splice sites, N can be: undef, -3, 0, 3, etc. (recommended stringent threshold: 3)\n";
        print "--reference: alias of the reference species, example: hg38, ponAbe2, rheMac3, mm10, etc. as in the input MAP_FILE\n";
	print "--blast_address: a file containing two rows divided by a tab: first row with species alias (hg38, mm10, etc.) and second row with full path to index files for BLASTN\n"; 
	print "--blast_threshold: e-value I threshold for filtering BLASTN hits\n"; 
	print "--blast_thread: T number of threads to run BLASTN\n";
	print "--db_reg_exp: regular expression EXP used to select lines from the --input MAP_FILE based on the --db REF_BED (e.g. ENST for human Gencode transcripts, ENSMUS for mouse, etc.)\n";
	print "--lower_seq_thresh: lower sequence length threshold L of input sequence to BLAST (anything shorter than L is not submitted to BLAST and marked as invalid reconstruction)\n";
	print "--upper_seq_thresh: upper sequence length threshold U of input sequence to BLAST (anything longer than U is not submitted to BLAST and marked as invalid reconstruction)\n";
	die "\n";
}

open ("map_file", $map_file) || die "It was not possible to open file $map_file\n";
open ("ref_file", $ref_file) || die "It was not possible to open file $ref_file\n";
open ("blast_file", $blast_file) || die "It was not possible to open file $blast_file\n";

#Saving blast address paths to internal Perl variables
my %blastAddresses;
while(<blast_file>){
	chomp;
	my @tmp = split /\t/;
	$blastAddresses{$tmp[0]} = $tmp[1];
}
close blast_file;

my %transcript2ids;
my %id2mapLines;
#Processing each line of the SpliceMap file to store the transcript_names and the map_ids into hashes: %transcript2ids and %ids2lines
while (<map_file>){
	chomp;
	my @line = split /\t/;
	#All transcript_names of the current map line go to @transcriptNames
	my @transcriptNames = split /,/, $line[12];
	#Foreach of the transcript_names of the current map line, store map_ids and lines
	foreach (@transcriptNames){
		#If transcript_name equals our IDs of interest (ex: ENST for ENSEMBL transcript IDs, store the map_ids)
		if (/$DB_reg_exp/){
			if(exists $transcript2ids{$_}){
				my @a = @{$transcript2ids{$_}};
				push @a, $line[0];
				#Hash key is the transcript_name and hash value is an array with the map_ids
				@transcript2ids{$_} = [@a];
			} else {
				my @a = ($line[0]);
				@transcript2ids{$_} = [@a];
			}
		}
	}

	#For the current map_id, store the line into hash %id2mapLines
	my $idKey = $line[0];
	if(exists $id2mapLines{$idKey}){
		my @a = @{$id2mapLines{$idKey}};
		push @a, [@line];
		#Hash key is the map_id and hash value is an array with the lines corresponding to this map_id
		@id2mapLines{$idKey} = [@a];
	} else {
		my @a = [@line];
		@id2mapLines{$idKey} = [@a];
	}


}

#Going through the BED-format reference dataset, for example GENCODE dataset
#Working through every line of the reference dataset, to work with each transcript_id
while (<ref_file>){	
	chomp;
	my @line = split /\t/;
	my $transcript = $line[3];
	my @ReferenceLimitBlocks = (); #For eventual reconstruction of orthologs' global Start and End
	my %Guide;
	#The following arrays contain all the map_lines corresponding to the current transcript_id
	my @mapLines = transcript2mapLines($transcript);
	my @sortedMapLines = sortMapLinesBySpeciesContigs(\@mapLines);

	unless (@sortedMapLines == 0){
		#We want in the following block to separate the big array of @sortedMapLines into smaller arrays according to the species+contig+strand
		#Setting the first reference of: species, contig and strand to variables below. They refer to the first element of @sortedMapLines
		my $species = $sortedMapLines[0][2];
		my $contig = $sortedMapLines[0][3];
		my $strand = $sortedMapLines[0][4];
		my @pre_SDAN = (); #@pre_SDAN refers to the ortholog BED line
		my @fixed_pre_SDAN;
		my %DAN;
		my $BEDline;
		my @ref_BED;
		my $reconstructedSignal;

		#Storing every map_line that belongs to this contig to array @pre_SDAN 
		foreach(@sortedMapLines){
			#If the elements, match, store the map_line to @pre_SDAN
			if ($species eq $_->[2] && $contig eq $_->[3] && $strand eq $_->[4]){
				push @pre_SDAN, [@{$_}];
				$reconstructedSignal = "no";
			#When the current line do not match to the reference, we go to else, print the @SDAN and re-start the reference to the current one
			} else {
				#Greedy algorithm
				@fixed_pre_SDAN = fixMinusStrand(\@pre_SDAN);
				%DAN = prepareBEDline(\@fixed_pre_SDAN,$transcript);

				#The following if transforms invalid into valid transcripts by getting the start and end based on the reference
				if($DAN{'status'} eq 'invalid'){
					$Guide{'guide'} = $DAN{'guide'};
					$Guide{'contig'} = $DAN{'contig'};
					$Guide{'strand'} = $DAN{'strand'};
					$Guide{'start'} = $DAN{'start'}; 
					$Guide{'end'} = $DAN{'end'};
					$Guide{'transcript'} = $DAN{'transcript'};

					my %missingSites = @{$DAN{'missing'}};
					my @reconstructed_pre_SDAN = @fixed_pre_SDAN;
					if (exists $missingSites{'start'}){
						my %referenceSequence = getSequence($ref_BED['start'], $reference_species);
						#sub args: (i) query sequence, (ii) reference for search genome, (iii) search genome, (iv) evalue cutoff for blast filter and (v) guide for best hit choice
						my %BLASTresult = BLAST2bestHit ($referenceSequence{'sequence'},$species,$cutoff,\%Guide,$blastThread,$lowerThresholdBlastSequenceInput,$upperThresholdBlastSequenceInput);

						if (%BLASTresult){
							@reconstructed_pre_SDAN = addSDAN (\@{$DAN{'missing'}},\@reconstructed_pre_SDAN,\%BLASTresult);
							%DAN = prepareBEDline (\@reconstructed_pre_SDAN,$transcript);
							$reconstructedSignal = "yes";
						}
					}

					if (exists $missingSites{'end'}){
						my %referenceSequence = getSequence($ref_BED['end'], $reference_species);
						#sub args: (i) query sequence, (ii) reference for search genome, (iii) search genome, (iv) evalue cutoff for blast filter and (v) guide for best hit choice
						my %BLASTresult = BLAST2bestHit ($referenceSequence{'sequence'},$species,$cutoff,\%Guide,$blastThread,$lowerThresholdBlastSequenceInput,$upperThresholdBlastSequenceInput);

						if (%BLASTresult){
							@reconstructed_pre_SDAN = addSDAN (\@{$DAN{'missing'}},\@reconstructed_pre_SDAN,\%BLASTresult);
							%DAN = prepareBEDline (\@reconstructed_pre_SDAN,$transcript);
							$reconstructedSignal = "yes";
						}
					}
				}

				writeBEDline(\%DAN,'printYes',$reconstructedSignal);

				#Get the reference start and end blocks for eventual future orthologs starts/ends reconstruction
				#The reconstruction will be different if the reference has only one exon block or more
				if($DAN{'species'} eq $reference_species){
					my @blockNumber = @{$DAN{'blockSizes'}};

					if (@blockNumber == 1){ #If reference has only one block, stores its BED line in @ref_BED
						push @ReferenceLimitBlocks, $DAN{'blockSizes'};
						my $curr_ref = writeBEDline(\%DAN,'printNo',$reconstructedSignal);
						push @ref_BED, $curr_ref;
					} elsif (@blockNumber > 1) { #If reference has more than one block, takes the first and end blocks, constructs BED lines for each and store it in @ref_BED
						push @ReferenceLimitBlocks, $DAN{'blockSizes'}[0];
						push @ReferenceLimitBlocks, $DAN{'blockSizes'}[-1];

						my %subsetDAN = extremesBEDlines (\%DAN); 

						my $curr_ref = writeBEDline(\%{$subsetDAN{'startBlock'}},'printNo',$reconstructedSignal);
						push @ref_BED, $curr_ref;
						my $curr_ref = writeBEDline(\%{$subsetDAN{'endBlock'}},'printNo',$reconstructedSignal);
						push @ref_BED, $curr_ref;
					}
				}
				#Lines below reset the variables so that the if above re-starts with the new species or contig
				@pre_SDAN = ();
				push @pre_SDAN, [@{$_}];
				$species = $_->[2];
				$contig = $_->[3];
				$strand = $_->[4];
			}
		}
		#Run the greedy algorithm to the last element and print it
		@fixed_pre_SDAN = fixMinusStrand(\@pre_SDAN);
		%DAN = prepareBEDline(\@fixed_pre_SDAN,$transcript); 

		#The following if transforms invalid into valid transcripts by getting the start and end based on the reference
		if($DAN{'status'} eq 'invalid'){
			$Guide{'guide'} = $DAN{'guide'};
			$Guide{'contig'} = $DAN{'contig'};
			$Guide{'strand'} = $DAN{'strand'};
			$Guide{'start'} = $DAN{'start'}; 
			$Guide{'end'} = $DAN{'end'};
			$Guide{'transcript'} = $DAN{'transcript'};

			my %missingSites = @{$DAN{'missing'}};
			my @reconstructed_pre_SDAN = @fixed_pre_SDAN;
			if (exists $missingSites{'start'}){
				my %referenceSequence = getSequence($ref_BED['start'], $reference_species);
				#sub args: (i) query seq, (ii) ref for search genome, (iii) search genome, (iv) evalue cutoff for blast filter, (v) guide for best hit choice
				my %BLASTresult = BLAST2bestHit ($referenceSequence{'sequence'},$species,$cutoff,\%Guide,$blastThread,$lowerThresholdBlastSequenceInput,$upperThresholdBlastSequenceInput);

				if (%BLASTresult){
					@reconstructed_pre_SDAN = addSDAN (\@{$DAN{'missing'}},\@reconstructed_pre_SDAN,\%BLASTresult);
					%DAN = prepareBEDline (\@reconstructed_pre_SDAN,$transcript);
					$reconstructedSignal = "yes";
				}
			}

			if (exists $missingSites{'end'}){
				my %referenceSequence = getSequence($ref_BED['end'], $reference_species);
				#sub args: (i) query seq, (ii) ref for search genome, (iii) search genome, (iv) evalue cutoff for blast filter, (v) guide for best hit choice
				my %BLASTresult = BLAST2bestHit ($referenceSequence{'sequence'},$species,$cutoff,\%Guide,$blastThread,$lowerThresholdBlastSequenceInput,$upperThresholdBlastSequenceInput);

				if (%BLASTresult){
					@reconstructed_pre_SDAN = addSDAN (\@{$DAN{'missing'}},\@reconstructed_pre_SDAN,\%BLASTresult);
					%DAN = prepareBEDline (\@reconstructed_pre_SDAN,$transcript);
					$reconstructedSignal = "yes";
				}
			}
		}
		writeBEDline(\%DAN,'printYes',$reconstructedSignal);
	}
}

close ref_file;
close map_file;

################### Subroutines #####################

#This subroutine receives as (input) the transcript_id, (i) gets all the map_ids of the reference referring to this transcript_is, and returns as (output) all the lines referring to the map_ids, in array @mapLines 
sub transcript2mapLines {
	#@ids is an array with all the map_ids related to the current transcript_id
	my $transcript = $_[0];
	my @mapLines;
	unless ($transcript2ids{$transcript}){
		@mapLines = ();
		return @mapLines;
	}
	my @ids = @{$transcript2ids{$transcript}};

	#Retrieving all map_lines for the current transcript_id
	foreach(@ids){
		push (@mapLines, @{$id2mapLines{$_}});
	}
	return @mapLines;
}

#This subroutine will crescent sort all mapLines according to the criteria in hierarchical order: (i) species, (ii) chromosome/contig, (iii) strand and (iv) position
sub sortMapLinesBySpeciesContigs {
	my @mapLines = @{$_[0]};

	#The first sorted species must always be the reference! Because the reference is used to build reference blocks for posterior ortholog missing block reconstruction

	my @firstBlockSort_onlyRef;
	my @secondBlockSort_others;

	#Separate mapLines into reference lines and other species lines
	foreach (@mapLines){
		if ($_->[2] eq "$reference_species"){
			push (@firstBlockSort_onlyRef, $_);	
		} else {
			push (@secondBlockSort_others, $_);
		}
	}

	#Sort map lines according to blocks
	my @sortedFirstBlock = sort { ($a->[3] cmp $b->[3]) || ($a->[4] cmp $b->[4]) || ($a->[5] <=> $b->[5]) } @firstBlockSort_onlyRef;
	my @sortedSecondBlock = sort { ($a->[2] cmp $b->[2]) || ($a->[3] cmp $b->[3]) || ($a->[4] cmp $b->[4]) || ($a->[5] <=> $b->[5]) } @secondBlockSort_others;

	#Combine all map lines in sorted order, with reference species first
	my @sortedMapLines = @sortedFirstBlock;
	push @sortedMapLines, @sortedSecondBlock;

	return @sortedMapLines;
}

#If strand is negative, we simply substitute N->S, A->D, D->A, S->N, so that S to N increases the coordinates  
sub fixMinusStrand {
	my @pre_SDAN = @{$_[0]};
	if ($pre_SDAN[0][4] eq '-'){
		foreach(@pre_SDAN){	
			if ($_->[7] =~ s/S/N/){
				next;
			} elsif ($_->[7] =~ s/D/A/){
				next;
			} elsif ($_->[7] =~ s/A/D/){
				next;
			} elsif ($_->[7] =~ s/N/S/){
				next;
			}
		}
	}
	return @pre_SDAN;
}

#This subroutine does the greedy algorithm to build the transcript
sub prepareBEDline{
	my @fixed_pre_SDAN = @{$_[0]};
	my $transcript = $_[1];
	my $countS = 0;
	my $countN = 0;
	my $start = 0;
	my $end = 0;
	my @blockStarts;
	my @blockSizes;
	my $guide = 'undef';
	my @missing;

	#Check if this data (species, contig and strand) has only one Start and only one End, if not, report an error for the current data
	foreach (@fixed_pre_SDAN){
			if ($_->[7] =~ /S/){
				$countS++;
				$start = $_->[5];
				$guide = ();
				$guide = $start;
			} elsif ($_->[7] =~ /N/){
				$countN++;
				$end = $_->[5];
				if($guide eq 'undef'){
					$guide = ();
					$guide = $end;
				}
			} elsif ($_->[7] =~ /D/){
				if($guide eq 'undef'){
					$guide = ();
					$guide = $end;
				}
			} elsif ($_->[7] =~ /A/){
				if($guide eq 'undef'){
					$guide = ();
					$guide = $end;
				}
			}	
	}

	if ($start == 0){
		push @missing, 'start';
	} elsif ($end == 0) {
		push @missing, 'end';
	}

	#If the data has start and end, build the transcript structure, if missing, output an error report
	if ($countS == 1 && $countN == 1 && $fixed_pre_SDAN[0][7] eq 'S' && $fixed_pre_SDAN[-1][7]){
		push @blockStarts, $start;
		my @blockSizes;
		my $lastSite;
		my $lastPosition;
		foreach (@fixed_pre_SDAN){
				if ($_->[7] =~ /S/){
					#chromStart was already set the block before, and blockStarts has already been initialized
					$lastSite = 'S';
					$lastPosition = $_->[5];
				} elsif ($_->[7] =~ /D/){
					#Independent on the last position, we want to set the last one to a 'D' after seeing it, since it will only be valid if the next one is an 'A'
					$lastSite = 'D';
					$lastPosition = $_->[5];
				} elsif ($_->[7] =~ /A/){
					if($lastSite =~ /A/){
						#Do nothing, since the first seen 'A' yields the longest transcript
					} elsif ($lastSite =~ /S/){
						#Do nothing, since this 'A' is not valid, since it is missing its previous 'D' pair
					} elsif ($lastSite =~ /D/) {
						#This is now a valid intron, since it has a pair 'A' preceeded by a 'D'
						#The last 'D' was the end of the last exon
						my $lengthLastExon = $lastPosition - $blockStarts[-1];
						push @blockSizes, $lengthLastExon;
						#This intron is valid, therefore it means that the current position is the start of a new exon
						push @blockStarts, $_->[5];
						$lastSite = 'A';
						$lastPosition = $_->[5];
					}
				} elsif ($_->[7] =~ /N/){
					#The End has already been processed on the block before
					if ($lastSite =~ /S/){
						#The current transcript structure has only one exon
						my $lengthExon = $end - $start + 1;
						push @blockSizes, $lengthExon;
					} elsif ($lastSite =~ /A/){
						#The last 'A' was the start of the current block
						my $lengthExon = $end - $lastPosition + 1;
						push @blockSizes, $lengthExon;
					} elsif ($lastSite =~ /D/){
						#The last 'D' is not valid, since it does not have its 'A' pair, so I have to keep track of the last exon, starting at the last @blockStarts position
						my $lengthExon = $end - $blockStarts[-1] + 1;
						push @blockSizes, $lengthExon;
					}

				}
		}
		
		$blockStarts[0] -= 1;
		#The @blockStarts we calculated on the previous block of the greedy algorithm is in reference to the contig, not the chromStart as in the BED file, so in the next block, we will make the blockStarts in reference to chromStart
		foreach(@blockStarts){
			$_ = $_ - $start + 1;
		}
		my %BEDline;
		$BEDline{'status'} = 'valid';
		$BEDline{'start'} = $start;
		$BEDline{'end'} = $end;
		$BEDline{'blockStarts'} = \@blockStarts;
		$BEDline{'blockSizes'} = \@blockSizes;
		$BEDline{'species'} = $fixed_pre_SDAN[0][2];
		$BEDline{'contig'} = $fixed_pre_SDAN[0][3];
		$BEDline{'strand'} = $fixed_pre_SDAN[0][4];
		$BEDline{'transcript'} = $transcript;
		$BEDline{'guide'} = $guide;
		return %BEDline;
	} else {
		my %BEDline;
		$BEDline{'status'} = 'invalid';
		$BEDline{'species'} = $fixed_pre_SDAN[0][2];
		$BEDline{'contig'} = $fixed_pre_SDAN[0][3];
		$BEDline{'strand'} = $fixed_pre_SDAN[0][4];
		$BEDline{'transcript'} = $transcript;
		$BEDline{'guide'} = $guide;
		$BEDline{'missing'} = \@missing;
		return %BEDline;
	}
}



#Subroutine writes and prints the BED line 
sub writeBEDline {
	my %BEDline = %{$_[0]};
	my $outputFormat = $_[1];
	my $reconstructedSignal = $_[2];

	if ($BEDline{'status'} eq 'valid'){
		my @tempStarts = @{$BEDline{'blockStarts'}};
		my @tempSizes = @{$BEDline{'blockSizes'}};
		my $blockNumber = @tempStarts;
		my $blockStart = join (",", @tempStarts);
		my $blockSizes = join (",", @tempSizes) . ",";

		if($outputFormat eq 'printYes'){
			if ($reconstructedSignal eq 'yes'){
				print "#$BEDline{'contig'}\t$BEDline{'start'}\t$BEDline{'end'}\t$BEDline{'transcript'}\t500\t$BEDline{'strand'}\t$BEDline{'start'}\t$BEDline{'end'}\t.\t$blockNumber\t$blockSizes\t$blockStart\t$BEDline{'species'}\n";
			} elsif ($reconstructedSignal eq 'no') {
				print "$BEDline{'contig'}\t$BEDline{'start'}\t$BEDline{'end'}\t$BEDline{'transcript'}\t500\t$BEDline{'strand'}\t$BEDline{'start'}\t$BEDline{'end'}\t.\t$blockNumber\t$blockSizes\t$blockStart\t$BEDline{'species'}\n";
			}
		} elsif ($outputFormat eq 'printNo'){
			my $BEDline;
			$BEDline = "$BEDline{'contig'}\t$BEDline{'start'}\t$BEDline{'end'}\t$BEDline{'transcript'}\t500\t$BEDline{'strand'}\t$BEDline{'start'}\t$BEDline{'end'}\t.\t$blockNumber\t$blockSizes\t$blockStart\t$BEDline{'species'}\n";
			return $BEDline;
		}
	} else {
		if($outputFormat eq 'printYes'){
			print "#invalidTranscript\t$BEDline{'species'}\t$BEDline{'transcript'}\t$BEDline{'contig'}\t$BEDline{'strand'}\n";
		} elsif ($outputFormat eq 'printNo'){
			my $BEDline;
			$BEDline = "#$BEDline{'species'}\t$BEDline{'transcript'}\t$BEDline{'contig'}\t$BEDline{'strand'}\n";
			return $BEDline;
		}
	}
}

#Subroutine returns a hash with (i) sequence from a BED line and (ii) start position
#Call: my $referenceSequence = getSequence($ref_BED, $reference_species);
sub getSequence {
	my $ref_BED = $_[0];
	my @line = split /\t/, $ref_BED;
	my $reference_species = $_[1];
	my $tmp = File::Temp->new(TEMPLATE => 'bed_XXXX', SUFFIX => '.tmp');
	my $fn = $tmp->filename;
	print $tmp $ref_BED;

#print Dumper $ref_BED;
#print Dumper $reference_species;

	my %myFasta;
	my $rf_temp = `./retrieve-fasta $reference_species $fn`;
#print "Status of retrieve fasta: $? (-1 eq 'failed to execute: https://stackoverflow.com/questions/11451680/what-is-the-meaning-of-the-built-in-variable-in-perl')\n";
	my $rf_seq;

	if($? == -1){
		print "ERROR: retrieve-fasta failed to be executed!\n";
		die;
	}

#print $tmp;
#print Dumper $rf_temp;
#die;

	if ($rf_temp =~ /\n(\w+\n)/) {
		$rf_seq = $1;
	}

	$myFasta{'sequence'} = $rf_seq;
	$myFasta{'contig'} = $line[0];
	$myFasta{'start'} = $line[1];
	$myFasta{'strand'} = $line[5];
	close $tmp;

	return %myFasta;
}

#Subroutine receives (i) query sequence and guide, (ii) reference for genome to search in, (iii) genome to search in, (iv) cutoff for blast hits. It returns a hash with BLAST best hit (i) contig, (ii) start, (iii) end, (iv) sequence. Best hit was chosen from a pool of BLAST hits within the e-value threshold and that was within the shortest distance to reference guide
#Hash %blastAddresses is required for this subroutine (as a global structure to the main program)
sub BLAST2bestHit {
	my $referenceSequence = $_[0];
	my $species = $_[1];
	#Cutoff for BLAST results
	my $cutoff = $_[2];
	my %guide = %{$_[3]}; 
	my $blastThread = $_[4];
	my $lowerThresholdBlastSequenceInput = $_[5];
	my $upperThresholdBlastSequenceInput = $_[6];
	my $length = length($referenceSequence);

	#Condition created to restrict BLAST runs, if the input is nonsense (e. g. input sequence too small, or too big)
	unless ($length <= $lowerThresholdBlastSequenceInput || $length >= $upperThresholdBlastSequenceInput){
		#Creating a temp file to contain the BLAST's query sequence 
		my $tmp = File::Temp->new(TEMPLATE => 'seq_XXXX', SUFFIX => '.tmp');
		my $fn = $tmp->filename;
		print $tmp $referenceSequence;

#		print Dumper $blastAddresses{$species};
#		print Dumper $referenceSequence;

		#Running BLASTN with the parameters givn to the subroutine (query sequence and database/species)
		my $blastResult = `blastn -query $tmp -db $blastAddresses{$species} -num_threads=$blastThread`;

		#Pretending that variable $blastResult is a filehandle $fh
		my $fh;
		open ($fh, "<", \$blastResult);

		#Passing the filehandle $fh to function 'new' of module Bio::SearchIO
		#$object is the object containing the blast results
		my $object = Bio::SearchIO->new(-format => 'blast', -fh => $fh, -signif => $cutoff);
		my $start;
		my $end;
		my $strand;
		my $contig;
		my $seq;
		my $markerDistance = 'undef';
		my $marker = 1;
		my $currDistance;
		my $curr_strand;
		my $comparator_strand;
		my $counter = 0;

		#Loop goes through $result from one to the next
		while ( my $result = $object->next_result() ) {
			#Loop gets the hits

			$counter++;

			while ( my $hit = $result->next_hit ) {
				while ( my $hsp = $hit->next_hsp ){

					my $evalue = $hsp->evalue();

					#Filter hits that are lower than the cutoff
					if($evalue < $cutoff){

						$curr_strand = $hsp->strand('hit');
						if($curr_strand eq '1'){
							$comparator_strand = '+';
						} elsif ($curr_strand eq '-1'){
							$comparator_strand = '-';
						}

						if ($markerDistance eq 'undef'){
							#If first hit, it will automatically be the best hit
							#$marker=0 marks if I want to substitue the best hit 
							$marker = 0;	
							#If same contig and strand!
							if ($hit->name() eq $guide{'contig'} && $comparator_strand eq $guide{'strand'}){
								$markerDistance = abs($hsp->start('hit') - $guide{'guide'});
							}
						}	

						#If same contig and strand!
						if ($hit->name() eq $guide{'contig'} && $comparator_strand eq $guide{'strand'}){
							$currDistance = abs($hsp->start('hit') - $guide{'guide'});
							if ($currDistance < $markerDistance){
								$marker = 0;
								$markerDistance = $currDistance;
							}
						}

						#Starts the best hit with the top hit and substitutes it IF evalue keeps under cutoff AND distance between start and guide is lower than the old best hit
						if ($marker == 0) {
							$start = $hsp->start('hit');
							$contig = $hit->name();
							$end = $hsp->end('hit');
							$strand = $hsp->strand('hit');
							$seq = $hsp->hit_string();
							$marker = 1;
						}
					}
				}
			}
		}

		my $newStrand;
		if($strand eq '1'){
			$newStrand = '+';
		} elsif ($strand eq '-1'){
			$newStrand = '-';
		}

		$seq =~ s/\-//g;
		my %hit;

		if($guide{'contig'} eq $contig && $guide{'strand'} eq $newStrand){
			%hit = (
			   	   'contig' => $contig,
				    'start' => $start,
				      'end' => $end,
				   'strand' => $newStrand,
				 'sequence' => $seq,
			);
		}

		close $fh;
		close $tmp;

		return %hit;

	}
}

#my @full_pre_SDAN = addSDAN ($_,\@pre_SDAN,\%BLASTresult);
#Subroutine checks which site is missing and adds it to @fixed_pre_SDAN. It resturns the full @fixed_pre_SDAN. 
sub addSDAN {
	my @missingSite = @{$_[0]};
	my @fixed_pre_SDAN = @{$_[1]};
	my %BLASTresult = %{$_[2]};
	my @missingLine;

	unless($BLASTresult{'contig'} eq ''){
		foreach (@missingSite){
			my $siteType;
			my $position;
			if ($_ eq 'start'){ #If the missing site is a start, the BLAST coordinate is formatted to BED coordinate
				$siteType = 'S';
				$position = $BLASTresult{$_} - 1;
				@missingLine = ("undef","*","$fixed_pre_SDAN[0][2]","$BLASTresult{'contig'}","$BLASTresult{'strand'}","$position","B","$siteType","0","?","?","?","?","undef","no");
				unshift @fixed_pre_SDAN,\@missingLine;
			} elsif ($_ eq 'end'){ #If the missing site is a start, the BLAST coordinate is formatted to BED coordinate
				$siteType = 'N';
				$position = $BLASTresult{$_} + 1;
				@missingLine = ("undef","*","$fixed_pre_SDAN[0][2]","$BLASTresult{'contig'}","$BLASTresult{'strand'}","$position","B","$siteType","0","?","?","?","?","undef","no");
				push @fixed_pre_SDAN,\@missingLine;
			}	
		}
	}

	return @fixed_pre_SDAN;
}

#This subroutine takes as input a processed and complete %DAN hash and returns a hash with a subset of %DAN
sub extremesBEDlines {
	my %DAN = %{$_[0]};
	my %subsetDAN;

	my $endStartBlock = $DAN{'start'} + $DAN{'blockSizes'}[0];
	my @blockStarts = ("0");
	my @blockSizesStartBlock = $DAN{'blockSizes'}[0];
	my %startBlock = (
		"status" => $DAN{'status'},
		"guide" => $DAN{'start'},
		"transcript" => $DAN{'transcript'},
		"species" => $DAN{'species'},
		"strand" => $DAN{'strand'},
		"contig" => $DAN{'contig'},
		"blockSizes" => \@blockSizesStartBlock,
		"blockStarts" => \@blockStarts,
		"end" => $endStartBlock,
		"start" => $DAN{'start'},
	);

	my $startEndBlock = $DAN{'end'} - $DAN{'blockSizes'}[-1];
	my @blockSizesEndBlock = $DAN{'blockSizes'}[-1];

	my %endBlock = (
		"status" => $DAN{'status'},
		"guide" => $startEndBlock,
		"transcript" => $DAN{'transcript'},
		"species" => $DAN{'species'},
		"strand" => $DAN{'strand'},
		"contig" => $DAN{'contig'},
		"blockSizes" => \@blockSizesEndBlock,
		"blockStarts" => \@blockStarts,
		"end" => $DAN{'end'},
		"start" => $startEndBlock,
	);

	$subsetDAN{'startBlock'} = \%startBlock;
	$subsetDAN{'endBlock'} = \%endBlock;
	
	return %subsetDAN;
}
