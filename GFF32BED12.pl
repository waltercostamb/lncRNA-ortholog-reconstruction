#!/usr/bin/perl -w

#Maria Beatriz Walter Costa

#To get the usage, run the script with the help option: perl SCRIPT --help

use strict;
use Data::Dumper;
use Getopt::Long qw(GetOptions);

my $gff_file;
my $type;
my $help;

#Get command line inputs
GetOptions(
        'gff=s' => \$gff_file,
        'type=s' => \$type,
        'help' => \$help,
) or die "Usage: perl $0 --help\n";

#Help page
if ($help){
	print "\nFormat transformation GFF3 to BED12\n";
	print "Maria Beatriz Walter Costa (bia.walter\@gmail.com)\n";
        print "\nUsage: perl $0 --gff GFF3_FILE --type TYPE\n";
        print "\nThis script transforms a GFF3 in BED12 for the user-specified line type (column 3 of the GFF).\n";
	print "If there is a UniProtKB ID, it gets this as feature name, otherwise it is the first ID.";
	print "\nOBS: script is set to no-intron types\n";
        die "\n";
}

#variables $blocks and $blockStartsFIXED are fixed now, because the script does NOT considers introns
#If you want to add introns, work on these variables
my $blocks = 1;
my $blockStartsFIXED = 0;

#Score is set to zero
my $score = 0;

#If not every obligatory input has been defined by the user in the command line, print warning
if ( ! defined($gff_file) || ! defined($type) ) {
	die "\nFor usage: $0 --help\n\n";
}

#Open input file
open ("gff_file", $gff_file) || die "It was not possible to open file $gff_file\n";

#Go through the genes file line by line, detect the line type and prints the correspondent BED12 line
while(<gff_file>){
	chomp;

	#Leave GFF after the genome is reached
	if (/##FASTA/) {
		last;
	}

	#Avoid comment lines
	unless (/^#/) {

		#Split line
		my @row = split (/\t/);
		my $line_type = $row[2];

		#See if line matches the desired type
		if ($line_type eq $type) {


			#Get BED12 fields, following the format specification from https://genome.ucsc.edu/FAQ/FAQformat#format1
			my $chrom = $row[0];
			#BED12 starts are zero based! GFF3 starts are one-based, so convert correctly
			my $chromStart = $row[3] - 1;
			my $chromEnd = $row[4];
			my $pre_name = $row[8];
			my $strand = $row[6];
			my $thickStart = $chromStart;
			my $thickEnd = $chromEnd;
			my $rgb = "-";
			my $blockCount = $blocks;
			my $blockStarts = $blockStartsFIXED;
		
			my $size = $chromEnd - $chromStart + 1;		
			my $blockSizes = $size.","; 

			#Get the name, either UniProtKB ID or ID
			my $name;
			if ($pre_name =~ /UniProtKB/){
				my @tmp1 = split (/UniProtKB\:/, $pre_name);
				my @tmp2 = split (/;/, $tmp1[1]);
				$name = $tmp2[0];			
			} else {
				my @tmp1 = split (/ID=/, $pre_name);
				my @tmp2 = split (/;/, $tmp1[1]);
				$name = $tmp2[0];
			}

			#print "START\n";
			#print "$_\n";
			print "$chrom\t$chromStart\t$chromEnd\t$name\t$score\t$strand\t$thickStart\t$thickEnd\t$rgb\t$blockCount\t$blockSizes\t$blockStarts\n";
			#print "END\n";
		}
	}
}

close gff_file;

