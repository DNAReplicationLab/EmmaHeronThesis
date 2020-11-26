#!/usr/bin/perl
#
# written by Conrad Nieduszynski, 2011
# questions on the code to: conrad.nieduszynski@nottingham.ac.uk
#
# wig_liftOver_wig.pl
#
# This programme converts wig files to bed, then puts through liftOver and converts back to wig.
# 
# this software is freeware.
# feel free to distribute and modify it.
#


use strict;
use warnings;
use IO::File;
use Getopt::Long;


### Get command line options and initialize values
my (
	$infile,
	$bedfile,
	$outfile,
	$chainfile,
);

# Command line options
GetOptions( 
	"in=s"      => \$infile, # the solexa data file
	'bed=s'     => \$bedfile, # name of output file
	'out=s'     => \$outfile, # name of output file
	'chain=s'   => \$chainfile, # the name of the data, goes into the type field of GFF
);

# initialize variables
my (
	# reusuable variables
	$refseq,
	$fixstart,
	$step,
	$span,
	$out_fh,
); 
my $count = 0;

open (DATA, "<$infile") or die "didn't open input file $infile $!\n";
open (OUT, ">$bedfile") or die "didn't open output file $bedfile $!\n";

while (my $line = <DATA>)
	{
	my $start; # specific line variables
	my $stop;
	my $score;
	
	# The wiggle file can have 3 different formats: BED format, variable step, 
	# and fixed step. We need to determine whether each line is a definition
	# line or a data line, based on the line's contents and/or number of 
	# elements. The definition lines will fill the reusable variables above
	# and help in filling out the specific variables.
	
	## check the line's contents
	$line =~ s/[\r\n]+$//; # strip all line endings
	my @data = split /\s+/, $line;
	
	# a track line
	if ($data[0] =~ /track/i) {
		# not much useable information in here for us
		# but we can use the track name as the type
		foreach (@data) {
			if (/name=(.+)/) {
				last;
			}
		}
	}
	
	# a variable step definition line
	elsif ($data[0] =~ /^variablestep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
		}
		next;
		
	} 
	
	# a fixed step definition line
	elsif ($data[0] =~ /^fixedstep$/i) { 
		foreach (@data) {
			if (/chrom=(\w+)/) {$refseq = $1}
			if (/span=(\w+)/) {$span = $1}
			if (/start=(\w+)/) {$fixstart = $1}
			if (/step=(\w+)/) {$step = $1}
		}
		next;
	} 
	
	# a BED data line
	elsif (scalar @data == 4) {
		$refseq = $data[0];
		$start = $data[1] + 1; # the BED line alone uses 0-based indexing
		$stop = $data[2];
		$score = $data[3];
	} 
	
	# a variable step data line
	elsif (scalar @data == 2) { 
		unless ($refseq) { 
			die "Bad formatting! variable step data but chromosome not defined!\n";
		}
		$start = $data[0];
		if ($span) {
			$stop = $start + $span;
		} 
		else {
			$stop = $start + 1;
		}
		$score = $data[1];
	} 
	
	# a fixed step data line
	elsif (scalar @data == 1) { 
		unless ($refseq) { 
			die "Bad formatting! fixed step data but chromosome not defined!\n";
		}
		unless ($fixstart) { 
			die "Bad formatting! fixed step data but start not defined!\n";
		}
		unless ($step) { 
			die "Bad formatting! fixed step data but step not defined!\n";
		}
		$start = $fixstart;
		$fixstart += $step; # prepare for next round
		if ($span) {
			$stop = $start + $span;
		} 
		else {
			$stop = $start + 1;
		}
		$score = $data[0];
	}

	# format the score value
	my $fscore;
	
	# add the gff data to the data structure
	print OUT $refseq, "\t",
		sprintf( "%.0f", $start), "\t",
		sprintf( "%.0f", $stop), "\t",
		sprintf( "%.3f", $score),
		"\n";
		
	$count++;
	
	}

close (DATA);
close (OUT);
