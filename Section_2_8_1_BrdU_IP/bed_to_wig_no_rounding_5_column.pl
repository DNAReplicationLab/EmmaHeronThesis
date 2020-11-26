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

my $chr_index    = 0;
my $start_index  = 1;
my $stop_index   = 2;
my $score_index  = 4;
my $step         = 'variable';

open (DATA, "<$infile") or die "didn't open input file $infile $!\n";
open (OUT, ">", "$outfile") or die "didn't open output file $outfile $!\n";

my $current_chr; # current chromosome
my $previous_pos; # previous position to avoid duplicates in wig file

while (my $line = <DATA>)
	{
	
	chomp $line;
	my @data = split /\t/, $line;

	# write definition line if necessary
	if ($data[$chr_index] ne $current_chr)
		{
		# new chromosome, new definition line
		
		if ($step eq 'variable')
			{			
			print OUT 'variableStep chrom=' . $data[$chr_index] . "\n";
			}
		
		# reset the current chromosome
		$current_chr = $data[$chr_index];
		$previous_pos = undef;
		}
	my $score = $data[$score_index];
	my $position = sprintf( "%.0f", $data[$start_index]);
	print OUT "$position\t$score\n" if ($position != $previous_pos);
	$previous_pos = $position;
	}

close (DATA);
close (OUT);
