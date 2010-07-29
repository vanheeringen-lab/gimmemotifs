#!@WHICHPERL@ -w
#
# $Id $
# $Log $
#
# FILE: priority2meme.pl
# AUTHOR: Philip Machanick
# CREATE DATE: 5/15/2009
# DESCRIPTION: Convert a file containing a TFBS motif matrix from 
# PRIORITY (Hartemink) to MEME output format.

# requires a background file in MEME format as command line arg
# Input format (# at end of line is a comment, not part of actual data):
# ###blank lines etc. skipped up to here
# Transcription factor: <TF name> ### <TF name> contains alphanumerics, _ or -
# ###skip any lines until one containing
# Number of sequences: <number>
# ###skip any lines until one containing
# Motif length: <number>
###skip any lines until one starting
# Phi:  ### followed by a line starting with 1 ending with a number
### that number must be = Motif length: <number>
## next 4 lines should be of form
# <base> <number> <number> ... # where base is one of a, c, g, t
## rest of file is ignored

use warnings;
use strict;
use Data::Dumper;

# Set up global variables. Assume uniform.
my $pseudo_total = 0; 		# default total pseudocounts
my @bases  = ("A", "C", "G", "T");
my $num_bases  = 4;
my %bg = ( 'A' => .25, 'C' => .25, 'G' => .25, 'T' => .25,);
my $STRANDS = 2;
my $use_numbers = 0;
my $num_seqs; # in data file
my $motif_length; # in data file

my $usage = "USAGE: priority2meme [options] <matrix file>

  Options: 
        -bg <background file>           MEME-format background file
	-numbers       		  	use numbers instead of strings as motif 
					IDs
	-pseudo <total pseudocounts>	add <total pseudocounts> times f_a to 
                                        each freq default: $pseudo_total\n
  Read a PRIORITY (Hartemink) matrix file and convert to MEME format.\n\n";

# Process command line arguments.
if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}

my $next_arg;
my $bg_file;
my %skips;
while (scalar(@ARGV) > 1) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-numbers") {
    $use_numbers = 1;
  } elsif ($next_arg eq "-bg") {
    $bg_file = shift(@ARGV);
  } elsif ($next_arg eq "-pseudo") {
    $pseudo_total = shift(@ARGV);
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}
my $matrix_file;
($matrix_file) = @ARGV;

my $a;
my $f;
my $total_bg;
# read the background file if defined, otherwise leave %bg on defaults.
# does no checking on contents of file
if  ($bg_file) {
    open(BG_FILE, $bg_file) || die("Can't open $bg_file.\n");
    $total_bg = 0;
    while (<BG_FILE>) {
	next if (/^#/);      # skip comments
	($a, $f) = split;
	# stop if we go over order 0 background
	last if length($a) > 1;
	$bg{uc $a} = $f;
	$total_bg += $f;
    }
# make sure they sum to 1
    foreach my $key (keys %bg) {
	$bg{$key} /= $total_bg;
    }
    close BG_FILE;
}

# Open the matrix file for reading.
open(MATRIX_FILE, "<$matrix_file") || die("Can't open $matrix_file.\n");

# Print the MEME header.
print("MEME version 3.0\n\n");
print("ALPHABET= ACGT\n\n");
print("strands: + -\n\n");

my $bg_source = ($bg_file) ? "`$bg_file'" : "defaults";
print(
  "Background letter frequencies " .
  "(from $bg_source with add-one prior applied):\n"
);
printf("A %f C %f G %f T %f\n\n",  $bg{"A"}, $bg{"C"}, $bg{"G"}, $bg{"T"});

# Read the input file.
my $num_skipped = 0;
my (%motifs, $motif_name, %motif_freqs, $tf_name);
# skip blank lines and get TF name
while (<MATRIX_FILE>) {
  next if (/^\s*$/); 		# skip blank lines
  if (/^\s*Transcription factor:\s*([\w\-]+)/) {
  	$tf_name = $1;
  	last;
  } else {
  	die "`$_': should contain TF name";
  }
}
die "premature end of file, no TF name" unless (defined($tf_name));

#get number of sequences
while (<MATRIX_FILE>) {
  next if (/^\s*$/); 		# skip blank lines
  if (/\s*Number of sequences:\s*(\d+)/) {
  	$num_seqs = $1;
  	last;
  }
}
die "Number of sequences: line missing in motif file" if (!defined($num_seqs));
#get motif length
while (<MATRIX_FILE>) {
  next if (/^\s*$/); 		# skip blank lines
  if (/^\s*Motif length:\s*(\d+)\s*/) {
  	$motif_length = $1;
  	last;
  }
}
die "Motif length: line missing in motif file" if (!defined($motif_length));

my $second_last = 0;
# skip to 'Phi:' line and line containing numbers of sequence positions
while (<MATRIX_FILE>) {
  next if (/^\s*$/); 		# skip blank lines
  if ($second_last) {
  	if (/^\s*1\s*.*(\d+)\s*$/) {	# this line should be 1 2 ... motif length
  		die "motif length in motif file inconsistent" if ($motif_length != $1);
  		last;
  	}
  }
  if (/^\s*Phi/) {	# stop one past Phi line
  	$second_last = 1;
  }
}

while (<MATRIX_FILE>) {
  last if (/^\s*$/); 		# blank line at end
  # Split the line into base name and everything else.
  my ($base, @data) = split;
  if ($base =~ /^\s*([acgtACGT])\s*$/) {
    # frequency line
    $base = $1;
    $motif_freqs{uc($base)} = [ @data ];
    die "line `$base' wrong length ".scalar(@data) unless (scalar @data == $motif_length);
    #exit(1);
  } else {
  	# anything else is the end of useful data: don't bother to check
  	last;
  }
}


my $num_motifs = 0;

 ###### Decide whether to print the motif.

 # If no criteria are given, then print it.
 my $print_it = 1;
 
 $motif_name = "1";
   
 # Print the motif.
 if ($print_it) {
   $num_motifs++;
   print(STDERR "Printing motif $motif_name.\n");
   if ($use_numbers) {
     print("MOTIF $num_motifs $motif_name\n\n");
     print("BL   MOTIF $num_motifs width=$motif_length seqs=$num_seqs\n");
   } else { 
     print("MOTIF $motif_name\n\n");
     print("BL   MOTIF $motif_name width=$motif_length seqs=$num_seqs\n");
   }

   # PSSM for MAST
   if (0) {
     print("log-odds matrix: alength= $num_bases w= $motif_length n= 0 bayes= 0 E= 0\n");
     for (my $i = 0; $i < $motif_length; $i++) {
foreach my $base (@bases) {
  printf(
    "%7.3f ", 
    (log($motifs{$motif_name}->{freqs}->{$base}->[$i])/ $bg{$base}) / log(2.0)
  );
}
print("\n");
     }
   }
   # PSFM for Meta-MEME
   print("letter-probability matrix: ");
   print("alength= $num_bases w= $motif_length nsites= $num_seqs E= 0\n");
   my @temp = keys %motif_freqs;
   for (my $i = 0; $i < $motif_length; $i++) {
   	for my $freq (sort keys %motif_freqs) {
   		print "$motif_freqs{$freq}[$i] ";
   	}
   	print "\n";
   }
   print("\n");
 } else {
   $num_skipped++;
 }


print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");


close(MATRIX_FILE);
