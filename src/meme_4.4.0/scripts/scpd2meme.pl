#!@WHICHPERL@ -w
#
# $Id $
# $Log $
#
# FILE: scpd2meme
# AUTHOR: William Stafford Noble, Timothy L. Bailey, Charles Grant
# CREATE DATE: 1/06/2006
# DESCRIPTION: Convert a file containing list of TFBS motif matrices from SCPD
# to MEME output format.

use warnings;
use strict;

# Set up global variables. Assume uniform.
my $pseudo_total = 0; 		# default total pseudocounts
my @bases  = ("A", "C", "G", "T");
my $num_bases  = 4;
my %bg = ( 'A' => .25, 'C' => .25, 'G' => .25, 'T' => .25,);
my $STRANDS = 2;
my $use_numbers = 0;

my $usage = "USAGE: scpd2meme [options] <matrix file>

  Options: 
	-skip <scpd ID> 		skip this ID (may be repeated)
	-numbers       		  	use numbers instead of strings as motif IDs
	-bg <background file>		set of f_a
	-pseudo <total pseudocounts>	add <total pseudocounts> times f_a to 
                                        each freq default: $pseudo_total\n\n";

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
  if ($next_arg eq "-skip") {
    $skips{shift(@ARGV)} = 1;
  } elsif ($next_arg eq "-numbers") {
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

# read the background file
if (defined($bg_file)) {
  my $a;
  my $f;
  my $total_bg;
  open(BG_FILE, $bg_file) || die("Can't open $bg_file.\n");
  $total_bg = 0;
  while (<BG_FILE>) {
    next if (/^#/);      # skip comments
    ($a, $f) = split;
    if ($a eq "A" || $a eq "a") {
      $bg{"A"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "C" || $a eq "c") {
      $bg{"C"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "G" || $a eq "g") {
      $bg{"G"} = $f; 
      $total_bg += $f;
    } elsif ($a eq "T" || $a eq "t") {
      $bg{"T"} = $f; 
      $total_bg += $f;
    }
  }
  # make sure they sum to 1
  my $key;
  foreach $key (keys %bg) {
    $bg{$key} /= $total_bg;
    #printf STDERR "$key $bg{$key}\n";
  }
  close BG_FILE;
}  # background file

# Open the matrix file for reading.
open(MATRIX_FILE, "<$matrix_file") || die("Can't open $matrix_file.\n");

# Print the MEME header.
print("MEME version 3.0\n\n");
print("ALPHABET= ACGT\n\n");
print("strands: + -\n\n");
print(
  "Background letter frequencies " .
  "(from dataset with add-one prior applied):\n"
);
printf("A %f C %f G %f T %f\n\n",  $bg{"A"}, $bg{"C"}, $bg{"G"}, $bg{"T"});

# Read the input file.
my $num_skipped = 0;
my $line_number = 0;
my %motifs;
my $line;
while ($line = <MATRIX_FILE>) {

  $line_number++;

  # Split the line into identifier and everything else.
  (my $id, my @data) = split('\t', $line);

  # An ID should start with '>'
  if (not defined($id) || (substr($id, 0, 3) ne ">")) {
    die ("No matrix id found at line $line_number in file $matrix_file.\n");
  }

  my $motif_name = substr($id,1);
  my $width;
  my @totals;
  $motifs{$motif_name} = {};
  $motifs{$motif_name}->{counts} = {};
  $motifs{$motif_name}->{freqs} = {};
  my $motif_counts = $motifs{$motif_name}->{counts};
  my $motif_freqs = $motifs{$motif_name}->{freqs};

  # Read the motif base counts. One line for each nucleotide.
  for(my $base_index = 0; $base_index < $num_bases; $base_index++) {
    
    $line = <MATRIX_FILE>;
    chomp($line);

    if (!defined $line ) {
      die ("Missing motif counts at line $line_number in file $matrix_file.\n");
    }
    $line_number++;

    (my $base, my @counts) = split('\s+', $line);
    $motif_counts->{$base} = \@counts; 
    my $width = scalar @counts;
  
    if ($base_index == 0) {
      # Initialize totals for each motif position.
      for (my $i = 0; $i < $width; $i++) {
        $totals[$i] = 0;
      }
     $motifs{$motif_name}->{width} = $width;
    }
    if ($width != $motifs{$motif_name}->{width}) {
      die(
        "Inconsistent number of positions in motif at line $line_number "
        . " in file $matrix_file"
      );
    }
    # Track total counts at each motif position.
    for (my $i = 0; $i < $width; $i++) {
      $totals[$i] += $motif_counts->{$base}->[$i];
    }
  }
  
  $motifs{$motif_name}->{totals} = \@totals;

  # Convert the base counts at each motif position into frequencies.
  for (my $i = 0; $i < $motifs{$motif_name}->{width}; $i++) {
    for (my $base_index = 0; $base_index < $num_bases; $base_index++) {
      my $base = $bases[$base_index];
      $motif_freqs->{$base}->[$i] = ($motif_counts->{$base}->[$i] 
        + ($pseudo_total * $bg{$base})) / ($totals[$i] + $pseudo_total);
    }
  }
}

my $num_motifs = 0;
for my $motif_name (keys %motifs) {

  ###### Decide whether to print the motif.

  # If no criteria are given, then print it.
  my $print_it = 1;

  # Were we explicitly asked to skip this one?
  if (defined($skips{$motif_name})) {
    $print_it = 0;
  } 
    
  # Print the motif.
  if ($print_it) {
    $num_motifs++;
    my $width = $motifs{$motif_name}->{width};
    my $num_seqs = $motifs{$motif_name}->{totals}->[0];
    print(STDERR "Printing motif $motif_name.\n");
    if ($use_numbers) {
      print("MOTIF $num_motifs $motif_name\n\n");
      print("BL   MOTIF $num_motifs width=$width seqs=$num_seqs\n");
    } else { 
      print("MOTIF $motif_name\n\n");
      print("BL   MOTIF $motif_name width=$width seqs=$num_seqs\n");
    }

    # PSSM for MAST
    if (0) {
      print("log-odds matrix: alength= $num_bases w= $width n= 0 bayes= 0 E= 0\n");
      for (my $i = 0; $i < $width; $i++) {
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
    print("alength= $num_bases w= $width nsites= $num_seqs E= 0\n");
    for (my $i = 0; $i < $width; $i++) {
      foreach my $base (@bases) {
        printf("%7.3f ", $motifs{$motif_name}->{freqs}->{$base}->[$i]);
      }
      print("\n");
    }
    print("\n");
  } else {
    $num_skipped++;
  }

}
print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");


close(MATRIX_FILE);
