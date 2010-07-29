#!@WHICHPERL@ -w
#
# FILE: uniprobe2meme
# AUTHOR: Timothy L. Bailey
# CREATE DATE: 18/12/2009
# DESCRIPTION: Convert a file containing list of TFBS motif matrices from 
# CHEN to MEME output format.

use warnings;
use strict;

# requires
push(@INC, split(":", $ENV{'PATH'}));   # look in entire path
require "motif2meme.pm";

# Set up global variables.
my $pseudo_total = 0; 		# default total pseudocounts
my $print_logodds = 0;
my $alph_type = "DNA";
my $alphabet = "ACGT";

my $usage = "USAGE: chen2meme [options]

  Options: 
	-skip <ID> 	   		skip this ID (may be repeated)
	-bg <background file>		set of f_a
	-pseudo <total pseudocounts>	add <total pseudocounts> times f_a to 
                                        each freq default: $pseudo_total\n
	-logodds			print log-odds matrix, too; 
					default: print frequency matrix only
	-h				print usage message

  Read a CHEN (concatenated) matrix file and convert to MEME format.

  Reads standard input.
  Writes standard output\n\n";

my $next_arg;
my $bg_file;
my %skips;
while (scalar(@ARGV) > 0) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-h") {
    printf(STDERR $usage);
    exit(0);
  } elsif ($next_arg eq "-skip") {
    $skips{shift(@ARGV)} = 1;
  } elsif ($next_arg eq "-bg") {
    $bg_file = shift(@ARGV);
  } elsif ($next_arg eq "-pseudo") {
    $pseudo_total = shift(@ARGV);
  } elsif ($next_arg eq "-logodds") {
    $print_logodds = 1;
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    printf(STDERR $usage);
    exit(1);
  }
}

my @residues = split //, $alphabet;
my $output = "";

# Get the background model.
my %bg = &read_background_file($alphabet, $bg_file);

# Print the MEME header.
$output = &print_meme_header($alph_type, $alphabet, \@residues, \%bg);

# Read the input file.
my $num_skipped = 0;
my (%motifs, $motif_name, $motif_freqs, $col);
while (<>) {

  next if (/^\s*$/); 		# skip blank lines

  if (/^>(\S+)/) {
    # motif id line
    my $id = $1;
    $motif_name = $id;
    # create the motif entry
    $motifs{$motif_name} = {};
    $motifs{$motif_name}->{name} = $motif_name;
    $motifs{$motif_name}->{descr} = "";
    $motifs{$motif_name}->{freqs} = {};
    $motif_freqs = $motifs{$motif_name}->{freqs};
    $col = 0;
  } else {
    # motif "column" line 
    my @col_freqs = split;
    my $i = 0;
    # map the line to the current column of the matrix
    foreach my $base (@residues) {
      $motif_freqs->{$base}[$col] = $col_freqs[$i++]; 
    }
    $col++;		# next column
  }

}

my $num_motifs = 0;
for my $motif_name (sort keys %motifs) {

  # Print the motif.
  if ( defined($skips{$motif_name})) {
    $num_skipped++;
  } else {
    $motifs{$motif_name}->{width} = scalar @{$motifs{$motif_name}->{freqs}->{$residues[0]}};
    my $motif_href = $motifs{$motif_name};
    &convert_counts_to_freqs($motif_href, \%bg, \@residues, $pseudo_total);
    $output .= &print_meme_motif($motif_href, \%bg, \@residues, $print_logodds);
    $num_motifs++;
  }

}
print(STDERR "Converted $num_motifs motifs.\n");
print(STDERR "Skipped $num_skipped motifs.\n");

print $output unless ($num_motifs == 0);
