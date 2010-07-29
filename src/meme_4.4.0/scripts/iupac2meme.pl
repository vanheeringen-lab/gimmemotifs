#!@WHICHPERL@ -w
#
# $Id:$
#
# FILE: iupac2meme
# AUTHOR: Timothy L. Bailey
# CREATE DATE: 3/06/2010
# DESCRIPTION: Convert a DNA IUPAC motif to MEME format

use warnings;
use strict;

# requires
push(@INC, split(":", $ENV{'PATH'}));   # look in entire path
require "motif2meme.pm";

# Set up global variables.
my $alph_type = 'DNA';
my $pseudo_total = 0; 		# default total pseudocounts
my $num_seqs = 20;
my $print_logodds = 0;
my $dna_alphabet = "ACGT";
my %dna_codes = (
  'A' => ['A'],
  'C' => ['C'],
  'G' => ['G'],
  'T' => ['T'],
  'U' => ['T'],
  'R' => ['A','G'],
  'Y' => ['C','T'],
  'M' => ['C','A'],
  'K' => ['T','G'],
  'W' => ['T','A'],
  'S' => ['C','G'],
  'B' => ['C','T','G'],
  'D' => ['A','T','G'],
  'H' => ['A','T','C'],
  'V' => ['A','C','G'],
  'N' => ['A','C','G','T']
);
my $protein_alphabet = "ACDEFGHIKLMNPQRSTVWY";
my %protein_codes = (
  'A' => ['A'],
  'C' => ['C'],
  'D' => ['D'],
  'E' => ['E'],
  'F' => ['F'],
  'G' => ['G'],
  'H' => ['H'],
  'I' => ['I'],
  'K' => ['K'],
  'L' => ['L'],
  'M' => ['M'],
  'N' => ['N'],
  'P' => ['P'],
  'Q' => ['Q'],
  'R' => ['R'],
  'S' => ['S'],
  'T' => ['T'],
  'V' => ['V'],
  'W' => ['W'],
  'Y' => ['Y'],
  'B' => ['D','N'],
  'Z' => ['E','Q'],
  'X' => ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
);

my $usage = "USAGE: iupac2meme [options] [<iupac_motif>]+

  Options: 
	-alpha dna|protein		IUPAC alphabet; default: $alph_type
	-numseqs <numseqs>		assume frequencies based on this many
					sequences; default: $num_seqs
	-bg <background file>		file with background frequencies of 
					letters, {f_a}; default: uniform background
	-pseudo <total pseudocounts>	add <total pseudocounts> times f_a to 
                                        each frequency; default: $pseudo_total
	-logodds			output the log-odds (PSSM) and frequency 
					(PSPM) motifs; default: PSPM motif only

  Convert an IUPAC motif to MEME format.

  Example IUPAC DNA motif: ACGGWNNYCGT
  Example IUPAC PROTEN motif: IKLVBZYXXHG\n\n";

# Process command line arguments.
if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}

my $next_arg;
my $bg_file;
my @motif_strings;
while (scalar(@ARGV) > 0) {
  $next_arg = shift(@ARGV);
  if ($next_arg eq "-numseqs") {
    $num_seqs = shift(@ARGV);
  } elsif ($next_arg eq "-alpha") {
    $alph_type = shift(@ARGV);
    $alph_type =~ y/a-z/A-Z/;
  } elsif ($next_arg eq "-bg") {
    $bg_file = shift(@ARGV);
  } elsif ($next_arg eq "-pseudo") {
    $pseudo_total = shift(@ARGV);
  } elsif ($next_arg eq "-logodds") {
    $print_logodds = 1;
  } elsif ($next_arg !~ /^-/){
    push(@motif_strings, $next_arg);
  } else {
    print(STDERR "Illegal argument ($next_arg)\n");
    exit(1);
  }
}

my $iupac_motif = $motif_strings[0];
unless (scalar @motif_strings > 0) {
  print(STDERR "You must provide a motif.\n");
  exit(1);
}

# set alphabet and codes
my ($alphabet, %codes);
if ($alph_type eq 'DNA') {
  $alphabet = $dna_alphabet;
  %codes = %dna_codes;
} elsif ($alph_type eq 'PROTEIN') {
  $alphabet = $protein_alphabet;
  %codes = %protein_codes;
} else {
  die ("Unknown alphabet type: $alph_type\n");
}
my @residues = split(//, $alphabet);
my $num_residues = scalar @residues;

# get the background model
my %bg = &read_background_file($alphabet, $bg_file);

my (@letters, $width, $motif_name, %motifs, %motif_freqs, $motif_freqs);

#
# convert the IUPAC motifs to frequencies.
#
foreach $iupac_motif (@motif_strings) {
  # set up the dictionary of motifs
  $motif_name = $iupac_motif;
  $iupac_motif =~ y/a-z/A-Z/;		# convert motif to upper case
  @letters = split(//, $iupac_motif);
  $width = scalar @letters;
  $motifs{$motif_name} = {};
  $motifs{$motif_name}->{name} = $motif_name;
  $motifs{$motif_name}->{descr} = "";
  $motifs{$motif_name}->{width} = $width;
  $motifs{$motif_name}->{num_seqs} = $num_seqs;
  $motifs{$motif_name}->{freqs} = {};
  $motif_freqs = $motifs{$motif_name}->{freqs};

  # initialize the counts with the pseudocounts
  for (my $i=0; $i<$width; $i++) {
    foreach my $residue (@residues) {
      $motif_freqs->{$residue}->[$i] = $bg{$residue}*$pseudo_total;
    }
  }

  # fill in the frequencies
  for (my $i=0; $i<$width; $i++) {
    my $c = $letters[$i];
    my $equivs = $codes{$c};
    die ("Unknown letter '$c' in $alph_type IUPAC motif.\n") unless defined $equivs;
    my $num_equivs = scalar @$equivs;
    foreach my $residue (@$equivs) {
      $motif_freqs->{$residue}->[$i] += $num_seqs/$num_equivs;
    }
    foreach my $residue (@residues) {
      $motif_freqs->{$residue}->[$i] /= $num_seqs + $pseudo_total;
    }
  }
} # convert motif to frequencies

# Print the MEME header
print(&print_meme_header($alph_type, $alphabet, \@residues, \%bg));

#
# output each motif in MEME format
#
my $num_motifs = 0;
for my $motif_name (sort keys %motifs) {
  $num_motifs++;
  my $motif_href = $motifs{$motif_name};
  print(&print_meme_motif($motif_href, \%bg, \@residues, $print_logodds));
}

#print(STDERR "Converted $num_motifs motifs.\n");
