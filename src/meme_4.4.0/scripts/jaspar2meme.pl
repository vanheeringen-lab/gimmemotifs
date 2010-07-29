#!@WHICHPERL@ -w

#**********************************************************************/
#* Copyright (c) University of Washington,                            */
#* Department of Genome Sciences, 2006. Written by Shobhit Gupta      */
#* and Timothy L. Bailey     					      */
#* All rights reserved.                                               */
#**********************************************************************/

# requires
push(@INC, split(":", $ENV{'PATH'}));   # look in entire path
require "motif2meme.pm";

#defaults:
my $alph_type = 'DNA';
my $pseudo_total = 0; 		# default total pseudocounts
my $alen = 4;
my $alphabet = "ACGT";
my $ext = "sites";
my $strands = 2;
my $use_numbers = 0;
my $descr = "";
my $print_logodds = 0;

my @letters = ('A', 'C', 'G', 'T');

my $usage = "USAGE: jaspar2meme [options] <Jaspar directory>

  Options: 
     -pfm	      read JASPAR count files (.pfm); 
		      default: site files (.sites)
     -cm	      read count file with line labels 'A|' etc. (.cm); 
		      default: site files (.sites)
     -numbers         use numbers instead of strings as motif IDs
     -strands [1|2]   print '+ -' '+' on the MEME strand line;
                      default: 2 (prints '+ -')
     -bg <bfile>      file with background frequencies in MEME
                      -bfile format; default: uniform frequencies
     -pseudo <A>      add <A> times background frequency to
                      each count when computing letter frequencies
                      default: $pseudo_total
     -descr <descr>   description of motif(s); printed if -numbers given
     -logodds         print log-odds matrix as well as frequency matrix;
                      default: frequency matrix only


  Convert a directory of JASPAR files into a MEME version 3 formatted 
  file suitable for use with MAST and other MEME Suite programs.

  A JASPAR '.sites' file describes a motif in terms of a multiple
  alignment of sites.  It contains a multiple alignment in modified 
  FASTA format.  Only capitalized sequence letters are part of the alignment.

  A JASPAR count file ('.pfm') contains a count matrix where the rows
  correspond to A, C, G and T, respectively.  

  A CM count file ('.cm') prefixes the rows with 'A| ', 'C| ', 'G| ' and 'T| '.

  A log-odds matrix and a probability matrix are output for each
  motif ('.sites') file.  The probability matrix is computed using
  pseudo-counts consisting of the background frequency (see -bg, above)
  multiplied by the total pseudocounts (see -pseudo, above).
  The log-odds matrix uses the background frequencies in the denominator
  and is log base 2.

  If a matrix_list.txt file exists and -pfm is given, the JASPAR names of the
  motifs are read from that file and included in the output.
  
  Writes standard output.
\n";

if (scalar(@ARGV) == 0) {
  printf(STDERR $usage);
  exit(1);
}

while($#ARGV > 0) {
  $arg = shift;
  if ($arg eq "-pfm") {
    $ext = "pfm";
  } elsif ($arg eq "-cm") {
    $ext = "cm";
  } elsif ($arg eq "-numbers") {
    $use_numbers = 1;
  } elsif ($arg eq "-strands") {
    $strands = shift @ARGV;
  } elsif ($arg eq "-bg") {
    $bg_file = shift @ARGV;
  } elsif ($arg eq "-pseudo") {
    $pseudo_total = shift @ARGV;
  } elsif ($arg eq "-logodds") {
    $print_logodds = 1;
  } elsif ($arg eq "-descr") {
    $descr = shift @ARGV;
  } else {
    print(STDERR "Illegal argument ($arg)\n");
    exit(1);
  }
}

while (scalar(@ARGV) > 0) {
  $JAS_DIR = shift(@ARGV);
}
$JAS_DIR =~ s/\/$//;

## Print usage if no argument
if ($JAS_DIR eq '') {
  print STDERR "\nJaspar directory missing\n";
  print STDERR "\nUsage:\n\n";
  print STDERR "jaspar2meme <Jaspar directory containing .sites files>\n\n";
  print STDERR "Meme version 3 format output file is written to stdout\n\n";
  exit (1);
}

if ($strands == 2) {
  $strand_string = " \+ \-"
} elsif ($strands == 1){
  $strand_string = " \+"
} else {
  print STDERR "value of \-strands not valid: should be 1 or 2 default is 2\n";
  exit (2);
}

#
# get the directory listing
#
my @files = `ls -1 $JAS_DIR`;

# get the background model
my %bg = &read_background_file($alphabet, $bg_file);

# get the matrix_list.txt file if there is one and
# get the motif names from it.
my %motif_names;
if ($ext eq "pfm") {
  foreach my $file (@files) {
    chop $file;
    #printf(STDERR "$file\n");
    if ($file eq "matrix_list.txt") {
      printf(STDERR "Found matrix_list.txt file.  Reading in motif names.\n");
      open(IN, "<$JAS_DIR/$file") || die "Can't open file $JAS_DIR/$file\n";
      while (<IN>) {
        my @words = split;
        #printf(STDERR "id %s name %s\n", $words[0], $words[2]);
        $motif_names{$words[0]} = $words[2];
      }
      close IN;
      last;
    }
  }
}

# 
# set up file header
#
my $file_header = &print_meme_header($alph_type, $alphabet, \@letters, \%bg);

my %motif = ();

#
# read in motif files and print motifs
#
my $id = 0;				# motif ID number
foreach my $file (@files) {

  next unless ($file =~ m/(\S*)\.$ext/);

  my $jaspar_id = $1;
  chomp $file;
  open(IN, "<$JAS_DIR/$file") || die "Can't open file $JAS_DIR/$file\n";

  my $motif_length;
  my %counts = ();
  my $n_sites = 0;
  my $header = "";

  if ($ext eq "pfm" || $ext eq "cm") {	# read a counts file (.pfm, .cm)
    my $letter_index = -1;
    while (<IN>) {
      next if (/^#/ || /^\s*$/); 		# skip comment, blank lines
      my @words = split;
      shift @words if ($ext eq "cm");	# discard line labels
      $motif_length = @words;
      $letter_index++;
      my $curr_letter = $letters[$letter_index];
      for (my $i=0; $i<$motif_length; $i++) {
        $counts{$i}{$curr_letter} = $words[$i]; 
        $n_sites += $words[$i];
      }
    }
    $n_sites /= $motif_length; 
  } elsif ($ext eq "sites") { 		# read a sites file
    my @seq = ();
    while (my $l = <IN>) {
      chomp $l;
      $l =~ s/\r//g;
      if ($l =~ m/^\>/) {
        $header = $l;
        $seq[$n_sites] = <IN>;
        chomp $seq[$n_sites];
        $seq[$n_sites] =~s /\r//g;
        $seq[$n_sites] =~ s/[a-z]//g;
        $n_sites ++;
      }
    } # End of file

    # convert sites to counts
    $motif_length = length $seq[0];
    for (my $query_child_index = 0; 
      $query_child_index < $n_sites;
      $query_child_index++
    ) {
      my $length_index;
      for ($length_index = 0; 
        $length_index < $motif_length;
        $length_index++
      ) {
        my $curr_letter = substr($seq[$query_child_index],
          $length_index,1);
        $counts{$length_index}{$curr_letter}++;
      }
    }
  } # read a file 

  #
  # print the MEME formatted output
  #
  print $file_header if ($id++ == 0);		# bump ID number

  $motif{width} = $motif_length;
  $motif{num_seqs} = $n_sites;
  $motif{freqs} = {};

  if (defined $motif_names{$jaspar_id}) {
    $header = $motif_names{$jaspar_id};
  } else {
    if (defined $motif_names{$jaspar_id}) {
      $header = $motif_names{$jaspar_id};
    } else {
      $header =~ s/(\d+)$//;
      $header =~ s/^([^\s]*\s)//;
    }
  }
  if ($use_numbers) {
    $motif{name} = $id;
    $motif{descr} = "$header\($jaspar_id\)";
  } else {
    $motif{name} = $jaspar_id;
    $motif{descr} = "$header $descr";
  }

  #
  # convert counts to frequencies
  #
  foreach (my $key1=0; $key1<$motif_length; $key1++) {
    my %vals = ();
    my $sum_row = 0;
    foreach my $key2 (@letters) {
      if (exists($counts{$key1}{$key2})){
        $sum_row += $vals{$key2} = $counts{$key1}{$key2} + ($pseudo_total * $bg{$key2});
      } else {
        $sum_row += $vals{$key2} = $pseudo_total * $bg{$key2};
      }
    }

    foreach my $key (@letters) {
      $motif{freqs}{$key}->[$key1] = $vals{$key}/$sum_row;
    }
  }

  #
  # print motif
  #
  print(&print_meme_motif(\%motif, \%bg, \@letters, $print_logodds));

  close IN;
} # foreach file
