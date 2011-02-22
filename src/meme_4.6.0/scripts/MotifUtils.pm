#File: MotifUtils.pm
#Project: generic
#Author: James Johnson (though used parts of exising scripts)
#Created: June 2010 (though the repackaged scripts are much older)

package MotifUtils;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(seq_to_intern matrix_to_intern intern_to_iupac intern_to_meme read_background_file meme_header motif_id motif_width parse_double round);

use POSIX qw(strtod);

##############################################################################
# Internal format description
##############################################################################
#
# The internal format for the motif is as follows:
#
# hash {
#   bg => hash {
#     <residue 1> => scalar (background prior probability of residue 1)
#     <residue 2> => scalar (background prior probability of residue 2)
#      ...
#     <residue N> => scalar (background prior probability of residue N)
#     dna => scalar (true if the bg is dna)
#     source => scalar (source description of bg)
#   }
#   strands => scalar (number of strands used in a dna motif)
#   id => scalar (name of motif)
#   alt => scalar (alternate name of motif)
#   url => scalar (url of motif)
#   width => scalar (width of the motif)
#   sites => scalar (number of sites that made up motif)
#   pseudo => scalar (total pseudocount added)
#   evalue => scalar (evalue of motif)
#   pspm => hash {
#     <residue 1> => array [
#       scalar (probability of residue 1 at position 1)
#       scalar (probability of residue 1 at position 2)
#       ...
#       scalar (probability of residue 1 at position M)
#     ]
#     <residue 2> => array [
#       scalar (probability of residue 2 at position 1)
#       scalar (probability of residue 2 at position 2)
#       ...
#       scalar (probability of residue 2 at position M)
#     ]
#   }
#   pssm (optional for meme files that have a better pssm) => hash {
#     <residue 1> => array [
#       scalar (probability of residue 1 at position 1)
#       scalar (probability of residue 1 at position 2)
#       ...
#       scalar (probability of residue 1 at position M)
#     ]
#     <residue 2> => array [
#       scalar (probability of residue 2 at position 1)
#       scalar (probability of residue 2 at position 2)
#       ...
#       scalar (probability of residue 2 at position M)
#     ]
#   }
# }
#
#
#
#


##############################################################################
# lookup tables
##############################################################################
my @dna_alphabet = qw(A C G T);
my @dna_iupac = qw(A C M G R S V T W Y H K D B N); #ordered by bit value
my %dna_bits = ('A' => 1, 'C' => 2, 'G' => 4, 'T' => 8, 'U' => 8, 'R' => 5, 
  'Y' => 10, 'M' => 3, 'K' => 12, 'W' => 9, 'S' => 6, 'B' => 14, 'D' => 13, 
  'H' => 11, 'V' => 7, 'N' => 15);
my $dna_ic = 2;

my @protein_alphabet = qw(A C D E F G H I K L M N P Q R S T V W Y);
my @protein_iupac = qw(A C D E F G H I K L M N P Q R S T V W Y B Z X);
my %protein_bits = ('A' => 1, 'C' => 2, 'D' => 4, 'E' => 8, 'F' => 16, 
  'G' => 32, 'H' => 64, 'I' => 128, 'K' => 256, 'L' => 512, 'M' => 1024, 
  'N' => 2048, 'P' => 4096, 'Q' => 8192, 'R' => 16384, 'S' => 32768, 
  'T' => 65536, 'V' => 131072, 'W' => 262144, 'Y' => 524288, 
  'B' => 2052, # D or N
  'Z' => 8200, # E or Q
  'X' => 1048575 # any
);
my %protein_lookup = (1 => 'A', 2 => 'C', 4 => 'D', 8 => 'E', 16 => 'F', 
  32 => 'G', 64 => 'H', 128 => 'I', 256 => 'K', 512 => 'L', 1024 => 'M', 
  2048 => 'N', 4097 => 'P', 8192 => 'Q', 16284 => 'R', 32768 => 'S', 
  65536 => 'T', 131072 => 'V', 262144 => 'W', 524288 => 'Y', 2052 => 'B',
  8200 => 'Z', 1048575 => 'X');

my $protein_ic = log(20) / log(2);

my $LN2 = log(2);

##############################################################################
# public functions
##############################################################################

#
# Takes one or more iupac sequences with posssible bracket expressions 
# and constructs a probability matrix in this modules internal format. 
#
# Usage:
#  seq_to_intern(
#   \%background,
#   $sequence_samples,
#   $sites_per_sample,
#   $pseudo_total,
#   id => $optional_id,
#   alt => $optional_alt,
#   strands => $optional_strands,
#   url => $optional_url
#   );
#
# Parameters
#  bg                 - a reference to the background model also defines if the motif is dna or protein
#  sequence samples   - one or more iupac sequences seperated by newlines.
#  sites per sample   - the number of sites attributed to each sample;
#  pseudo total       - a pseudocount that is distributed according to the background model
#
# Options
#  id                 - use value for motif name
#  alt                - use value for motif alternate name 
#  strands            - only relevant to dna; the strands used to create the motif
#  url                - use value for motif information website
#
sub seq_to_intern {
  my ($bg, $seq_lines, $sites_per_sample, $pseudo_total, %opt) = @_;
  my ($id, $alt, $strands, $url) = ($opt{id}, $opt{alt}, $opt{strands}, $opt{url});
  $strands = ($bg->{dna} ? 2 : 1) unless defined $strands;
  $alt = '' unless defined($alt);
  # split the sequence lines into many sequences
  my @seqs = split(/\n/, $seq_lines);

  #lookups
  my %nums = ($bg->{dna} ? %dna_bits : %protein_bits);
  my @residues = ($bg->{dna} ? @dna_alphabet : @protein_alphabet);
  my $alpha = join('', keys(%nums));

  #keep track of errors related to parsing seq
  my @errors = ();

  #keep track of counts from multple samples
  my @matrix = ();
  my $name = '';

  # get the dimensions
  my $asize = scalar(@residues);
  my $width = 0;

  my $sample_count = 0;
  for (my $line = 0; $line < scalar(@seqs); $line++) {
    my ($bitfields, $convert_errors) = seq_to_bitfield($seqs[$line], $alpha, \%nums, $line + 1);
    next if (scalar(@{$bitfields}) == 0 && !@{$convert_errors}); #skip empty lines
    push(@errors, @{$convert_errors});
    if ($width == 0) {
      $width = scalar(@{$bitfields});
      if ($width > 0) { #usable line found
        #initilize matrix to all zeros
        my $matrix_size = $asize * $width;
        for (my $i = 0; $i < $matrix_size; $i++) {
          $matrix[$i] = 0;
        }
        #create a name (more accurate than intern_to_iupac for protein)
        $name = bitfield_to_seq($bg->{dna}, $bitfields) unless $id;
      }
    }
    if ($width == 0 || scalar(@{$bitfields})  != $width) { #skip bad lines
      push(@errors, errmsg($line + 1, $seqs[$line], "Bad width. Sample skipped.", 0));
      next;
    }
    for (my $pos = 0; $pos < $width; $pos++) {
      my $bits = $bitfields->[$pos];
      my $fract = $sites_per_sample / bits_set($bits);
      for (my $a = 0; $a < $asize; $a++) {
        $matrix[$pos * $asize + $a] += $fract if ($bits & (1 << $a)); # apply seq
      }
    }
    $sample_count++;
  } # parse line

  if ($width == 0) {
    push(@errors, "Motif has no width.\n");
    return (undef, \@errors);
  }

  my $sites = $sample_count * $sites_per_sample;
  #initilize the internal motif datastructure
  my %motif = (id => $id , alt => $alt, url => $url, strands => $strands, width => $width, 
      sites => $sites, pseudo => $pseudo_total, bg => $bg, evalue => 0, pspm => {});
  my $motif_pspm = $motif{pspm};

  my $total = $sites + $pseudo_total;
  for (my $a = 0; $a < $asize; $a++) {
    my $residue = $residues[$a];
    my $pseudo_count = $bg->{$residue} * $pseudo_total;
    $motif_pspm->{$residue} = [];
    for (my $pos = 0; $pos < $width; $pos++) {
      $motif_pspm->{$residue}->[$pos] = ($pseudo_count + $matrix[$pos * $asize + $a]) / $total;
    }
  }

  unless ($id) {
    if ($sample_count > 1) { # if multiple samples were used then convert to iupac
      $name = intern_to_iupac(\%motif);
    }
    $motif{id} = $name;
  }

  return (\%motif, \@errors);
}


#
# Converts a motif in matrix form into the internal format.
#
# Usage:
#  matrix_to_intern(
#   \%bg,
#   $matrix,
#   $orientation,
#   $site_count,
#   $pseudo_total,
#   alt => $alt
#   );
#
# Parameters
#  bg            - a reference to the background model also defines if motif is dna or protein
#  matrix        - a matrix; columns are space separated and rows are newline separated
#  orientation   - the matrix orientation; allowed values are 'auto', 'col' and 'row'
#  site count    - the number of sites in the matrix; ignored for count matrices
#  pseudo total  - a pseudocount that is distributed according to the background model
#
# Options
#  id            - use value for motif name
#  alt           - use value for motif alternate name 
#  strands       - only relevant to dna; the strands used to create the motif
#  url           - use value for motif information website
#  rescale       - force the use of the site count, even when the matrix is already counts
#
sub matrix_to_intern {
  my ($bg, $matrix, $orientation, $site_count, $pseudo_total, %opt) = @_;
  my ($id, $alt, $strands, $url, $rescale) = ($opt{id}, $opt{alt}, $opt{strands}, $opt{url}, $opt{rescale});
  $strands = ($bg->{dna} ? 2 : 1) unless defined $strands;
  $alt = '' unless defined($alt);
  $rescale = 0 unless defined($rescale);

  #initilize the internal motif datastructure
  #NB width refers to the motif width not the matrix width
  my %motif = (id => $id, alt => $alt, url => $url, strands => $strands, width => 0, 
      sites => 0, pseudo => $pseudo_total, bg => $bg, evalue => 0, pspm => {});
  my $motif_pspm = $motif{pspm};
  
  # lookups
  my @residues = ($bg->{dna} ? @dna_alphabet : @protein_alphabet);
  my $asize = scalar(@residues);

  # dimensions (set defaults)
  my $height = 0; 
  my $width = ($orientation eq 'row' ? $asize : -1); 

  # parse matrix into an array
  my $total = 0;
  my @matrix_array = ();
  my @errors = ();
  my @lines = split(/\n/, $matrix);
  for (my $line_i = 0; $line_i < scalar(@lines); $line_i++) {
    my $line = $lines[$line_i];
    $line =~ s/^\s+//;#trim left
    $line =~ s/\s+$//;#trim right
    # skip empty lines
    next if $line eq '';
    my @nums = split(/\s+/, $line);
    if ($width == -1) {
      $width = scalar(@nums); #expected width
    } elsif (scalar(@nums) != $width) {
      push(@errors, "Error: Expected $width elements on line $line_i but got " . scalar(@nums) . ".");
      return (undef, \@errors);
    }
    # count this row
    $height++;
    for (my $num_i = 0; $num_i < scalar(@nums); $num_i++) {
      my $num = parse_double($nums[$num_i]);
      if (not defined $num) {
        push(@errors, "Warning: Element $num_i is not a number.");
        $num = 0; #specify default
      }
      push(@matrix_array, $num);
      $total += $num;
    }
  }
  # check for an empty matrix
  if ($height == 0) {
    push(@errors, "Error: Empty matrix.");
    return (undef, \@errors);
  }
  # check if the dimensions are plausible
  if ($orientation eq 'col' && $height != $asize) {
    push(@errors, "Error: Expected $asize rows but got $height.");
    return (undef, \@errors);
  }
  if ($orientation eq 'auto' && $width != $asize && $height != $asize) {
    push(@errors, "Error: Expected either the row or column count to be the alphabet size $asize but got $width by $height.");
    return (undef, \@errors);
  }
  # now determine the orientation and transform the matrix to a row matrix
  if ($orientation eq 'auto') {
    if ($width == $asize && $height == $asize) {
      # need to guess orientation so use variance
      my $avg = $total / $asize;
      my $row_variance = 0;
      my $col_variance = 0;
      for (my $i = 0; $i < $asize; $i++) {
        my $row_sum = 0;
        my $col_sum = 0;
        for (my $j = 0; $j < $asize; $j++) {
          $row_sum += $matrix_array[$i * $width + $j];
          $col_sum += $matrix_array[$j * $width + $i];
        }
        $row_variance += ($row_sum - $avg) ** 2; # row sum to the power of 2
        $col_variance += ($col_sum - $avg) ** 2; # col sum to the power of 2;
      }
      $row_variance /= $asize;
      $col_variance /= $asize;
      if ($row_variance > $col_variance) { # probably a column matrix
        @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
      }
    } elsif ($height == $asize) { # dimensions suggest column matrix
      @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
      $height = $width;
      $width = $asize;
    }
  } elsif ($orientation eq 'col') {
    @matrix_array = transpose_matrix(\@matrix_array, $width, $height);
    $height = $width;
    $width = $asize;
  }
  #should now have a row matrix
  if ($rescale) {
    if (not defined($site_count)) {
      push(@errors, "Error: Rescale option requires a site count. Can't continue without a site count.");
      return (undef, \@errors);
    }
  } else {
    #detect probability matrix
    my $row_avg = $total / $height;
    if (int($row_avg + 0.5) > 1) { # probably a count matrix
      $site_count = int($row_avg + 0.5); #int truncates toward zero so this rounds for positive numbers
    } elsif (not defined($site_count)) {
      push(@errors, "Error: Expected count matrix but got probability matrix. Can't continue without site count.");
      return (undef, \@errors);
    }
  }

  #set sites and width
  $motif{sites} = $site_count;
  $motif{width} = $height;

  # initialize the counts with the pseudocounts
  foreach my $residue (@residues) {
    my $pseudocount = $bg->{$residue} * $pseudo_total;
    for (my $i = 0; $i < $width; $i++) {
      $motif_pspm->{$residue}->[$i] = $pseudocount;
    }
  }

  # normalise each row individually
  my $row_total = $site_count + $pseudo_total;
  for (my $y = 0; $y < $height; ++$y) {
    my $sum = 0; #calculate this row's sum
    for (my $x = 0; $x < $width; ++$x) {
      $sum += $matrix_array[$y*$width + $x];
    }
    if ($sum == 0) {
      push(@errors, "Error: Motif position $y summed to zero.");
      return (undef, \@errors);
    }
    # normalise
    for (my $x = 0; $x < $width; ++$x) { 
      $matrix_array[$y*$width + $x] /= $sum;
    }
    # copy to motif
    for (my $x = 0; $x < $width; ++$x) {
      my $residue = $residues[$x];
      # add counts
      $motif_pspm->{$residue}->[$y] += $matrix_array[$y*$width + $x] * $site_count;
      # convert to probabilities
      $motif_pspm->{$residue}->[$y] /= $row_total;
    }
  }

  unless (defined $id) {
    $id = intern_to_iupac(\%motif);
    $motif{id} = $id;
  }

  # now make motif
  return (\%motif, \@errors);
}

#
# intern_to_iupac(
#   $motif
#   );
#
# Measures the eucledian distance between each motif position and each iupac letter
# and selects the closest letter.
#
sub intern_to_iupac {
  my ($motif) = @_;
  my %bg = %{$motif->{bg}};
  my $sites = $motif->{sites};
  my $pseudo = $motif->{pseudo};
  my @residues = ($bg{dna} ? @dna_alphabet : @protein_alphabet);
  my @iupac = ($bg{dna} ? @dna_iupac : @protein_iupac);
  my %bits = ($bg{dna} ? %dna_bits : %protein_bits);
  my $alpha_ic = ($bg{dna} ? $dna_ic : $protein_ic);

  
  #calculate iupac code ic's assuming same sites and background as motif
  my @iupac_data = ();
  foreach my $code (@iupac) {
    my $set = $bits{$code};
    my $counts = $sites / &bits_set($set);
    my $ic = $alpha_ic;
    my %probs;
    my %fractions;
    for (my $r = 0; $r < scalar(@residues); $r++) {
      my $residue = $residues[$r];
      my $p = (($set & 1<<$r ? $counts : 0) + $bg{$residue} * $pseudo) / ($sites + $pseudo);
      $probs{$residue} = $p;
      $ic -= ($p ? ($p * (log($p) / $LN2)) : 0);
    }
    foreach my $residue (@residues) {
      $fractions{$residue} = $probs{$residue} * $ic;
    }
    push(@iupac_data, {code => $code, fractions => \%fractions});
  }


  my $out = '';
  for (my $i = 0; $i < $motif->{width}; $i++) {

    # calculate stack ic
    my $stack_ic = $alpha_ic;
    for (my $a = 0; $a < scalar(@residues); $a++) {
      my $residue = $residues[$a];
      my $p = $motif->{pspm}->{$residue}->[$i];
      next if $p == 0;
      $stack_ic -= ($p * (log($p) / $LN2));
    }

    my $stack_fraction = $stack_ic / $alpha_ic;
    my $best_ed = 100;
    my $best_code = '?';
    foreach my $data (@iupac_data) {
      my $fractions = $data->{fractions};
      my $ed = 0;
      foreach my $residue (@residues) {
        my $f1 = $motif->{pspm}->{$residue}->[$i] * $stack_fraction;
        my $f2 = $fractions->{$residue};
        $ed += ($f1 -$f2) ** 2
      }
      $ed = $ed ** 0.5;
      if ($ed < $best_ed) {
        $best_ed = $ed;
        $best_code = $data->{code};
      }
    }
    $out .= $best_code;
  }

  return $out;
}

#
# intern_to_meme(
#   $motif,
#   $add_pssm,
#   $add_pspm,
#   $add_header
#   );
#
#   Returns a motif in minimal meme format with the header optimal so multiple motifs can be concatenated.
#   Warning: MEME does special stuff with pseudo counts on proteins to generate the log odds matrix which 
#            this approach does not do!
#
sub intern_to_meme {
  my ($motif, $add_pssm, $add_pspm, $add_header) = @_;

  my %bg = %{$motif->{bg}};
  my @residues = ($bg{dna} ? @dna_alphabet : @protein_alphabet);

  my $output = "";

  #
  # get the text for the PSPM (and PSSM)
  #
  my $log_odds = "log-odds matrix: alength= ".scalar(@residues).
      " w= ".$motif->{width}." E= ".$motif->{evalue}."\n";
  my $letter_prob = "letter-probability matrix: alength= ".scalar(@residues).
      " w= ".$motif->{width}." nsites= ".$motif->{sites}." E= ".$motif->{evalue}."\n";
  for (my $i = 0; $i < $motif->{width}; $i++) {
    foreach my $residue (@residues) {
      #get the probability
      my $prob = $motif->{pspm}->{$residue}->[$i];
      #append to the letter probability matrix
      $letter_prob .= sprintf("%10.6f\t", $prob);
      my $score;
      if ($motif->{pssm}) {
        $score = $motif->{pssm}->{$residue}->[$i];
      } else {
        #calculate the log odds as log of zero is undefined give it a value of -10000 (TODO check does this make sense)
        $score = $prob > 0 ? round((log($prob/$bg{$residue}) / log(2.0)) * 100) : -10000;
      }
      #append to the log odds matrix
      $log_odds .= sprintf("%6d\t", $score);
    }
    $letter_prob .= "\n";
    $log_odds .= "\n";
  }
  $letter_prob .= "\n";
  $log_odds .= "\n";

  # Print the motif.
  $output .= meme_header($motif->{bg}, $motif->{strands}) if $add_header;
  $output .= "MOTIF ".$motif->{id}." ".$motif->{alt}."\n\n";
  $output .= $log_odds if ($add_pssm);
  $output .= $letter_prob if ($add_pspm);
  $output .= "URL ".$motif->{url}."\n\n" if ($motif->{url});

  return($output);
}

#
# read_background_file()
#
# Read a background file in fasta-get-markov format.
# If $bg_file is not defined, sets background to uniform.
# Contains the key/value 'source' that has the background 
# source description.
#
# Returns background as a hash.
#
sub read_background_file {
  my ($is_dna, $bg_file) = @_;

  # get the letters in the alphabet 
  my @letters = ($is_dna ? @dna_alphabet : @protein_alphabet);
  my $alph_length = scalar @letters;
  my (%bg, $a, $f);

  $bg{dna} = $is_dna;

  # initialize the background to uniform if no file given
  if (! defined $bg_file) {
    $bg{source} = "uniform background";
    foreach $a (@letters) {
      $bg{$a} = 1.0/$alph_length;
    }
    return(%bg);
  }
  
  # read the background file
  $bg{source} = "file `$bg_file'";
  open(BG_FILE, $bg_file) || die("Can't open $bg_file.\n");
  my $total_bg = 0;
  while (<BG_FILE>) {
    next if (/^#/);      		# skip comments
    ($a, $f) = split;
    next unless (length($a) == 1);	# skip higher order model
    $a =~ y/a-z/A-Z/;
    $bg{$a} = $f;
    $total_bg += $f;
  }
  close BG_FILE;

  # make sure they sum to 1 by normalizing
  foreach my $key (@letters) {
    die ("The letter '$key' was not given a value in the background file: $bg_file\n")
      unless (defined $bg{$key});
    $bg{$key} /= $total_bg;
  }

  return(%bg);
}  # background file

#
# meme_header(
#   \%background,
#   $strands
#   );
#
#   Returns the header part of minimal meme format.
#   Strands is only used if the motif is dna and defaults to 2.
#
#
sub meme_header {
  my ($bg, $strands) = @_;

  $strands = 2 unless defined $strands and $strands == 1;

  my @residues = ($bg->{dna} ? @dna_alphabet : @protein_alphabet);

  # Print the file MEME header.
  my $output;
  $output = "MEME version 4.4\n\n";
  $output .= "ALPHABET= ".join('',@residues)."\n\n";
  $output .= "strands: +\n\n" if ($bg->{dna} && $strands == 1);
  $output .= "strands: + -\n\n" if ($bg->{dna} && $strands == 2);
  $output .= "Background letter frequencies (from " . ($bg->{source} ? $bg->{source} : "unknown source") . "):\n";
  foreach my $residue (@residues) {
    $output .= sprintf("%s %.5f ", $residue, $bg->{$residue});
  }
  $output .= "\n\n";
  return($output);
}

#
# motif_id(
#   $motif
#   );
#
# Returns the motif id
#
sub motif_id {
  my $motif = shift;
  return $motif->{id};
}

#
# motif_width(
#   $motif
#   );
#
# Returns the width of the motif
#
sub motif_width {
  my $motif = shift;
  return $motif->{width};
}

#
# parse_double
#
# Parses a number using strtod
#
sub parse_double {
  my $str = shift;
  $str =~ s/^\s+//;#trim left
  $str =~ s/\s+$//;#trim right
  $! = 0; #clear errno
  my($num, $unparsed) = strtod($str);
  if (($str eq '') || ($unparsed != 0) || $!) {
      return undef;
  } else {
      return $num;
  } 
} 

sub round {
  # int truncates which is the same as always rounding toward zero
  if ($_[0] > 0) {
    return int($_[0] + 0.5);
  } else {
    return int($_[0] - 0.5);
  }
}

#
# transpose_matrix
#
# returns a transposed copy of the matrix.
#
sub transpose_matrix {
  my ($matrix_ref, $width, $height) = @_;
  my @matrix = @$matrix_ref;
  die("Dimensions don't match matrix size.") unless scalar(@matrix) == $width * $height;
  my @copy = @matrix;

  my $copy_width = $height;
  my $copy_height = $width;
  for (my $x = 0; $x < $width; ++$x) {
    for (my $y = 0; $y < $height; ++$y) {
      my $copy_x = $y;
      my $copy_y = $x;
      $copy[$copy_x + $copy_y * $copy_width] = $matrix[$x + $y * $width];
    }
  }
  return @copy;
}

#
# seq_to_bitfield
#
# returns an array of bitfields representing the sequence.
#
sub seq_to_bitfield {
  my ($seq, $alpha, $nums, $lineno) = @_;

  #error positions for making error msg
  my @unrecognised = ();
  my @unexpected_open = ();
  my @unexpected_close = ();
  my @nocontent = ();

  #parse the iupac codes and bracket expressions into an array of bitfields
  my @bitfields = ();

  my @chars = split(//, $seq);
  my $inside_bracket = 0;
  my $bracket_value;
  my $bracket_start;
  for (my $i = 0; $i < scalar(@chars); $i++) {
    my $char = $chars[$i];
    if (!$inside_bracket) { #not inside a bracket expression
      if ($char eq '[') { #found the start of a bracket expression
        $inside_bracket = 1;
        $bracket_value = 0;
        $bracket_start = $i;
      } elsif ($char =~ m/^([$alpha])$/i) { #found a normal character
        push(@bitfields, $nums->{uc($1)}); #add it to the bitfields
      } elsif ($char =~ m/^\s$/i) { #found whitespace
         #ignore
      } elsif ($char eq ']') { #misplaced bracket 
        push(@unexpected_close, $i); #add error
      } else { #unrecognised character
        push(@unrecognised, $i); #add error
      }
    } else { #inside a bracket expression
      if ($char eq ']') { #found the end of a bracket expression
        $inside_bracket = 0;
        if ($bracket_value > 0) { #bracket had content so add it to the bitfields
          push(@bitfields, $bracket_value);
        } else { #no valid characters within the bracket expression
          push(@nocontent, $bracket_start);
        }
      } elsif ($char =~ m/^([$alpha])$/i) { #found a normal character
        $bracket_value |= $nums->{uc($1)}; #add to the set of options that the bracket expression could equal
      } elsif ($char =~ m/^\s$/i) { #found whitespace
        #ignore
      } elsif ($char eq '[') { #misplaced bracket
        push(@unexpected_open, $i); #add error
      } else { #unrecognised character
        push(@unrecognised, $i); #add error
      }
    }
  }
  if ($inside_bracket) { #missing close bracket
    if ($bracket_value > 0) { #bracket had content so add it to the bitfields
      push(@bitfields, $bracket_value);
    } else { #no valid characters within the bracket expression
      push(@nocontent, $bracket_start); #add error
    }
  }
  my @errors = ();
  push(@errors, errmsg($lineno, $seq, "Unrecognised", @unrecognised));
  push(@errors, errmsg($lineno, $seq, "Unexpected ']'", @unexpected_close));
  push(@errors, errmsg($lineno, $seq, "Unexpected '['", @unexpected_open));
  push(@errors, errmsg($lineno, $seq, "No content in '[...]'", @nocontent));
  push(@errors, errmsg($lineno, $seq, "Missing ']'", $bracket_start)) if $inside_bracket;
  return (\@bitfields, \@errors);
}

#
# errmsg
#
# outputs two lines, the first showing the line where the problem occured,
# the second showing the positions of the error with ^ and the error message.
#
sub errmsg {
  my ($lineno, $line, $msg, @errors) = @_;
  return () unless @errors;
  my $pos = 0;
  my $indent = ' ' x length("$lineno: ");
  my $str = "$lineno: $line\n$indent";
  foreach my $i (@errors) {
    next if $i < $pos;
    $str .= ' ' x ($i - $pos);
    $str .= '^';
    $pos = $i + 1;
  }
  $str .= "  $msg";
  return ($str);
}

#
# bitfield_to_seq
#
# Converts an array of bitfields into a equalivent sequence making use of IUPAC codes and regular expression bracket expressions.
#
sub bitfield_to_seq {
  my ($is_dna, $bitfields) = @_;
  my $width = scalar(@{$bitfields});
  # generate a seq
  my $seq = "";
  if ($is_dna) {
    # it's easy for dna because there's an iupac code for every possible combination
    for (my $i = 0; $i < $width; $i++) {
      $seq .= $dna_iupac[$bitfields->[$i] - 1]; #arrays index from zero...
    }
  } else {
    my $asize = scalar(@protein_alphabet);
    # for protein there are only 3 special codes which do not cover every combination 
    for (my $i = 0; $i < $width; $i++) {
      my $bits = $bitfields->[$i];
      my $guess = $protein_lookup{$bits}; # attempt to find a single code
      if ($guess) {
        $seq .= $guess;
      } else { #some combination will be needed
        my $hasB = !(~$bits & $protein_bits{'B'});
        my $hasZ = !(~$bits & $protein_bits{'Z'});
        $bits &= ~$protein_bits{'B'} if $hasB; #turn off bits for B
        $bits &= ~$protein_bits{'Z'} if $hasZ; #turn off bits for Z
        $seq .= '[';
        for (my $j = 0; $j < $asize; $j++) {
          $seq .= $protein_alphabet[$j] if ($bits & (1 << $j));
        }
        $seq .= 'B' if $hasB;
        $seq .= 'Z' if $hasZ;
        $seq .= ']';
      }
    }
  }
  return $seq;
}

#
# bits
#
# returns the number of set bits
#
sub bits_set {
  my $num = shift;
  my $count = 0;
  for (; $num > 0; $num = $num>>1) {
    ++$count if ($num & 1);
  }
  return $count;
}
