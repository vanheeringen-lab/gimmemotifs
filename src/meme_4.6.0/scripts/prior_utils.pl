# AUTHOR: Philip Machanick
# CREATE DATE: 28 July 2009

# FIXME: some hacks to work in this directory hierarchy
# including removing hypergeometric prior

# common routines for creating priors

use Scalar::Util qw(looks_like_number);
use strict;
use Data::Dumper;

#requires
push(@INC, split(":", $ENV{'PATH'}));  # look in entire path
require "read_fasta_file.pl"; # for read_fasta_sequence
#require "hypergeometric-dynamicprog.pl"; # for hypergeometric

my ($NAMEPOS, $SEQPOS, $COMMENTPOS) = (0, 1, 2); # order in read_fasta_sequence
my @ALLOWED_MODELS = ('oops', 'zoops', 'tcm');

################################################################################
#
#  normalize_size
#
#  returns as an array of 2 answers sum of normalized sizes of negative and
#  positive sequences respectively
#  Each sequence array is passed as reference to and array of
#  [name, sequencedata, comment]
#  approach: calculate add up mean length / length of each sequence
#
################################################################################

sub normalize_size {
    my ($neg_sequences,    # ref to array of ref to array
        $pos_sequences,    # ref to array of ref to array
        $L_mean            # average over all sequence lengths
    ) = @_;
    
    my ($N_neg_normalized, $N_pos_normalized) = (0,0);
  # recallibrate no. negative sequences to maximum value
  # of adding all normalized scores
  foreach my $seq_data (@$neg_sequences) {
    $N_neg_normalized += $L_mean/length($$seq_data[$SEQPOS]);
  }
  # same for positive sequences
  foreach my $seq_data (@$pos_sequences) {
    $N_pos_normalized += $L_mean/length($$seq_data[$SEQPOS]);
  }
  return ($N_neg_normalized, $N_pos_normalized);
}

################################################################################
#
#  log_background
#
#  convert background model to logs and as a bonus return and array of the
#  minimum and maximum letter probabilies
#
################################################################################

sub log_background {
    my ($backgroundmodel    # reference to hash of letter probabilies, indexed on letter
    ) = @_;
    my ($min_prob, $max_prob);
    foreach my $letter (keys %$backgroundmodel) {
        my $letter_prob = $$backgroundmodel{$letter};
        if (!defined ($min_prob) || $letter_prob < $min_prob) {
            $min_prob = $letter_prob;
        }
        if (!defined ($max_prob) || $letter_prob > $max_prob) {
            $max_prob = $letter_prob;
        }
        $$backgroundmodel{$letter} = log($letter_prob);
    }
    return ($min_prob, $max_prob);
}

################################################################################
#
#  make_background
#
#  returns a cumulative log background model as an array of references to
#  arrays, each of the referenced arrays the log cumulative background for
#  the sequence in the same position in the sequence array.
#  Datatypes:
#  input:   $alphabet: string of allowed letters
#           $sequences: reference to array of references to arrays
#                each of the contained arrays contains
#                [name, sequencedata, comment]
#  returns: reference to array of references to cumulative background
#            models in same order as sequence data, offset by 1 so that
#            position i+1 is sumulative probability for up to letter i
#            in the original sequence
#  Basic method:
#          For each letter in each sequence, look up its background
#          probability as read in by read_background_file and store
#          the log of that probability added to previous logs to give
#          a cumulative total up to that position in the sequence.
#          Allows the background for a w-mer to be computed using
#          one subtraction and one exponential.
#
################################################################################

sub make_background {
    my ($backgroundmodel,    # reference to hash of letter probabilies, indexed on letter
        $sequences           # reference to array of references 
    ) = @_;
    
    my @cumul_back_model = ();

    foreach my $seq_data (@$sequences) {
        $_ = $$seq_data[$SEQPOS];
        my @seqchars = split(//);
        my $last = 0;
        # make thefirst element 0; position j+1 represents cumulative probability up to j
        my @log_cumul_background = ( 0 );
        foreach my $char (@seqchars) {
            die "`$char' not in background" unless defined($$backgroundmodel{$char});
            $last += $$backgroundmodel{$char};
            push (@log_cumul_background, $last);
        }
        push (@cumul_back_model, \@log_cumul_background);
    }
#    print STDERR "cumul background model: ".Dumper(\@cumul_back_model);

  return \@cumul_back_model;
}


################################################################################
#
#  normalized_classifier_scores
#
#  returns a reference to a hash indexed on sequence name containing
#  scores for each name, represented as reference to an array of numbers
#  Each sequence array is passed as reference to and array of
#  [name, sequencedata, comment]
#  approach: calculate add up mean length / length of each sequence
#
#  for triples, we also need the scores for the last $width
#
################################################################################
sub normalized_classifier_scores {
    my ($pos_sequences,     # ref to array each element a sequence (string)
        $pos_dictionary,    # ref to hash indexed on w-mer: count of w-mer in pos set
        $neg_dictionary,    # ref to hash indexed on w-mer: count of w-mer in neg set
        $width,             # current motif width under consideration
        $N_pos,             # size of positive set (# sequences)
        $N_neg,             # size of negative set (# sequences)
        $N_pos_normalized,  # normalized size of positive set
        $N_neg_normalized,  # normalized size of negative set
        $min_score,         # ref to scalar: to return min score from score_w_mer
        $max_score,         # ref to scalar: to return max score from score_w_mer
        $scale,             # if defined scale value to standard range
        $absolute,          # if defined, don't weight for number of sequences
        $type,              # if defined, score discriminative or hypergeometric prior
        $triples,           # if defined score spaced triples rather than w-mers
        $fixed_start,       # if defined score spaced triples with start position fixed
        $ambigs,            # ambiguous characters: score any subsequence with these as 0
        $lastscores         # scores for $width-1 used to avoid recomputing triples scores
    ) = @_;
    
    my %prior_values;
    my $i = 0;
    my %score_cache;

    foreach my $seq_data (@$pos_sequences) {
      my @scores;
      # scores are stored in a hash table hashed on sequence name
      # sequence names and data are stored in an array indexed by $NAMEPOS and $SEQPOS
      my $last_i_scores = defined ($lastscores) ? $$lastscores{@$seq_data[$NAMEPOS]} : undef;
      if ($type && $type eq "D") { # FIXME: only D priors work for protein
          score_discrim_w_mer (\@scores, @$seq_data[$SEQPOS], $pos_dictionary, $neg_dictionary, $width,
            $N_pos, $N_neg, $N_pos_normalized, $N_neg_normalized, $min_score, $max_score, $scale,
            $absolute, $triples, $fixed_start, $ambigs, $last_i_scores, \%score_cache);
      } elsif (defined ($type) && $type eq "hyper") { # also covers case of D-hyper 
          &score_w_mer_hypergeometric (\@scores, @$seq_data[$SEQPOS], $pos_dictionary,
            $neg_dictionary, $width,
            $N_pos, $N_neg, $N_pos_normalized, $N_neg_normalized, $min_score, $max_score, $scale,
            $absolute);  #FIXME: doesn't do triples
      } else { # FIXME: triples not updated to latest stuff in score_discrim_w_mer
          &score_w_mer (\@scores, @$seq_data[$SEQPOS], $pos_dictionary, $neg_dictionary, $width,
            $N_pos, $N_neg, $N_pos_normalized, $N_neg_normalized, $min_score, $max_score, $scale,
            $absolute, $triples, $fixed_start, $last_i_scores, \%score_cache);
      }
      my $name = $$seq_data[$NAMEPOS];
#       print STDERR "scores (total # =". @scores. @$seq_data[$NAMEPOS].":\n".Dumper(\@scores);
      $prior_values{@$seq_data[$NAMEPOS]} = [@scores];
    }
    return \%prior_values;
}

################################################################################
#
#  printScoreAsPrior
#
#  convert scores to prior values; OOPS, ZOOPS and TCM models
#  Datatypes:
#  input:  $allscores: reference to hash of scores indexed on FASTA name
#          each element a reference to array of numbers
#          $model: should be "zoops", "oops" or "tcm"
#          $width: width of motif for which prior computed
#          $scale_min, $scale_max: values used to rescale scores before
#          converting to priors (from [0..1] to [$scale_min .. $scale_max])
#          $outfile: optional path to output file; defaults to STDOUT
#          $names: reference to array of sequence names, use to order
#          output if not undef
#  output: to stdout unless $outfile defined
#          >name 
#  Basic method (oops, zoops):
#          For each sequence, rescale scores to form a probability distribution.
#          For Zoops, K_0 = 1, OOPS K_0 = 0 representing a score
#          for finding no motif in a given sequence.
#          For j = 1..L-w+1: alpha = sum K_j = score_j/(1-score_j)
#          Each prior value p_j = K_j / alpha
#
#          TCM variation: save rescaling until all scores seen and rescale to
#          form a single distirbution
#
#          FIXME: OOPS, TCM not tested
#  
#
################################################################################
sub printScoreAsPrior {
    my ($allscores,
        $model,
        $width,
        $scale_min,
        $scale_max,
        $outfile,
        $names) = @_;
    my $alpha = 0;    # for TCM, don't reset in counting loop
    my %allK;         # only used for TCM
    my $OK;    # test if model string OK:
    foreach my $m (@ALLOWED_MODELS) {
        if ($m eq $model) {
            $OK = 1;
            last;
        }
    }
    $outfile = "-" unless ($outfile); # STDOUT if not defined
    open(OUTFILE, ">$outfile") or die ("failed to open $outfile\n");
    die "invalid model `$model' in printScoreAsPrior" unless defined $OK;
    $names = [sort keys %$allscores] if !$names;
    foreach my $name (@$names) {
      if ($model ne 'tcm') {    # keep accumulating over all seqs for TCM
          $alpha = $model eq 'zoops' ? 1 : 0; # represents K_0 in following summation
      }
      my @K = ();
      my $L = @{$$allscores{$name}}; #$nonzeros = @{$$allscores{$name}} - $width + 1;
      for (my $j = 0; $j < $L; $j++) {
        my $P_ij = ${$$allscores{$name}}[$j];
        push(@K,$P_ij/(1-$P_ij));
        $alpha += $K[-1];  # add on last K
      }
      
      if ($model ne 'tcm') {
          &print_prior_seq ($name, $width, $scale_min, $scale_max, \@K, $alpha, *OUTFILE);
      } else  { #tcm
          $allK{$name} = \@K;
      }
    }
    # we now have a basis for scaling all values over all sequences to a
    # single probability distribution
    if ($model eq 'tcm') {
        foreach my $name (@$names) {
            &print_prior_seq ($name, $width, $scale_min, $scale_max, $allK{$name}, $alpha, *OUTFILE);
        }
    }
}

################################################################################
#
#  print_prior_seq
#  Datatypes:
#  input:  $name: FASTA name of original sequence
#          $width: motif width for prior computed
#          $epsilon: adjustment of original range of scores
#          $K: ref to array of numbers -- priors before final
#          scaling to form probability distribution
#          $alpha: to divide $K values to form probabilities
#          $outfile: optional path to output file; defaults to STDOUT
#  output: stdout PSP file format
#    >name $width epsilon = $epsilon
#    number number ... number
#  Basic method: divides $K values to form probabilities
#  Complexity: O(L) where L is sequence length
#
################################################################################

sub print_prior_seq {
    my ($name,    # name from original FASTA file
    $width,       # width or motif
    $scale_min,   # adjustment to original value range
    $scale_max,
    $K,           # reference to array of partially computed priors
    $alpha,       # value to divide K values to convert to probabilities
    $outfile      # file handle if defined
    ) = @_;
    local *OUTFILE;
    if ($outfile) {
        *OUTFILE = $outfile;
    } else {
        *OUTFILE = \*STDOUT;
    }
    printf OUTFILE '>'.$name." $width scaledmin = %G scaledmax = %G\n", $scale_min, $scale_max;
    foreach my $k_val (@$K) {
      printf OUTFILE "%8.7G ", $k_val / $alpha; # was 20.19lG in older versions
    }
    # now print terminating w-1 zeros
    for (my $i = 0; $i < $width - 2; $i++) {
      print OUTFILE "0 ";
    }
    print OUTFILE "0\n";
}

################################################################################
#
#  print_scores
#
#  print scores after rescaling to the given range
#
################################################################################
sub print_scores {
    my ($allscores,     # reference to hash indexed on names, each element its scores
        $scale_min,     # min scaled value
        $scale_max,     # max scaled value
        $min_score,
        $max_score,
        $width,         # width of motif
        $scale          # if defined, scaled to standard range
    ) = @_;

    # sort the names here to ensure output is in a consistent order
    foreach my $name (sort keys %$allscores) {
      print ">".$name."\n";
      my $scores = $$allscores{$name};
      foreach my $sc (@$scores) {
        if (defined($scale) && $sc > $scale_max) {
          print STDERR "printing $sc > $scale_max max_score = $max_score; min_score = $min_score\n";
        }
        printf "%.19G ", $sc;
      }
      for (my $i = 0; $i < $width - 1; $i++) {
        print "0";
        print " " unless ($i == $width - 2);
      }
      printf "\n";
    }
    
}

################################################################################
#
#  scale_scores
#
#  scale scores to the standard range
#  in this case $allscores is used to return results as well as pass in
#
################################################################################
sub scale_scores {
    my ($allscores,     # reference to hash indexed on names, each element its scores
        $scale_min,     # min scaled value
        $scale_max,     # max scaled value
        $min_score,     # previously determined minimum score
        $max_score,     # previously determined maximum score
        $verbose        # if defined, stats to STDERR
    ) = @_;

    my %scorestats;

    # sort the names here to ensure output is in a consistent order
    foreach my $name (sort keys %$allscores) {
      my $scores = $$allscores{$name};
      for (my $j = 0; $j < @$scores; $j++) {
          my $sc = $$scores[$j];
          die "scaling scores but $max_score <= $min_score" unless ($max_score >= $min_score);
          if ($max_score == $min_score) {
              $sc = $scale_max;
          } else {
              $sc = ($sc - $min_score)/($max_score-$min_score)*($scale_max-$scale_min)+$scale_min;
          }
          die "scaled score $sc out of range [$scale_min..$scale_max]"
              unless ($sc >= $scale_min && $sc <= $scale_max);
        if (defined $verbose) {
          if (exists $scorestats{$sc}) {
            $scorestats{$sc}++;
            } else {
              $scorestats{$sc} = 1;
            }
        }
        $$scores[$j] = $sc;    # save the scaled score
      }
    }
    
    if (defined $verbose) {
      print STDERR "\n";
      foreach my $sc (sort keys(%scorestats)) {
        print STDERR "score $sc occurred $scorestats{$sc} times\n";
      }
    }
}


################################################################################
#
#  class_prior_prob
#
#  use classification of each w-mer in terms of the number of positive set
#  sequences in which it occurs plus the number of negative set sequences in
#  which it does not occur as a fraction of the number of sequences.
#  Datatypes:
#  input:  $pos_sequences: reference to array of positive sequences, each 
#          element a string of letters of the allowed alphabet
#          $neg_sequences: same type as $pos_sequences representing negative
#          sequences
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $width: width of motif under consideration
#          normalizing not required
#  output: $pos_dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#          $neg_dictionary: same for negative sequences
#  Basic method: Normalizes count as average seq length / length of seq in which
#  w-mer found (unless mean length not defined, in which case presence count
#  is 1).
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  Detail at count_all_norm.
#
################################################################################
sub class_prior_prob {
  my ($pos_sequences,   # reference to array of positive sequence data
      $neg_sequences,   # reference to array of negative sequence data
      $revcomp,         # if defined, also count reverse complement of w-mer
      $width,           # motif width
      $pos_dictionary,  # ref to hash count positive seqs where each w-mer seen
      $neg_dictionary,  # ref to hash count negative seqs where each w-mer seen
      $pos_log_cuml_bg, # ref to log cumulative background probabilities for pos seqs
      $neg_log_cuml_bg, # ref to log cumulative background probabilities for pos seqs
      $min_bg_prob      # lowest background probability of any letter
  ) = @_;
    # count the number of sequences in which a w-mer occurs first in the positive
    # then in the negative sequences; the count is weighted for the probabilty
    # of not finding the w-mer based on the background model
    return (&count_all_prob($pos_sequences, $revcomp, $width, $pos_dictionary,
        $pos_log_cuml_bg, $min_bg_prob),
        &count_all_prob($neg_sequences, $revcomp, $width, $neg_dictionary,
        $neg_log_cuml_bg, $min_bg_prob));
}


################################################################################
#
#  count_all_prob FIXME DOCUMENTATION NOT CORRECT
#  FIXME: doesn't do later features like spaced triples
#
#  for each w_mer found in a given sequence, pass in to count_w_mer
#  returns sum of all available scores / number of bases visited
#  
#  Datatypes:
#  input:  $sequences: reference to array of sequences, each a string
#          of letters of the allowed alphabet
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#
#  Basic method: for each sequence, visit each w-mer once and if it
#  has not been encountered before in the current sequence, increment
#  its count (initialize if never seen before), and store the position
#  of the sequence as "last" to remember where it was last seen.
#  Count is 1 in absolute mode or (1 - p_a)^(L-w+1) for background
#  probabiliy of letter a p_a, sequence length L, motif width w) in
#  probabilistic normalization mode.
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  There is variation for w because there is some caching of precomputed
#  matches so a larger w with more unique w-mers will result in slower
#  execution.
#
################################################################################

sub count_all_prob {
  my ($sequences,     # reference to array of sequence data
      $revcomp,       # if defined, also count reverse complement of w-mer
      $width,         # motif width
      $dictionary,    # count of sequences where each w-mer seen
      $log_cuml_bg,   # log cumulative bg for all sequences in this set
      $min_bg_prob    # lowest background probability for any letter
  ) = @_;

  # $i used to ensure we only count a w-mer once per sequence plus index background
  # seq_bg[$j+$w-1] - @seq_bg[$j-1]
  my $i = 0;
  my $maxscore = 0;
  my $min_prob_for_w = 1 - ($min_bg_prob ** $width);
#  print STDERR "in count_all_prob".defined($version)?$version." ":""."\n";
  foreach my $seqdata (@$sequences) {
    my $seq = $$seqdata[$SEQPOS];
    my $bg = @$log_cuml_bg[$i];
    my $L_w_1 = length($seq) - $width + 1; # L-w+1
    for (my $j = 0; $j < $L_w_1; $j++) {
      # subtract cumulative log prob before start of this w-mer then
      # use exp to get probability of this w-mer; cumulative probability
      # vector has a 0 at the front so we don't have to special-case dealing
      # with item 0, i.e. all indexes offset by +1. Cumulative probability for
      # position j+w-1 in sequence is at j+w in cumulative probability
      my $prob = exp ($$bg[$j+$width] - $$bg[$j]);
      my $count = (1-$prob) ** $L_w_1;
      # print STDERR "p = $prob, '(1-$prob)' = (1-$prob), L-w+1 = $L_w_1, s = $count;";
      my $w_mer = substr($seq, $j, $width);
      &count_w_mer ($w_mer, $revcomp, $i, $count, $dictionary);
    }
    # probability for 
    $maxscore += $min_prob_for_w ** $L_w_1;
    $i++;
  }
  return $maxscore;
}


################################################################################
#
#  D_prior
#
#  use discrimination of each w-mer in times it appears in all positive set
#  sequences versus the number of times it occurs in each negative set sequences,
#  as a fraction of the number of sequences.
#  Datatypes:
#  input:  $pos_sequences: reference to array of positive sequences, each 
#          element a string of letters of the allowed alphabet
#          $neg_sequences: same type as $pos_sequences representing negative
#          sequences
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $meanL: number, mean over all sequence lengths (pos+neg); undef if
#          normalizing not required
#  output: $pos_dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#          $neg_dictionary
#  Basic method: counts number of occurrence of each word in positive and negative
#  sets, and calculates a score for each word that is the fraction of all those
#  occurrences that are in the positive set. Here, we do the count; the fraction is
#  calculated in normalized_classifier_scores (so named because other methods may
#  do normalizing).
#  
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  Data structures for counting described at count_all_w_mers.
#
################################################################################
sub D_prior {
  my ($pos_sequences,   # reference to array of positive sequence data
      $neg_sequences,   # reference to array of negative sequence data
      $revcomp,         # if defined, also count reverse complement of w-mer
      $width,           # motif width
      $pos_dictionary,  # ref to hash count positive seqs where each w-mer seen
      $neg_dictionary,  # ref to hash count negative seqs where each w-mer seen
      $triples,         # if defined score spaced triples rather than w-mers
      $alternates       # if defined dictionary of alternatives for given character
  ) = @_;

    # count the number of times each w-mer occurs first in the positive
    # then in the negative sequences
    &count_all_w_mers($pos_sequences, $revcomp, $width, $pos_dictionary, $triples, $alternates);
    &count_all_w_mers($neg_sequences, $revcomp, $width, $neg_dictionary, $triples, $alternates);
} # D_prior


################################################################################
#
#  class_prior_norm
#
#  use classification of each w-mer in terms of the number of positive set
#  sequences in which it occurs plus the number of negative set sequences in
#  which it does not occur as a fraction of the number of sequences.
#  Datatypes:
#  input:  $pos_sequences: reference to array of positive sequences, each 
#          element a string of letters of the allowed alphabet
#          $neg_sequences: same type as $pos_sequences representing negative
#          sequences
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $meanL: number, mean over all sequence lengths (pos+neg); undef if
#          normalizing not required
#  output: $pos_dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#          $neg_dictionary
#  Basic method: Normalizes count as average seq length / length of seq in which
#  w-mer found (unless mean length not defined, in which case presence count
#  is 1).
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  Detail at count_all_norm.
#
################################################################################
sub class_prior_norm {
  my ($pos_sequences,   # reference to array of positive sequence data
      $neg_sequences,   # reference to array of negative sequence data
      $revcomp,         # if defined, also count reverse complement of w-mer
      $width,           # motif width
      $meanL,           # previously computed average sequence length
      $pos_dictionary,  # ref to hash count positive seqs where each w-mer seen
      $neg_dictionary,  # ref to hash count negative seqs where each w-mer seen
      $triples,         # if defined score spaced triples rather than w-mers
      $alternates       # if defined dictionary of alternatives for given character
  ) = @_;

    # count the number of sequences in which a w-mer occurs first in the positive
    # then in the negative sequences; the count is weighted for the length of the
    # sequence in which the w-mer is found
    &count_all_norm($pos_sequences, $revcomp, $width, $meanL, $pos_dictionary, $triples, $alternates);
    &count_all_norm($neg_sequences, $revcomp, $width, $meanL, $neg_dictionary, $triples, $alternates);
}

sub hamming_prior {
  my ($pos_sequences,   # reference to array of positive sequence data
      $neg_sequences,   # reference to array of negative sequence data
      $revcomp,         # if defined, also count reverse complement of w-mer
      $width,           # motif width
      $meanL,           # previously computed average sequence length
      $pos_dictionary,  # ref to hash count positive seqs where each w-mer seen
      $neg_dictionary,  # ref to hash count negative seqs where each w-mer seen
      $triples,         # if defined score spaced triples rather than w-mers
      $alternates,      # if defined dictionary of alternatives for given character
      $min_score,       # reference: return value through this
      $max_score        # reference: return value through this
  ) = @_;

    my %prior_values;
    # score each word for how similar it is to other words in the positive and
    # negative sets. An exact match scores 1; any other match scores matches/w
    for (my $i = 0; $i <  @$pos_sequences; $i++) { # for all positive sequences
      my @scores = ();
      my $seq_data = $$pos_sequences[$i];
      $_ = $$seq_data[$SEQPOS];
      my @seq = split(//);
      my $L_wP1 = @seq - $width + 1;
      # for all words width w in a given sequence
      for (my $j = 0; $j < $L_wP1; $j++) {
          my $j_end = $j + $width - 1;
          my ($pos_count, $neg_count) = (0, 0);
          # find the comparison word and score
          for (my $p = 0; $p <  @$pos_sequences; $p++) {    # for all positive sequences again
              my $comparison;
              if ($p == $i) {    # don't split the array again, reuse
                  $comparison = \@seq;
              } else {
                  my $comparisonseqdata = $$pos_sequences[$p];
                  $_ = $$comparisonseqdata[$SEQPOS];
                  $comparison = [split(//)];
              }
              # compare original + comparison (OK to count self as 1 == all match)
              my $LComp_wP1 = @$comparison - $width + 1;
              my $stop = $LComp_wP1<$L_wP1?$LComp_wP1:$L_wP1;
              # find all words to compare with given word
              my $max_matches = 0;
              for (my $k = 0; $k < $stop; $k++) { # stop based on shorter of 2 sequences
                my $matches = 0;
                my ($k_pos, $k_end) = ($k, $k + $width - 1);
                for (my $pos = $j; $pos < $j_end; $pos++) { 
                  $matches ++ if ($seq[$j] eq $$comparison[$k_pos]);
                  $k_pos++;
                }
                # record the best match count for this word
                $max_matches = $matches if $matches > $max_matches;
              }
              $pos_count += $max_matches/$width;
          }

          for (my $n = 0; $n <  @$neg_sequences; $n++) {    # for all negative sequences now
              my $comparisonseqdata = $$neg_sequences[$n];
              $_ = $$comparisonseqdata[$SEQPOS];
              my $comparison = [split(//)];
              # compare original + comparison (OK to count self as 1 == all match)
              my $matches = 0;
              my $LComp_wP1 = @$comparison - $width + 1;
              my $stop = $LComp_wP1<$L_wP1?$LComp_wP1:$L_wP1;
              # find all words to compare with given word
              my $max_matches = 0;
              for (my $k = 0; $k < $stop; $k++) { # stop based on shorter of 2 sequences
                my $matches = 0;
                my ($k_pos, $k_end) = ($k, $k + $width - 1);
                for (my $pos = $j; $pos < $j_end; $pos++) { 
                  $matches ++ if ($seq[$j] eq $$comparison[$k_pos]);
                  $k_pos++;
                }
                # record the best match count for this word
                $max_matches = $matches if $matches > $max_matches;
              }
              $neg_count += 1 - ($max_matches/$width);
          }
          my $score = ($pos_count + $neg_count) / (@$pos_sequences + @$neg_sequences);
          push (@scores, $score);
          $$min_score = $score if (!$$min_score || $score < $$min_score);
          $$max_score = $score if (!$$max_score || $score > $$max_score);
      }
          $prior_values{@$seq_data[$NAMEPOS]} = [@scores];
    }
    return \%prior_values;
}

################################################################################
#
#  count_all_w_mers
#
#  for each w_mer found in a given sequence, pass in to count_w_mer
#  Datatypes:
#  input:  $sequences: reference to array of sequences, each a string
#          of letters of the allowed alphabet
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $meanL: number, mean over all sequence lengths (pos+neg)
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#  Basic method: for each sequence, visit each w-mer once and increment
#  its count (initialize if never seen before) in a hash table.
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  There is variation for w because there is some caching of precomputed
#  matches so a larger w with more unique w-mers will result in slower
#  execution.
#
################################################################################

sub count_all_w_mers {
  my ($sequences,     # reference to array of sequence data
      $revcomp,       # if defined, also count reverse complement of w-mer
      $width,         # motif width
      $dictionary,    # count of sequences where each w-mer seen
      $triples,       # if defined score spaced triples rather than w-mers
      $alternates     # alternatives for ambiguous characters
  ) = @_;

  my $i = 0; # used to ensure we only count a w-mer once per sequence
  foreach my $seqdata (@$sequences) {
    my $seq = $$seqdata[$SEQPOS];
    my $L_wP1 = length($seq) - $width + 1;
    my $count = 1; # always increment by 1
    for (my $j = 0; $j < $L_wP1; $j++) {
      my $w_mer = substr($seq, $j, $width);
      if ($triples) {
          # FIXME: counting alternative representations not yet implemented
          &count_w_wide_triples ($w_mer, $i, $count, $dictionary, undef);# $alternates
      } else { # count_w_mer_triples for don't cares in any position FIXME: put in alternates
          &count_w_mer ($w_mer, $revcomp, $i, $count, $dictionary);
      }
      $i++;    # finagle the "last seen" index so we count every instance
    }
  }
}


################################################################################
#
#  count_all_norm
#
#  for each w_mer found in a given sequence, pass in to count_w_mer
#  Datatypes:
#  input:  $sequences: reference to array of sequences, each a string
#          of letters of the allowed alphabet
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $meanL: number, mean over all sequence lengths (pos+neg)
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#  Basic method: for each sequence, visit each w-mer once and if it
#  has not been encountered before in the current sequence, increment
#  its count (initialize if never seen before), and store the position
#  of the sequence as "last" to remember where it was last seen.
#  Count is 1 in absolute mode or (mean sequence length) / (this sequence
#  length) in normalized mode.
#  Complexity: O(LN) where L is sequence length and N # sequences.
#  There is variation for w because there is some caching of precomputed
#  matches so a larger w with more unique w-mers will result in slower
#  execution.
#
################################################################################

sub count_all_norm {
  my ($sequences,     # reference to array of sequence data
      $revcomp,       # if defined, also count reverse complement of w-mer
      $width,         # motif width
      $meanL,         # previously computed average sequence length
      $dictionary,    # count of sequences where each w-mer seen
      $triples,       # if defined count spaced triples rather than w-mers
      $alternates     # alternatives for ambiguous characters
  ) = @_;

  my $i = 0; # used to ensure we only count a w-mer once per sequence
  foreach my $seqdata (@$sequences) {
    my $seq = $$seqdata[$SEQPOS];
    my $count = defined($meanL) ? $meanL/length($seq) : 1;
    for (my $j = 0; $j < length($seq) - $width + 1; $j++) {
      my $w_mer = substr($seq, $j, $width);
      if ($triples) {
          # FIXME: counting alternative representations not yet implemented
          &count_w_wide_triples ($w_mer, $i, $count, $dictionary, undef);# $alternates
      } else { # count_w_mer_triples for don't cares in any position FIXME: put in alternates
          &count_w_mer ($w_mer, $revcomp, $i, $count, $dictionary);
      }
    }
    $i++;
  }
}


################################################################################
#
#  count_w_mer
#
#  for each w_mer found in a given sequence, increment dictionary count
#  (or initialize if not already present) unless already seen in the same
#  sequence: relies on processing the sequences in order. 
#  Datatypes:
#  input:  $w_mer: string containing w-mer to be counted
#          $revcomp: undef if reverse complement not required, any other
#          value otherwise
#          $i: index of current sequence to avoid double-counting
#          $count: amount to increment score
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#  Basic method: Value of each count is 1 if abs mode otherwise normalized
#  by seq length.
#  Whichever of the w-mer or reverse complement of a w-mer is found first
#  will set the dictionary entry; from there on, only that entry will be
#  the same (if $revcomp defined) because it is stored as a reference in
#  bother entries each time (again if $revcomp defined).
#
################################################################################

sub count_w_mer {
  my (
      $w_mer,        # the w_mer we are currently counting
      $revcomp,      # if defined, also count if revcomp($w_mer present)
      $i,            # index unique to a sequence to enforce single counts
      $count,        # 1 or scaled score depending on score method
      $dictionary    # indexed on w-mer, score and last i per entry
  ) = @_;
  my $entry;
  my $found;
  my $reverse = &reverse_complement($w_mer);
  # score 1 or scaled score for each time present in $seq
  if (exists ($dictionary->{$w_mer})) {
    $entry = $dictionary->{$w_mer};
    &count_present ($entry, $i, $count);
    $found = 1;
  # for revcomp case count if exists as reverse complement
  } elsif (defined($revcomp)) {
    if (exists ($dictionary->{$reverse})) {
      $entry = $dictionary->{$reverse};
      &count_present ($entry, $i, $count);
      $found = 1;
    }
  }
  if (!defined($found)) {
    my %new_entry;
    $new_entry{"count"} = $count;
    $new_entry{"last"} = $i;
    $entry = \%new_entry;
  }
  # here we don't care if we found it as forward or revcomp;
  # store a reference to the entry, so it's updated in either
  # case; we want it stored at the location where the actual
  # w-mer was found so it will be found again when using the
  # scores. When we hit the first revcomp instance it will be
  # copied there too.
  $dictionary->{$w_mer} = $entry;
  if (defined($revcomp)) {
      $dictionary->{$reverse} = $entry;
  }
}


################################################################################
#
#  count_w_mer_triples
#
#  for each w-mer found in a given sequence, increment dictionary count
#  (or initialize if not already present) unless already seen in the same
#  sequence: relies on processing the sequences in order. Intended for protein,
#  this variation counts all triples within a w-mer, ignoring letters between
#  the triples, e.g., for VDFMVQSPNEKP, with w=4, would count triples
#  VDF. VD.M V.FM .DFM in the first w-mer where "." == don't care.
#  Datatypes:
#  input:  $w_mer: string containing w-mer to be counted
#          $i: index of current sequence to avoid double-counting
#          $count: amount to increment score
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#  Basic method: Value of each count is 1 if abs mode otherwise normalized
#  by seq length.
#  Whichever of the w-mer or reverse complement of a w-mer is found first
#  will set the dictionary entry; from there on, only that entry will be
#  the same (if $revcomp defined) because it is stored as a reference in
#  bother entries each time (again if $revcomp defined).
#
################################################################################

sub count_w_mer_triples {
  my (
      $w_mer,        # the w_mer we are currently counting
      $last,         # index unique to a sequence to enforce single counts
      $count,        # 1 or scaled score depending on score method
      $dictionary    # indexed on w-mer, score and last i per entry
  ) = @_;
  my @letters = split(//,$w_mer);
  my $width = @letters;
  my ($w_1, $w_2) = ($width-1, $width-2);
  
  for (my $i = 0; $i < $w_2; $i++) {
    for (my $j = $i+1; $j < $w_1; $j++) {
        for (my $k = $j+1; $k < $width; $k++) {
            my ($entry, $found);
            # separate i,j,k with "," to avoid need to worry about digits/number
            my $key = $i.",".$j.",".$k.":$letters[$i].$letters[$j].$letters[$k]";
            if (exists($dictionary->{$key})) {
                $entry = $dictionary->{$key};
                &count_present ($entry, $last, $count);
                $found = 1;          
            }
            if (!defined($found)) {
              my %new_entry;
              $new_entry{"count"} = $count;
              $new_entry{"last"} = $last;
              $entry = \%new_entry;
            }
            $dictionary->{$key} = $entry;
        }
    }
  }

}

################################################################################
#
#  count_w_wide_triples
#
#  for each w-wide triple in a given sequence, increment dictionary count
#  (or initialize if not already present) unless already seen in the same
#  sequence: relies on processing the sequences in order. Intended for protein,
#  counts all triples within a w-mer starting and ending respectively at the
#  start and end point of the w-mer, ignoring letters between.
#  Example: for VDFMVQSPNEKP, with w=4, would count triples
#  VD.M and V.FM in the first w-mer where "." == don't care.
#  Datatypes:
#  input:  $w_mer: string containing w-mer to be counted
#          $i: index of current sequence to avoid double-counting
#          $count: amount to increment score
#  output: $dictionary: reference to hash indexed on w-mer, with entries
#          each also a hash containing values indexed on "count" (count
#          of that w-mer) and "last" (index of previous place w-mer seen)
#  Basic method: Value of each count is 1 if abs mode otherwise normalized
#  by seq length.
#  Whichever of the w-mer or reverse complement of a w-mer is found first
#  will set the dictionary entry; from there on, only that entry will be
#  the same (if $revcomp defined) because it is stored as a reference in
#  bother entries each time (again if $revcomp defined).
#
################################################################################

sub count_w_wide_triples {
  my (
      $w_mer,        # the w_mer we are currently counting
      $last,         # index unique to a sequence to enforce single counts
      $count,        # 1 or scaled score depending on score method
      $dictionary,   # indexed on w-mer, score and last i per entry
      $alternates    # alternative representations for some characters
  ) = @_;
  
  # FIXME: handling alternates not yet in place: changes them to an X
  # an X also should be handled as something other than an exact match
  my @letters = split(//,$w_mer);
  my $w_1 = @letters-1;
  my ($start, $end) = ($letters[0], $letters[-1]);
  for (my $j = 1; $j < $w_1; $j++) {
      my ($entry, $found);
      my $middle = $letters[$j];
      if ($alternates) {
           if(exists $$alternates{$start}) {
               $start = "X"; # FIXME $$alternates{$start};
           }
           if (exists $$alternates{$end}) {
               $end = "X"; # FIXME $$alternates{$end};
           }
           if (exists $$alternates{$middle}) {
               $middle = "X"; # FIXME $$alternates{$middle};
           }
      }
      my $key = "$j,$w_1:$start.$letters[$j].$end";
      if (exists($dictionary->{$key})) {
          $entry = $dictionary->{$key};
          &count_present ($entry, $last, $count);
          $found = 1;          
      }
      if (!defined($found)) {
        my %new_entry;
        $new_entry{"count"} = $count;
        $new_entry{"last"} = $last;
        $entry = \%new_entry;
      }
      $dictionary->{$key} = $entry;
  }
}

#&count_w_mer_triples ($w_mer, $last, $count, $dictionary);

################################################################################
#
#  count_present
#
#  Given an entry in the dictionary, count it if not seen previously in
#  this sequence. If you want to count all occurrences, imcrement $last
#  for each call.
#
################################################################################
sub count_present {
  my (
      $entry,    # dictionary entry to update, reference to hash
      $last,     # current sequence number
      $count     # amount to increment the score
  ) = @_;
  if ($$entry{"last"} != $last) { # only count once per seq
    $$entry{"last"} = $last;
    $$entry{"count"} += $count;
  }
}


################################################################################
#
#    score_w_mer_hypergeometric  FIXME: omitted from this version
#
#    for each position for the w_mer starting there, look its score up in the
#    positive and negative dictionaries and calculate the overall score  as
#    (log(p)) where
#    p = hypergeometric (score_pos, score_pos+score_neg, N_pos, N_pos+N_neg)
#    This classifies each w-mer as to how strongly it discriminates the
#    positive sequences.
#    FIXME: slow so not implemented for later features like spaced triples
#
################################################################################

sub score_w_mer_hypergeometric {
#  print STDERR "\n";
}

################################################################################
#
#    extend_w_mer
#
#    given scored triples, record pairs representing positions 1,2 and 1,3 with
#    the missing letter -- EXPERIMENTAL, NOT USED
#
################################################################################

sub extend_w_mer {
  my ($seq,              # string containing sequence
      $dictionary,       # for each w_mer in given set, how many sequences it's in
                         # key is w-mer, entries each a hash: "count" relevant here
      $width,            # width of w_mer under consideration
      $N_pos,            # number of positive sequences
      $N_neg,            # number of negative sequences
      $N_pos_normalized, # number of positive sequences normalized for maximum total count
      $N_neg_normalized, # number of negative sequences normalized for maximum total count
      $min_score,        # reference to scalar: highest score calculated here
      $max_score,        # reference to scalar: lowest score calculated here
      $scale,            # if defined, calculate $min_score, $max_score
      $absolute,         # if defined, don't weight for number of sequences
      $triples,          # if defined score spaced triples rather than w-mers
      $fixed_start       # if defined score spaced triples with start position fixed
    ) = @_;

    my %part_dictionary;
    foreach (keys %$dictionary) {
        # split the key into midps,endpos and A.B.C
        my @keyparts = split(/:/);
        $_ = $keyparts[1];
        my @letters = split(/\./);
        $_ = $keyparts[0];
        my @positions = split(/,/);
#print STDERR "key parts `@keyparts', letters = `@letters', positions = `@positions'\n";
#exit(1);
        my $key = "1,2#".$positions[0].":$letters[0].$letters[1]";
        $part_dictionary{$key} = exists($part_dictionary{$key}) ? $part_dictionary{$key}+1 : 1;
        $key = "1,3#".$positions[1].":$letters[0].$letters[2]";
        $part_dictionary{$key} = exists($part_dictionary{$key}) ? $part_dictionary{$key}+1 : 1;
    }
#print STDERR "extended dictionary: ".Dumper(\%part_dictionary);
}

################################################################################
#
#    score_discrim_w_mer
#
#    for D-prior score
#
#    for each position for the w_mer starting there, look its count up in the
#    positive and negative dictionaries and calculate the overall score as
#    (pos count + [1]) / (pos count + [1] + neg count + [1 + m/n]) where values
#    in [] are pseudocounts to avoid false peaks for rarely appearing w-mers
#    where m = no. negative sequences, n = no. positive sequences
#
#    for triples, use $last_i_scores, the scores for $width-1, to cut repetition
#    stored in a reference to an array
#
################################################################################

sub score_discrim_w_mer {
  my ($scores,           # reference to array to contain scores for this sequence
      $seq,              # string containing sequence
      $pos_dictionary,   # for each w_mer in positive set, how often it occurs in total
                         # key is w-mer, entries each a hash: "count" relevant here
      $neg_dictionary,   # for negative set, same data as positive set
      $width,            # width of w_mer under consideration
      $N_pos,            # number of positive sequences
      $N_neg,            # number of negative sequences
      $N_pos_normalized, # number of positive sequences normalized for maximum total count
      $N_neg_normalized, # number of negative sequences normalized for maximum total count
      $min_score,        # reference to scalar: highest score calculated here
      $max_score,        # reference to scalar: lowest score calculated here
      $scale,            # if defined, calculate $min_score, $max_score
      $absolute,         # if defined, don't weight for number of sequences
      $triples,          # if defined score spaced triples rather than w-mers
      $fixed_start,      # if defined score spaced triples with start position fixed
      $ambigs,           # ambiguous characters: if any score subsequence as 0
      $last_i_scores,    # scores for $width-1 for this sequence
      $score_cache       # previous scores for triples of this width (ref to hash)
    ) = @_;

  my $L_wP1 = length($seq) - $width + 1; # L-w+1
  my ($pseudocount_enum, $pseudocount_denom) = (1, 1 + $N_neg/$N_pos);
  for (my $j = 0; $j < $L_wP1 ; $j++) {
    my $w_mer = substr($seq, $j, $width);
    my @letters = split(//,$w_mer);
    my $w_1 = @letters-1;
    my $score;
    if ($triples) {
        my $w_mer_max_score = 0;
        if ($last_i_scores) {
            # if we have scores for $w-1 initialize max as max from
            # positions $j and $j+1; $j+1 not off end of data because
            # of loop endpoint
            $w_mer_max_score = @$last_i_scores[$j];
            # if fixedstart mode, don't consider tuples starting in previous
            # scores' j+1 position
            if (!$fixed_start && @$last_i_scores[$j+1] > $w_mer_max_score) {
                $w_mer_max_score = @$last_i_scores[$j+1];
            }
        }
        # base case: for $width = 3, $j loop will go only once
        # and $w_mer_max_score will not be defined at start
        my $i = 0;
        my $k = $width-1;
            for (my $j = $i+1; $j < $width-1; $j++) {
                    my $key = ($j-$i).",".($k-$i).":$letters[$i].$letters[$j].$letters[$k]";
                    if ($ambigs && $key =~ /$ambigs/) {
                        $score = 0;
                    } elsif (exists($pos_dictionary->{$key}) && exists($$score_cache{$key})) {
                        $score = $$score_cache{$key};
                    } else {
                        my $pos_total =
                            exists($pos_dictionary->{$key})?$pos_dictionary->{$key}{"count"}:0;
                        my $neg_total =
                            exists($neg_dictionary->{$key})?$neg_dictionary->{$key}{"count"}:0;
                        $score =  ($pos_total + $pseudocount_enum) / ($pos_total + $neg_total +
                            $pseudocount_denom);
                        $$score_cache{$key} = $score;
print STDERR "score $score > 1 : pos = $pos_total, neg =  $neg_total, pce = $pseudocount_enum, pcd = $pseudocount_denom \n" if ($score > 1);
                    }
                    if (! $w_mer_max_score || $score > $w_mer_max_score) {
                        $w_mer_max_score = $score;
print STDERR "score $score > 1 : pce = $pseudocount_enum, pcd = $pseudocount_denom \n" if ($score > 1);
                    }
                }
        $score = $w_mer_max_score;
    } else {
        if ($ambigs && $w_mer =~ /$ambigs/) {
            $score = 0;
        } else {
            my $pos_total = (exists($pos_dictionary->{$w_mer})?$pos_dictionary->{$w_mer}{"count"}:0);
            my $neg_total = (exists($neg_dictionary->{$w_mer})?$neg_dictionary->{$w_mer}{"count"}:0);
            # add pseudocounts to the pos_total / ($pos_total + $neg_total) score here
            $score =  ($pos_total + $pseudocount_enum) / ($pos_total + $neg_total + $pseudocount_denom);
        }
    }
    if (1) { #defined ($scale)
      if ((! defined($$min_score)) || $score < $$min_score) {
          $$min_score = $score;
      }
      if (! defined($$max_score) || $score > $$max_score) {
          $$max_score = $score;
      }
    }
    push(@$scores, $score);
  }
}

################################################################################
#
#    score_w_mer
#
#    for each position for the w_mer starting there, look its score up in the
#    positive and negative dictionaries and calculate the overall score  as
#    (pos count + ( #neg seqs - neg count)) / (# pos seqs + # neg seqs)
#    normalized for how long each sequence in which the w_mer was found and how
#    many sequences there are in the positive and negative sets.
#    This classifies each w-mer as to how strongly it discriminates the
#    positive sequences.
#
################################################################################

sub score_w_mer {
  my ($scores,           # reference to array to contain scores for this sequence
      $seq,              # string containing sequence
      $pos_dictionary,   # for each w_mer in positive set, how many sequences it's in
                         # key is w-mer, entries each a hash: "count" relevant here
      $neg_dictionary,   # for negative set, same data as positive set
      $width,            # width of w_mer under consideration
      $N_pos,            # number of positive sequences
      $N_neg,            # number of negative sequences
      $N_pos_normalized, # number of positive sequences normalized for maximum total count
      $N_neg_normalized, # number of negative sequences normalized for maximum total count
      $min_score,        # reference to scalar: highest score calculated here
      $max_score,        # reference to scalar: lowest score calculated here
      $scale,            # if defined, calculate $min_score, $max_score
      $absolute,         # if defined, don't weight for number of sequences
      $triples,          # if defined score spaced triples rather than w-mers
      $fixed_start,      # if defined score spaced triples with start position fixed
      $last_i_scores,    # FIXME: not used here: needed to speed up triples
      $score_cache       # FIXME: check if all OK
    ) = @_;

  my $Neg_weight = $N_pos/$N_neg;

  if (defined($absolute)) {
    $Neg_weight = 1;
  }
  my $L = length($seq);
  for (my $j = 0; $j < $L - $width + 1; $j++) {
    my ($neg_total, $pos_total) = (0,0);
    my $w_mer = substr($seq, $j, $width);
    my @letters = split(//,$w_mer);
    my $w_1 = @letters-1;
    my $w_mer_max_score;
    my $score;

#FIXME: triples do not include all the updates and performance enhancements of
#score_discrim_w_mer##########################################################

    if ($triples) {
        my $last = $fixed_start?1:$width-2;
        for (my $i = 0; $i < $last; $i++) {
            for (my $j = $i+1; $j < $width-1; $j++) {
                for (my $k = $j+1; $k < $width; $k++) {
                    my $key = ($j-$i).",".($k-$i).":$letters[$i].$letters[$j].$letters[$k]";
                    # match for spaced triples that can start with don't cares
#                    my $key = $i.",".$j.",".$k.":$letters[$i].$letters[$j].$letters[$k]";
                    if (exists($pos_dictionary->{$key}) && exists($$score_cache{$key})) {
                        $score = $$score_cache{$key};
                    } else {
                        $pos_total = exists($pos_dictionary->{$key})?$pos_dictionary->{$key}{"count"}:0;
                        $neg_total = exists($neg_dictionary->{$key})?$neg_dictionary->{$key}{"count"}:0;
                        $score =  ($pos_total + ($N_neg_normalized-$neg_total)*$Neg_weight) /
                            ($N_pos_normalized + $N_neg_normalized*$Neg_weight);
                        $$score_cache{$key} = $score;
                    }
                    if (! defined($w_mer_max_score) || $score > $w_mer_max_score) {
                        $w_mer_max_score = $score;
                    }
                    
                }
            }
        }
    } else {
        #    in the "absolute" case reduces to
        #    $score = ($pos_total+$N_neg - $neg_total) / ($N_pos +$N_neg);
        $pos_total = exists($pos_dictionary->{$w_mer})?$pos_dictionary->{$w_mer}{"count"}:0;
        $neg_total = exists($neg_dictionary->{$w_mer})?$neg_dictionary->{$w_mer}{"count"}:0;
        $w_mer_max_score =  ($pos_total + ($N_neg_normalized-$neg_total)*$Neg_weight) /
            ($N_pos_normalized + $N_neg_normalized*$Neg_weight);

    }
    if (defined ($scale)) {
      if (! defined($$min_score) || $w_mer_max_score < $$min_score) {
          $$min_score = $w_mer_max_score;
      }
      if (! defined($$max_score) || $w_mer_max_score > $$max_score) {
          $$max_score = $w_mer_max_score;
      }
    }
    push(@$scores, $w_mer_max_score);
  }
}


################################################################################
#
#  check_numeric
#
#  check if value is  defined and numeric; type, max and min optional
#  if supplied type "int" signals check if truncated == original
#  returns 0 if tests fail, 1 otherwise
#
################################################################################

sub check_numeric {
  my ($data,    # string to check
      $type,        # if defined and "int" check that no fraction otherwise ignore
      $min,         # if defined, check that value >= $min
      $max          # if defined, check that value <= $max
  ) = @_;
  if (defined($data)) {
        if (!looks_like_number($data)) {
          return 0;
        }
        if (defined($max) && $data > $max) {
          return 0;
        }
        if (defined($min) && $data < $min) {
          return 0;
        }
        if (defined($type) && $type eq "int") {
          return int($data) == $data;
        }
        # passed all tests here:
        return 1;
  } else {
    return 0;
  }
}

################################################################################
#
#  reverse_complement
#
#  return reverse complement of DNA string: assumes from right alphabet
#  and all in caps (any other letters reversed without translation).
#
################################################################################

sub reverse_complement {
  ($_) = @_;
  tr/ACGT/TGCA/;
  $_ = reverse;
  return $_;
}

################################################################################
#
#  roundup
#
#  return smallest integer >= given value
#
################################################################################

sub roundup {
    my ($n) = @_;
    return(($n == int($n)) ? $n : int($n + 1))
}


################################################################################
#
#    read_seq_file
#
#    read a sequence file and return the sequences converted to upper case plus
#    names and comments after names $sequences. If alphabet given die if letter
#    not in given alphabet string encountered.
#
################################################################################

sub read_seq_file {
  my ($infilename,   # string containing path to file
      $sequences,    # reference to array to contain the result (uppercased)
      $total_length, # reference to scalar to contain length of new sequence
      $alphabet,     # if defined used to check for valid letters
      $translations, # if each item in array ref a translation
      $verbose       # if defined print stats to STDERR
  ) = @_;
  $$total_length = 0;
  if (defined $verbose) { print STDERR "$infilename: "; }
  open(INFILE, "<$infilename") or die ("failed to open $infilename\n");
  my $firsttime = 1; # tell read_fasta_sequence new file started

  while (my @seq = read_fasta_sequence(\*INFILE, $firsttime)) {
    $firsttime = undef; # tell read_fasta_sequence file not new anymore
    $_ = $seq[$SEQPOS];
    $seq[$SEQPOS] = uc;
    if (defined $alphabet) {
      die "sequence data `$seq[$SEQPOS]' not from alphabet `$alphabet'"
        unless ($seq[$SEQPOS] =~ m/^[$alphabet]+$/);
    }
    if ($translations && @$translations) { # array ref points to nonempty array
   	    $_ = $seq[$SEQPOS];
        foreach my $translation (@$translations) {
    	    eval "tr/$$translation[0]/$$translation[1]/";
    	}
   	    $seq[$SEQPOS] = $_;
    }
    push (@$sequences, \@seq);
    # seq[0] is name, $seq[1] is sequence, $seq[2] is comment
    $$total_length += length($seq[$SEQPOS]);
  }
  die "no sequence data in $infilename" unless $$total_length > 0;
  close INFILE or die "failed to close $infilename";
  if (defined $verbose) { print STDERR "$$total_length bases or amino acids\n"; }
}


################################################################################
#
#		list_w_mers
#
#       return all unique w_mers in a set of sequences as a list
#		optionally including reverse complements
#
################################################################################

sub list_w_mers {
	my ($sequences, $revcomp, $width) = @_;
	my %dictionary;
	my $N = @$sequences;
	my @w_mers;
	foreach my $seq (@$sequences) {
		my $reverse;
		if (defined($revcomp)) {
			$reverse = reverse_complement($seq);
		}
		my $M = length($seq);
		for (my $j = 0; $j < $M-$width+1; $j++) {
			$dictionary{substr($seq, $j, $width)} = 1;
			if (defined($revcomp)) {
				$dictionary{substr($reverse, $j, $width)} = 1;
			}
		}
	}
	for my $name (keys %dictionary) {
		push (@w_mers, $name);
	}
	return \@w_mers;
}

################################################################################
#
#       W-mer to PSPM
#
#  convert a w-mer to a PSP in MEME format,
#
################################################################################
sub wmer2pspm {
  my $letters;
  ($_, $letters) = @_;
  my @bases = split(//);
  my @matrix;
  my $lastposition = @$letters[scalar@{$letters}-1];
  foreach my $base (@bases) {
    my $line = "";
    foreach my $position (@$letters) {
      #print "considering `$base' at `$position'\n";
      $line .= (uc($base) eq $position)?1:0;
      $line .= ($position eq $lastposition)?"":" ";
    }
    push(@matrix,$line);
  }
  return \@matrix;
}

################################################################################
#  read_prior_file
#
#  Read a PRIOR file and return the values in an
#  associative array indexed by sequence name (based on read_fasta_file).
#  Values in the associative array are in a reference to an anonymous array.
#  Does NOT check validity of data to allow use for a range of prior types.
#  Only assumes a FASTA-like name line starting with >name followed by a
#  sequence of values that may contain spaces (unlike read_fasta_file, which
#  edits out spaces)
#
#  if $keepcomments is defined then return a hash containing references
#  to arrays of references to arrays of PSP, comments (from end of name line)
#
#  $file    open file pointer
#
################################################################################
sub read_prior_file {
  my ($file, $keepcomments) = @_;

  open(F, "<$file") || die("Can't open file $file\n");
  
  my ($name, $seq, $comments) = ("","","");
  my %seqs;

  # read the file
  while (<F>) {
    if (/^>(\S+)\s+(\S.*)$/) {        # found start of sequence
      if ($name) {        # save previous sequence in array
          $_ = $seq;
          if ($keepcomments) {
              $seqs{$name} = [[split], $comments];
          } else {
              $seqs{$name} = [split];
          }
      }
      $name = $1;        # save sequence name
      $comments = $2;
      $seq = "";
    } else {
      $seq .= $_;
    }
  } # read file
  if ($name && $seq) {        # save very last sequence in array
    $_ = $seq;
    if ($keepcomments) {
        $seqs{$name} = [[split], $comments];
    } else {
        $seqs{$name} = [split];
    }
  }

  %seqs;          # return associative array of names and arrays of values
} # read_prior_file

1;
