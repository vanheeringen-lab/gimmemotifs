#!@WHICHPERL@
# AUTHOR: Philip Machanick
# CREATE DATE: 23 July 2009

# Calculate discriminative prior based on Hartemink idea
# in format for use with MEME
# options (currently inactive) for other prior formats

# input: positive set (X) and negative set (Y) as FASTA file
# outout: stdout

use Scalar::Util qw(looks_like_number);
use strict;
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename;

my $PGM = $0;      # name of program
$PGM =~ s#.*/##;                # remove part up to last slash
#@args = @ARGV;      # arguments to program
my @args = ("-pos filename", "-neg filename");

#
# defaults
#
# sequence alphabet choices
#protein can include U (Sec) selenocysteine http://www.chem.qmul.ac.uk/iubmb/newsletter/1999/item3.html
my $PROT_ALPHABET = "ACDEFGHIKLMNPQRSTUVWY"; 
my $DNA_ALPHABET = "ACGT";    
# set the alphabet: change if protein specified
my $ALPHABET = $DNA_ALPHABET;

my $scale_max; # used to scale if -scale set to min..max
my $scale_min;
my $scale = 1; # turn off to stop scaling (DEBUG)

my $min_triple = 3; # for fixed endpoint triples, start counting loop uses this instead of min_w

$| = 1;        # flush after all prints
$SIG{'INT'} = \&cleanup;  # interrupt handler
# Note: so that interrupts work, always use for system calls:
#   if ($status = system("$command")) {&cleanup($status)}
my $script_loc = dirname(abs_path($0));
# requires
push(@INC, split(":", $ENV{'PATH'}));  # look in entire path
push(@INC, "@BINDIR@"); # add in location of installed binaries
push(@INC, $script_loc); # add in script directory
require 'prior_utils.pl';

my ($min_w, $max_w) = (4,20); # defaults

my $usage = <<USAGE;    # usage message
USAGE: $PGM @args [options]

  Options:
     -minw W        minimum width of motif to consider
                    default: $min_w
     -maxw W        maximum width of motif to consider
                    default: $max_w
     -alpha DNA|prot
                    choose whether DNA or protein sequences
                    default: DNA
     -triples       use spaced triples for matches
                    default: do exact matches of w-mers
     -fixedstart    for triples, anchor the start when scoring
                    triples of width < w
                    default: allow start to move
     -equiv "letters{-letters}"
                    each set of letters separated by a dash is
                    treated as equivalent in searches; sets
                    may not overlap
                    default: no letters treated as equivalent
     -arbitraryOK   allow N in DNA or X in protein
                    default: do not allow N or X
     -revcomp       count reverse complements in computing scores
                    default: only count forward matches
     -scalemin <number>
                    scale scores to mimumum <number>
                    default: 0.1 or 1-scalemax if set
     -scalemax <number>
                    scale scores to maximum of <number>
                    default: 0.9 or 1-scalemin if set
     -maxrange      instead of choosing W with maximum score choose W with
                    maximum difference between maximum and minimum scores
     -raw           output scores instead of priors
     -reportscores  report pos and neg file names, min and max scores,
                    min and max w : tab-separated to STDERR
     -verbose       send stats to stderr reporting frequency of each
                    score
                    default: do not report statistics
               
  Compulsory:
     -pos filename  sequences likely to contain binding sites
     -neg filename  sequences unlikely to contain binding sites

  Calculates a positional prior by classifying each position of width W
  as to how strongly it fits the positive set versus the negative set.
  For each sequence:
  >name scaledmin = 0.1 scaledmax = 0.9
  prior probability for each position in the sequence
  
  The prior probabilities are each a single number, one per letter of the
  sequence data. The actual values of the scaled minimum and maximum may
  vary if either or both -scalemin and -scalemax are set.
  
  Each input file should be in FASTA format. Anything after the name on
  the name line is echoed to output. The name has appended to it width
  W chosen for the sequences' prior.
    
  Reads -pos FASTA file and -neg FASTA file.
  Writes stdout. If -verbose is used, writes stderr.

        Copyright
        (2009) The University of Queensland
        All Rights Reserved.
        Author: Philip Machanick
USAGE

my $nargs = @args;      # number of required args
if ($#ARGV+1 < $nargs) { &print_usage("$usage", 1); }

my ($pos_filename, $neg_filename);
my ($priortype, $d_prior, $d_hyper_prior, $hypergeomtric_prior, $hamming_prior,
   $revcomp, $absolute, $maxrange, $raw, $report, $NorX, $triples, $fixed_start, $equivs);
my $verbose; # define to get stats out
# get input arguments
my $alpha_choice;

# if allowing a -type flag, set these based on command line
$d_prior = 1;
$priortype = "D";

while ($#ARGV >= 0) {
  $_ = shift;
  if ($_ eq "-h") {        # help
    &print_usage("$usage", 0);
  } elsif ($_ eq "-minw") {
    $min_w = shift; # smallest window size
    &print_usage("$usage", 1) unless (check_numeric($min_w, "int", 0));
  } elsif ($_ eq "-revcomp") {
    $revcomp = 1;
  } elsif ($_ eq "-maxw") {
    $max_w = shift; # biggest window size
    &print_usage("$usage", 1) unless (check_numeric($max_w, "int", 0));
  } elsif ($_ eq "-alpha") {
    $alpha_choice = shift; # alphabet type
    $alpha_choice = uc $alpha_choice; # don't care about capitalization
    if ($alpha_choice eq "PROT") {
        $ALPHABET = $PROT_ALPHABET;
    } elsif ($alpha_choice ne "DNA") {
        &print_usage("$usage", 1);
    }
  } elsif ($_ eq "-equiv") {
    $equivs = shift;
    &print_usage("$usage", 1) unless $equivs;
  } elsif ($_ eq "-arbitraryOK") {
    $NorX = "";
  } elsif ($_ eq "-triples") {
    $triples = 1;
  } elsif ($_ eq "-fixedstart") {
    $fixed_start = 1;
  } elsif ($_ eq "-pos") {
    $pos_filename = shift;
    &print_usage("$usage", 1) unless $pos_filename;
  } elsif ($_ eq "-neg") {
    $neg_filename = shift;
    &print_usage("$usage", 1) unless $neg_filename;
  } elsif ($_ eq "-scalemin") {
      $scale_min = shift;
      &print_usage("$usage", 1) unless (check_numeric($scale_min, "", 0,1));
  } elsif ($_ eq "-scalemax") {
      $scale_max = shift;
      &print_usage("$usage", 1) unless (check_numeric($scale_max, "", 0,1));
  } elsif ($_ eq "-maxrange") {
    $maxrange = 1;
  } elsif ($_ eq "-raw") {
    $raw = 1;
  } elsif ($_ eq "-reportscores") {
    $report = 1;
  } elsif ($_ eq "-verbose") {
    $verbose = 1;
  } else {
    &print_usage("$usage", 1);
  }
}
# give up if any of the compulsory command line args is not set
&print_usage("$usage", 1) unless (defined($pos_filename) && defined($neg_filename));
# check that W values make sense
&print_usage("$usage", 1) unless ($min_w <= $max_w);
## reverse complement and protein don't mix
&print_usage("$usage", 1) if (defined($revcomp) && $ALPHABET eq $PROT_ALPHABET);
# give up if any of the compulsory command line args is not set
&print_usage("$usage", 1) unless (defined($pos_filename) && defined($neg_filename));
# check that W values make sense
&print_usage("$usage", 1) unless ($min_w <= $max_w);
## reverse complement and protein don't mix
&print_usage("$usage", 1) if (defined($revcomp) && $ALPHABET eq $PROT_ALPHABET);

# give up if any of the compulsory command line args is not set
&print_usage("$usage", 1) unless (defined($pos_filename) && defined($neg_filename));
# check that W values make sense
&print_usage("$usage", 1) unless ($min_w <= $max_w);
# reverse complement and protein don't mix
&print_usage("$usage", 1) if (defined($revcomp) && $ALPHABET eq $PROT_ALPHABET);
# fixed_start only makes sense for triples
&print_usage("$usage", 1) if ($fixed_start && !$triples);

# defaults for scaling: if one not set but other is, set by
# scale_max = 1 - scale_min, otherwise default to [0.1..0.9]
if (!$scale_min && $scale_max) {
  $scale_min = 1 - $scale_max;
}
if (!$scale_max && $scale_min) {
  $scale_max = 1 - $scale_min;
}
$scale_min = 0.1 unless $scale_min;
$scale_max = 0.9 unless $scale_max;
die "invalid scale min,max: ($scale_min,$scale_max)"
    unless $scale_min > 0 && $scale_max < 1 && $scale_min < $scale_max;

die "only class and D priors implemented for protein" #FIXME -- protein not handled for all types yet
    if ($ALPHABET eq $PROT_ALPHABET && defined ($d_hyper_prior) ||
    defined ($hypergeomtric_prior));

$alpha_choice = "DNA" unless (defined ($alpha_choice));
my %alternates;
my $ambigs;
if (defined ($NorX)) { 
    $NorX .= $alpha_choice eq "DNA"?"N":"BJZX";
    $ambigs = '['.$NorX.']';
} else {
    $NorX = "";
}
# FIXME: ############# NOT IMPLEMENTED YET
# FIXME: should also do alternates for DNA
# FIXME: actual matches score 0 for any word or triple containing an ambig
# if you find one of the keys, check also if the stored value matches
if ($alpha_choice eq "PROT") {
    $alternates{"B"} = "DN";
    $alternates{"J"} = "IL";
    $alternates{"Z"} = "QE";
}

# the alphabet set from the command line or defaults to DNA
my @LETTERS = split(//, $ALPHABET);  # array of letters
my $ALENGTH = length($ALPHABET);  # length of alphabet

my @pos_sequences;
my @neg_sequences;
my ($L_pos, $L_neg);

my @translations;
if ($equivs) {
    my %equivsets;
    my @equivalents = split(/-/, $equivs);
    # check no duplicates within or accross equivalent sets
    foreach $_ (@equivalents) {
	my @letters = split(//);
	die "equivalent set `@letters' size < 2" if @letters < 2;
	foreach my $letter (@letters) {
	    die "`$letter' duplicated" if exists $equivsets{$letter};
	    $equivsets{$letter} = 1;
	}
	my $trtarget = pop(@letters);
	my $trsource = join("",@letters);
	
	push (@translations, [$trsource,$trtarget]);
    }
}

# when file is read, if equivalent sets specified, the returned sequences are translated to
# render all the alternatives as one of the equivalent letters
# read positive file storing names and FASTA comments and returning total sequence length
&read_seq_file ($pos_filename, \@pos_sequences, \$L_pos, $ALPHABET.$NorX, \@translations, $verbose);
# read negative file the same way
&read_seq_file ($neg_filename, \@neg_sequences, \$L_neg, $ALPHABET.$NorX, \@translations, $verbose);

my $N_pos  = @pos_sequences;
my $N_neg  = @neg_sequences;
my $L_mean = ($L_neg + $L_pos) / ($N_neg + $N_pos);

# count of number of positive, negative sequences that contain each w-mer
# indexed by w-mer, each entry is {"score", "last"} where "last" is the
# last position a w-mer was seen, to avoid counting twice for same sequence
# results returned in dictionaries

my $bestscores;

my ($last_min_score, $last_max_score, $max_score_w);

# need this to calculate number of sequences containing a w-mer or not,
# normalized for weighted count; if all sequences of same length
# each adds up to N_neg, N_pos
# absolute is an option for C priors, not settable from the command
# line in this version
my ($N_neg_normalized, $N_pos_normalized) =
    defined ($absolute) || defined ($d_prior) ?
        ($N_neg, $N_pos) :
        &normalize_size(\@neg_sequences, \@pos_sequences, $L_mean);

# fixed endpoint triples must score every width from 3 up
my $start = ($triples)?$min_triple:$min_w;
my $lastscores;    # remember scores for previous w

#foreach (1,2) {
#    (%pos_dictionary, %neg_dictionary) = ((), ());  # FIXME reset dictionaries before new prior type
for (my $width = $start; $width <= $max_w; $width++) {
    # new dictionaries each width for scoring: protein implicitly uses the old
    # values by reusing previous scores
    my (%pos_dictionary, %neg_dictionary);
    my $allscores;
    my ($min_score, $max_score) = (1, 0);
    if ($d_prior || $d_hyper_prior) {
        &D_prior(\@pos_sequences, \@neg_sequences, $revcomp, $width,
            \%pos_dictionary, \%neg_dictionary, $triples, \%alternates);
    } elsif ($hypergeomtric_prior) {    # don't normalize for seq length
        &class_prior_norm(\@pos_sequences, \@neg_sequences, $revcomp, $width, undef,
            \%pos_dictionary, \%neg_dictionary);
    } elsif ($hamming_prior) {
        $allscores = &hamming_prior(\@pos_sequences, \@neg_sequences, $revcomp, $width, $L_mean,
            \%pos_dictionary, \%neg_dictionary, $triples, \%alternates, \$min_score, \$max_score);
    } else {
        &class_prior_norm(\@pos_sequences, \@neg_sequences, $revcomp, $width, $L_mean,
            \%pos_dictionary, \%neg_dictionary, $triples, \%alternates);
    }
#print STDERR "min = $min_score max = $max_score\n";
    # do the classifier scores using the previously computed counts. Also O(LN)
    # normalized Ns = unnormalized for D-prior
    my ($N_X, $N_Y, $p_type);
    if (defined($priortype) && $priortype eq "D-hyper") {
        $N_X = $L_pos - $N_pos * ($width - 1);
        $N_Y = $L_neg - $N_neg * ($width - 1);
        $p_type = "hyper";
    } else {
        $N_X = $N_pos;
        $N_Y = $N_neg;
        $p_type = $priortype;
    }
    # hamming prior does its own thing
    $allscores = &normalized_classifier_scores (\@pos_sequences, \%pos_dictionary,
            \%neg_dictionary, $width, $N_X, $N_Y, $N_pos_normalized, $N_neg_normalized,
            \$min_score, \$max_score, $scale, $absolute, $p_type, $triples, $fixed_start,
            $ambigs, $lastscores)
    unless ($hamming_prior);
    # remember the last scores for the next iteration
    $lastscores = $allscores;
    # min_score, max_score now known for this w
    # record the maximum sequence set and its w
    # use >= so we get the widest motif achieving the max score
    # exception: for triples, this will always give us maxw so stop at first
    # width that finds the max score
    # if $maxrange defined (or != 0) then use $max_score-$min_score to choose best
    # w, otherwise use $maxscore
    if ($width >= $min_w) { # for triples, we may start at < minw
        my $newmax;
        if (!$last_max_score) {
            $newmax = 1;
        } elsif ($maxrange) {
            $newmax = $max_score?($max_score-$min_score) - ($last_max_score - $last_min_score):0;
        } else {
            $newmax = $max_score - $last_max_score;
        }
        # if not using triples, update scores if <= last maximum
        # for triples only update scores of beat last maximum
        if ((!$triples && $newmax >= 0) || ($newmax > 0)) {
            $bestscores = $allscores;
            $last_max_score = $max_score;
            $last_min_score = $min_score;
            $max_score_w = $width;
        }
        if (defined ($report)) {
            print STDERR $pos_filename."\t".$neg_filename."\t".$min_score."\t".$max_score."\t".
                $width."\t".$max_score_w."\n";
        }
    }
#    print STDERR "pos dictionary++++++++++: ".Dumper(\%pos_dictionary);
#    print STDERR "neg dictionary----------: ".Dumper(\%neg_dictionary);
}
#$d_prior = 1; # FIXME -- crude hack to force doing D prior second time
#$priortype = "D";
#}

#extend_w_mer ("", \%pos_dictionary);
#extend_w_mer ("", \%neg_dictionary);


#print STDERR "max score at width $max_score_w\n";
die "error finding best W" unless defined ($bestscores);
#print STDERR "best scores: ".Dumper($bestscores);

#rescale the scores linearly to standard range if required
if (defined ($scale)) {
    &scale_scores ($bestscores, $scale_min, $scale_max, $last_min_score, $last_max_score, $verbose);
}
# FIXME: model should be a variable; FIXME epsilon should be variable
if (defined($raw)) {
    &print_scores($bestscores, $scale_min, $scale_max, $last_min_score,
        $last_max_score, $max_score_w, $scale);
} else {
    &printScoreAsPrior ($bestscores, 'zoops', $max_score_w, $scale_min, $scale_max);
}

#&print_scores($allscores, $scale, $scale_min, $scale_max, $min_score,
#    $max_score, $min_w, $verbose);


################################################################################
#                       Subroutines                                            #
################################################################################


################################################################################
#
#       print_usage
#
#  Print the usage message and exit.
#
################################################################################
sub print_usage {
  my($usage, $status) = @_;
 
  if (-c STDOUT) {      # standard output is a terminal
    open(C, "| more");
    print C $usage;
    close C;
  } else {        # standard output not a terminal
    print STDERR $usage;
  }

  exit $status;
}
 
################################################################################
#       cleanup
#
#       cleanup stuff
#
################################################################################
sub cleanup {
  my($status, $msg) = @_;
  if ($status && "$msg") {print STDERR "$msg: $status\n";}
  exit($status);
}
