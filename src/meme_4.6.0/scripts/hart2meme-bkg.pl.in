#!@WHICHPERL@
# AUTHOR: Timothy L. Bailey
# CREATE DATE: 14/11/2008
#use strict;
use File::Basename;

$PGM = $0;			# name of program
$PGM =~ s#.*/##;                # remove part up to last slash
#@args = @ARGV;			# arguments to program
$| = 1;				# flush after all prints
$SIG{'INT'} = \&cleanup;	# interrupt handler
# Note: so that interrupts work, always use for system calls:
# 	if ($status = system("$command")) {&cleanup($status)}

# requires
push(@INC, split(":", $ENV{'PATH'}));	# look in entire path

# defaults

my $usage = <<USAGE;		# usage message
  USAGE:
	$PGM
	
	Convert a Hartemink style conditional probability file into
	a MEME (fasta-get-markov) style joint probability file.

	Currently only works up to order 5.

	Reads standard input.
	Writes standard output.

        Copyright
        (2008) The University of Queensland
        All Rights Reserved.
        Author: Timothy L. Bailey
USAGE

$nargs = 0;			# number of required args
if ($#ARGV+1 < $nargs) { &print_usage("$usage", 1); }

# get input arguments
while ($#ARGV >= 0) {
  $_ = shift;
  if ($_ eq "-h") {				# help
    &print_usage("$usage", 0);
  } else {
    &print_usage("$usage", 1);
  }
}

# Get tuples for up to order 5
$tuple[0] = 'a';
$tuple[1] = 'c';
$tuple[2] = 'g';
$tuple[3] = 't';
$start = 0;
$end = 3;
$n = $end - $start + 1;
for ($order=1; $order<=5; $order++) {
  $next = $end + 1;
  for ($i=$start; $i<=$end; $i++) {
    for ($j=0; $j<=3; $j++) {
      $new_tuple = $tuple[$i] . $tuple[$j];
      $tuple[$next++] = $new_tuple;
    }
  }
  $n *= 4;
  $start = $end+1;
  $end = $start + $n - 1;
}

# read the Hartemink style conditional probabilities
$i = 0;
while (<STDIN>) {
  next if (/^#/ || /^\s*$/);                    # skip comment, blank lines
  $p{$tuple[$i++]} = $_;		# save conditional with tuple name
};
$n = $i;

# print the MEME style joint probabilities
for ($i=0; $i<$n; $i++) {
  $t = $tuple[$i];
  if ($i>3) {
    $y = substr $t, 0, -1;	# tuple = p(X|Y)
    $p{$t} *= $p{$y};		# p(X,Y) = P(X|Y) * P(Y)
  }
  printf "%s %.5f\n", $t, $p{$t};
}

$status = 1;

# cleanup files
&cleanup($status, "");
 
################################################################################
#                       Subroutines                                            #
################################################################################
 
################################################################################
#
#       print_usage
#
#	Print the usage message and exit.
#
################################################################################
sub print_usage {
  my($usage, $status) = @_;
 
  if (-c STDOUT) {			# standard output is a terminal
    open(C, "| more");
    print C $usage;
    close C;
  } else {				# standard output not a terminal
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
