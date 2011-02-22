#!/usr/bin/perl
# AUTHOR: Timothy L. Bailey
# CREATE DATE: 4/16/2003
# MODIFIED: Philip Machanick
# MODIFIED DATE: 8/9/2009

#
# combinations of n choose m
# instead of efficient technique to calculate once-off,
# store factorials for reuse and use basic factorial formula
# n! / (m! x (n-m)!)

# factorial (0) = 1
my @__log_fact__ = (0);
sub comb_dp {
  my($n,                # n things
     $m                 # taken m at a time
  ) = @_;
  my($i, $t);
  $t = 0;
#  for($i=0; $i<=$m-1; $i++) { $t += log($n-$i) - log($i+1); }
  return(int(exp(&log_factorial($n) - &log_factorial($m) - 
		 &log_factorial($n-$m)) +0.5) + 0.0);
} # comb
sub log_comb_dp {
  my($n,                # n things
     $m                 # taken m at a time
  ) = @_;
# efficient formula for doing once: 
#  for($i=0; $i<=$m-1; $i++) { $t += log($n-$i) - log($i+1); }
# instead we use factorials because we can reuse values this way
  return &log_factorial($n) - &log_factorial($m) - &log_factorial($n-$m);
} # log_comb

# add any new values required to the saved factorial array
# and return the required log factorial
# uses global @__log_fact__ array to work from previous
# known values
sub log_factorial {
    my ($n  # log(factorial(n))
	) = @_;
    for (my $i = @__log_fact__; $i <= $n; $i++) {
	push (@__log_fact__, log($i) + $__log_fact__[-1]);
    }
    return $__log_fact__[$n];
}

1;
