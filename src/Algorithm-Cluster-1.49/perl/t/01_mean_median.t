use Test::More tests => 6;

use lib '../blib/lib','../blib/arch';

use_ok ("Algorithm::Cluster");
require_ok ("Algorithm::Cluster");


#########################


#------------------------------------------------------
# Tests
# 
my ($meanval, $medianval);
my $dummy;
my $dummy_ref = \$dummy;
my $dummy_sub_ref = sub {};

my $data1 = [ 34.3, 3, 2 ];
my $data2 = [ 5, 10 ,15, 20 ];


is (sprintf ("%7.4f", Algorithm::Cluster::mean($data1)), '13.1000');
is (sprintf ("%7.4f", Algorithm::Cluster::mean($data2)), '12.5000');

is (sprintf ("%7.4f", Algorithm::Cluster::median($data1)), ' 3.0000');
is (sprintf ("%7.4f", Algorithm::Cluster::median($data2)), '12.5000');
