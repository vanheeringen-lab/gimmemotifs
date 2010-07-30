use Test::More tests => 20;

use lib '../blib/lib','../blib/arch';

use_ok ("Algorithm::Cluster");
require_ok ("Algorithm::Cluster");


#########################


#------------------------------------------------------
# Data for Tests
# 

#----------
# dataset
#
my $gweight = [ 1,1,1,1,1 ];
my $eweight = [ 1,1,1,1 ];
my $data   = [
    [ 1.1, 2.2, 3.3, 4.4, 5.5, ], 
    [ 3.1, 3.2, 1.3, 2.4, 1.5, ], 
    [ 4.1, 2.2, 0.3, 5.4, 0.5, ], 
    [ 12.1, 2.0, 0.0, 5.0, 0.0, ], 
];
my $mask = [
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
];

#------------------------------------------------------
# Tests
# 
my $matrix;

#----------
# test dataset with transpose==0
#

$matrix = Algorithm::Cluster::distancematrix(
    transpose =>        0,
    dist      =>      'e',
    data      =>    $data,
    mask      =>    $mask,
    weight    => $gweight,
);


#----------
# Make sure that the length of $matrix matches the length of @data1
is (scalar @$data, scalar @$matrix);

#----------
# Test the values in the distance matrix

is (sprintf ("%7.3f", $matrix->[1]->[0] ), '  5.800');
is (sprintf ("%7.3f", $matrix->[2]->[0] ), '  8.800');
is (sprintf ("%7.3f", $matrix->[2]->[1] ), '  2.600');
is (sprintf ("%7.3f", $matrix->[3]->[0] ), ' 32.508');
is (sprintf ("%7.3f", $matrix->[3]->[1] ), ' 18.628');
is (sprintf ("%7.3f", $matrix->[3]->[2] ), ' 12.908');


#----------
# test dataset with transpose==1
#

$matrix = Algorithm::Cluster::distancematrix(
    transpose =>        1,
    dist      =>      'e',
    data      =>    $data,
    mask      =>    $mask,
    weight    => $eweight,
);

#----------
# Make sure that the length of $matrix matches the length of @data1
is (scalar @{$data->[0]}, scalar @$matrix );

#----------
# Test the values in the distance matrix


is (sprintf ("%6.2f", $matrix->[1]->[0] ), ' 26.71');
is (sprintf ("%6.2f", $matrix->[2]->[0] ), ' 42.23');
is (sprintf ("%6.2f", $matrix->[2]->[1] ), '  3.11');
is (sprintf ("%6.2f", $matrix->[3]->[0] ), ' 15.87');
is (sprintf ("%6.2f", $matrix->[3]->[1] ), '  6.18');
is (sprintf ("%6.2f", $matrix->[3]->[2] ), ' 13.36');
is (sprintf ("%6.2f", $matrix->[4]->[0] ), ' 45.32');
is (sprintf ("%6.2f", $matrix->[4]->[1] ), '  5.17');
is (sprintf ("%6.2f", $matrix->[4]->[2] ), '  1.23');
is (sprintf ("%6.2f", $matrix->[4]->[3] ), ' 12.76');
