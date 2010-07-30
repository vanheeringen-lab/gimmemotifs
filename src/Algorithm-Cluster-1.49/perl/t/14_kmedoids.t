use Test::More tests => 30;

use lib '../blib/lib','../blib/arch';

use_ok ("Algorithm::Cluster");
require_ok ("Algorithm::Cluster");


#########################


#------------------------------------------------------
# Data for Tests
# 

#----------
# dataset 1
#
my $matrix   =  [
    [],
    [ 3.4],
    [ 4.3, 10.1],
    [ 3.7, 11.5,  1.1],
    [ 1.7,  4.1,  3.4,  3.4],
    [10.1, 20.5,  2.5,  2.7,  9.8],
    [ 2.5,  3.7,  3.1,  3.6,  1.1, 10.1],
    [ 3.4,  2.2,  8.8,  8.7,  3.3, 16.6,  2.7],
    [ 2.1,  7.7,  2.7,  1.9,  1.8,  5.7,  3.4,  5.2],
    [ 1.6,  1.8,  9.2,  8.7,  3.4, 16.8,  4.2,  1.3,  5.0],
    [ 2.7,  3.7,  5.5,  5.5,  1.9, 11.5,  2.0,  1.7,  2.1,  3.1],
    [10.0, 19.3,  1.0,  3.7,  9.1,  1.2,  9.3, 15.7,  6.3, 16.0, 11.5]
];

#------------------------------------------------------
# Tests
# 

my ($clusters, $error, $found);

#------------------------------------------------------
# Test with repeated runs of the k-medoids algorithm
# 

my %params1 = (
        nclusters =>         4,
        distances =>   $matrix,
        npass     =>       100,
);
                                                                                
($clusters, $error, $found) = Algorithm::Cluster::kmedoids(%params1);

#----------
# Make sure that the length of @clusters matches the length of @data
is (scalar @$matrix, scalar @$clusters );

#----------
# Test the cluster assignments
is ($clusters->[ 0], 9);
is ($clusters->[ 1], 9);
is ($clusters->[ 2], 2);
is ($clusters->[ 3], 2);
is ($clusters->[ 4], 4);
is ($clusters->[ 5], 5);
is ($clusters->[ 6], 4);
is ($clusters->[ 7], 9);
is ($clusters->[ 8], 4);
is ($clusters->[ 9], 9);
is ($clusters->[10], 4);
is ($clusters->[11], 2);

# Test the within-cluster sum of errors
is (sprintf ("%7.3f", $error), ' 11.600');


#------------------------------------------------------
# Test the k-medoids algorithm with a specified initial clustering
# 

$initialid = [0,0,1,1,1,2,2,2,3,3,3,3];

my %params2 = (
    nclusters =>         4,
    distances =>   $matrix,
    npass     =>         1,
    initialid => $initialid,
);
                                                                                
($clusters, $error, $found) = Algorithm::Cluster::kmedoids(%params2);

#----------
# Make sure that the length of @clusters matches the length of @data
is (scalar @$matrix, scalar @$clusters );

#----------
# Test the cluster assignments
is ($clusters->[ 0], 9);
is ($clusters->[ 1], 9);
is ($clusters->[ 2], 2);
is ($clusters->[ 3], 2);
is ($clusters->[ 4], 4);
is ($clusters->[ 5], 2);
is ($clusters->[ 6], 6);
is ($clusters->[ 7], 9);
is ($clusters->[ 8], 4);
is ($clusters->[ 9], 9);
is ($clusters->[10], 4);
is ($clusters->[11], 2);

# Test the within-cluster sum of errors
is (sprintf ("%7.3f", $error), " 13.000");
