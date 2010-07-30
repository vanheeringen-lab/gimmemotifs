use Test::More tests => 6;

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
my $weight1 = [ 1,1,1,1,1 ];
my $data1   = [
    [ 1.1, 2.2, 3.3, 4.4, 5.5, ], 
    [ 3.1, 3.2, 1.3, 2.4, 1.5, ], 
    [ 4.1, 2.2, 0.3, 5.4, 0.5, ], 
    [ 12.1, 2.0, 0.0, 5.0, 0.0, ], 
];
my $mask1 = [
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
];

#----------
# dataset 2
#
my $weight2 = [ 1,1 ];
my $data2   = [
    [ 1.1, 1.2 ],
    [ 1.4, 1.3 ],
    [ 1.1, 1.5 ],
    [ 2.0, 1.5 ],
    [ 1.7, 1.9 ],
    [ 1.7, 1.9 ],
    [ 5.7, 5.9 ],
    [ 5.7, 5.9 ],
    [ 3.1, 3.3 ],
    [ 5.4, 5.3 ],
    [ 5.1, 5.5 ],
    [ 5.0, 5.5 ],
    [ 5.1, 5.2 ],
];
my $mask2 = [
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
    [ 1, 1 ],
];


#------------------------------------------------------
# Tests
# 
my ($clusterid);

my %params;

%params = (
    transpose =>         0,
    dist      =>       'e',
    data      =>    $data1,
    mask      =>    $mask1,
    weight    =>  $weight1,
    niter     =>       100,
    nxnodes   =>        10,
    nynodes   =>        10,
);

$clusterid = Algorithm::Cluster::somcluster(%params);

is (scalar(@$data1), scalar(@$clusterid) );
is (scalar(@{$clusterid->[0]}), 2 );

%params = (
    transpose =>         0,
    dist      =>       'e',
    data      =>    $data2,
    mask      =>    $mask2,
    weight    =>  $weight2,
    niter     =>       100,
    nxnodes   =>        10,
    nynodes   =>        10,
);

$clusterid = Algorithm::Cluster::somcluster(%params);

is (scalar(@$data2), scalar(@$clusterid) );
is (scalar(@{$clusterid->[0]}), 2 );
