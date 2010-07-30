use Test::More tests => 28;

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
my $weight1 =  [ 1,1,1,1,1 ];

my $data1   =  [
    [ 1.1, 2.2, 3.3, 4.4, 5.5, ], 
    [ 3.1, 3.2, 1.3, 2.4, 1.5, ], 
    [ 4.1, 2.2, 0.3, 5.4, 0.5, ], 
    [ 12.1, 2.0, 0.0, 5.0, 0.0, ], 
];

my $mask1 =  [
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
];

#----------
# dataset 2
#
my $weight2 =  [ 1,1 ];

my $data2   =  [
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

my $mask2 =  [
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
my ($clusters, $centroids, $error, $found);
my ($i,$j);

my %params = (
    nclusters =>         3,
    transpose =>         0,
    method    =>       'a',
    dist      =>       'e',
);

#----------
# test dataset 1
#
($clusters, $error, $found) = Algorithm::Cluster::kcluster(
    %params,
    data      =>    $data1,
    mask      =>    $mask1,
    weight    =>  $weight1,
    npass     =>       100,
);

#----------
# Make sure that the length of @clusters matches the length of @data
ok( scalar @$data1 == scalar @$clusters);

#----------
# Test the cluster coordinates
ok ( $clusters->[ 0] != $clusters->[ 1] );
ok ( $clusters->[ 1] == $clusters->[ 2] );
ok ( $clusters->[ 2] != $clusters->[ 3] );

# Test the within-cluster sum of errors
ok( sprintf ("%7.3f", $error) == '  1.300');


#----------
# test dataset 2
#
$i=0;$j=0;
($clusters, $error, $found) = Algorithm::Cluster::kcluster(
    %params,
    data      =>    $data2,
    mask      =>    $mask2,
    weight    =>  $weight2,
    npass     =>       100,
);


#----------
# Make sure that the length of @clusters matches the length of @data
ok (scalar @$data2 == scalar @$clusters);

#----------
# Test the cluster coordinates
ok ($clusters->[ 0] == $clusters->[ 3]);
ok ($clusters->[ 0] != $clusters->[ 6]);
ok ($clusters->[ 0] != $clusters->[ 9]);
ok ($clusters->[11] == $clusters->[12]);

# Test the within-cluster sum of errors
ok ( sprintf ("%7.3f", $error) == '  1.012');

#----------
# test kcluster with initial cluster assignments
#
$initialid = [0,1,2,0,1,2,0,1,2,0,1,2,0];

($clusters, $error, $found) = Algorithm::Cluster::kcluster(
    %params,
    data      =>     $data2,
    mask      =>     $mask2,
    weight    =>   $weight2,
    npass     =>          1, 
    initialid => $initialid, 
);

#----------
# Test the cluster coordinates
ok ( $clusters->[ 0] == 2 );
ok ( $clusters->[ 1] == 2 );
ok ( $clusters->[ 2] == 2 );
ok ( $clusters->[ 3] == 2 );
ok ( $clusters->[ 4] == 2 );
ok ( $clusters->[ 5] == 2 );
ok ( $clusters->[ 6] == 0 );
ok ( $clusters->[ 7] == 0 );
ok ( $clusters->[ 8] == 2 );
ok ( $clusters->[ 9] == 1 );
ok ( $clusters->[10] == 1 );
ok ( $clusters->[11] == 1 );
ok ( $clusters->[12] == 1 );

# Test the within-cluster sum of errors
ok ( sprintf ("%7.3f", $error) == '  3.036' );
   
ok ($found == 1 );
                 
__END__
