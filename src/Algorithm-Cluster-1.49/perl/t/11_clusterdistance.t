use Test::More tests => 8;

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
my $data1 = [
    [  1.1, 2.2, 3.3, 4.4, 5.5, ], 
    [  3.1, 3.2, 1.3, 2.4, 1.5, ], 
    [  4.1, 2.2, 0.3, 5.4, 0.5, ], 
    [ 12.1, 2.0, 0.0, 5.0, 0.0, ], 
];
my $mask1 = [
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
    [ 1, 1, 1, 1, 1, ], 
];
my $data1_c1 = [ 0 ];
my $data1_c2 = [ 1,2 ];
my $data1_c3 = [ 3 ];


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
my $data2_c1 = [ 0, 1, 2, 3 ];
my $data2_c2 = [ 4, 5, 6, 7 ];
my $data2_c3 = [ 8 ];


#------------------------------------------------------
# Tests
# 
my $distance;

#----------
# test dataset 1
#
my %params = (
    transpose =>         0,
    method    =>       'a',
    dist      =>       'e',
    data      =>    $data1,
    mask      =>    $mask1,
    weight    =>  $weight1,
    cluster1  => $data1_c1,
    cluster2  => $data1_c2,
);

$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), '  6.650');

$params{cluster1} = $data1_c1;
$params{cluster2} = $data1_c3;
$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), ' 32.508');

$params{cluster1} = $data1_c2;
$params{cluster2} = $data1_c3;
$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), ' 15.118');

#----------
# test dataset 2
#
%params = (
    transpose =>         0,
    method    =>       'a',
    dist      =>       'e',
    data      =>    $data2,
    mask      =>    $mask2,
    weight    =>  $weight2,
    cluster1  => $data2_c1,
    cluster2  => $data2_c2,
);
$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), '  5.833');

$params{cluster1} = $data2_c1;
$params{cluster2} = $data2_c3;
$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), '  3.298');

$params{cluster1} = $data2_c2;
$params{cluster2} = $data2_c3;
$distance = Algorithm::Cluster::clusterdistance(%params);
is (sprintf ("%7.3f", $distance ), '  0.360');
