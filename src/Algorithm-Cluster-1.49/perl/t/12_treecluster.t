use Test::More tests => 224;

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
    [ 0.8223, 0.9295 ],
    [ 1.4365, 1.3223 ],
    [ 1.1623, 1.5364 ],
    [ 2.1826, 1.1934 ],
    [ 1.7763, 1.9352 ],
    [ 1.7215, 1.9912 ],
    [ 2.1812, 5.9935 ],
    [ 5.3290, 5.9452 ],
    [ 3.1491, 3.3454 ],
    [ 5.1923, 5.3156 ],
    [ 4.7735, 5.4012 ],
    [ 5.1297, 5.5645 ],
    [ 5.3934, 5.1823 ],
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
my $tree;
my $node;

#----------
# test dataset 1
#

#--------------[PALcluster]-------
my %params = (
    transpose  =>         0,
    method     =>       'a',
    dist       =>       'e',
    data       =>    $data1,
    mask       =>    $mask1,
    weight     =>  $weight1,
);

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data1) - 1, $tree->length );

$node = $tree->get(0);
is ($node->left, 2);
is ($node->right, 1);
is (sprintf("%7.3f", $node->distance), "  2.600");

$node = $tree->get(1);
is ($node->left, -1);
is ($node->right, 0);
is (sprintf("%7.3f", $node->distance), "  7.300");

$node = $tree->get(2);
is ($node->left, 3);
is ($node->right, -2);
is (sprintf("%7.3f", $node->distance), " 21.348");


#--------------[PSLcluster]-------
$params{method} = 's';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data1) - 1, $tree->length );

$node = $tree->get(0);
is ($node->left, 1);
is ($node->right, 2);
is (sprintf("%7.3f", $node->distance), "  2.600");

$node = $tree->get(1);
is ($node->left, 0);
is ($node->right, -1);
is (sprintf("%7.3f", $node->distance), "  5.800");

$node = $tree->get(2);
is ($node->left, -2);
is ($node->right, 3);
is (sprintf("%7.3f", $node->distance), " 12.908");


#--------------[PCLcluster]-------
$params{method} = 'c';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data1) - 1, $tree->length );

$node = $tree->get(0);
is ($node->left, 1);
is ($node->right, 2);
is (sprintf("%7.3f", $node->distance), "  2.600");

$node = $tree->get(1);
is ($node->left, 0);
is ($node->right, -1);
is (sprintf("%7.3f", $node->distance), "  6.650");

$node = $tree->get(2);
is ($node->left, -2);
is ($node->right, 3);
is (sprintf("%7.3f", $node->distance), " 19.437");


#--------------[PMLcluster]-------
$params{method} = 'm';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data1) - 1, $tree->length );

$node = $tree->get(0);
is ($node->left, 2);
is ($node->right, 1);
is (sprintf("%7.3f", $node->distance), "  2.600");

$node = $tree->get(1);
is ($node->left, -1);
is ($node->right, 0);
is (sprintf("%7.3f", $node->distance), "  8.800");

$node = $tree->get(2);
is ($node->left, 3);
is ($node->right, -2);
is (sprintf("%7.3f", $node->distance), " 32.508");



#----------
# test dataset 2
#

#--------------[PALcluster]-------
%params = (
    transpose  =>         0,
    method     =>       'a',
    dist       =>       'e',
    data       =>    $data2,
    mask       =>    $mask2,
    weight     =>  $weight2,
);

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data2) - 1, $tree->length);


$node = $tree->get(0);
is ($node->left, 5);
is ($node->right, 4);
is (sprintf("%7.3f", $node->distance), "  0.003");

$node = $tree->get(1);
is ($node->left, 9);
is ($node->right, 12);
is (sprintf("%7.3f", $node->distance), "  0.029");

$node = $tree->get(2);
is ($node->left, 2);
is ($node->right, 1);
is (sprintf("%7.3f", $node->distance), "  0.061");

$node = $tree->get(3);
is ($node->left, 11);
is ($node->right, -2);
is (sprintf("%7.3f", $node->distance), "  0.070");

$node = $tree->get(4);
is ($node->left, -4);
is ($node->right, 10);
is (sprintf("%7.3f", $node->distance), "  0.128");

$node = $tree->get(5);
is ($node->left, 7);
is ($node->right, -5);
is (sprintf("%7.3f", $node->distance), "  0.224");

$node = $tree->get(6);
is ($node->left, -3);
is ($node->right, 0);
is (sprintf("%7.3f", $node->distance), "  0.254");

$node = $tree->get(7);
is ($node->left, -1);
is ($node->right, 3);
is (sprintf("%7.3f", $node->distance), "  0.391");

$node = $tree->get(8);
is ($node->left, -8);
is ($node->right, -7);
is (sprintf("%7.3f", $node->distance), "  0.532");

$node = $tree->get(9);
is ($node->left, 8);
is ($node->right, -9);
is (sprintf("%7.3f", $node->distance), "  3.234");

$node = $tree->get(10);
is ($node->left, -6);
is ($node->right, 6);
is (sprintf("%7.3f", $node->distance), "  4.636");

$node = $tree->get(11);
is ($node->left, -11);
is ($node->right, -10);
is (sprintf("%7.3f", $node->distance), " 12.741");



#--------------[PSLcluster]-------
$params{method} = 's';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data2) - 1, $tree->length );


$node = $tree->get(0);
is ($node->left, 4);
is ($node->right, 5);
is (sprintf("%7.3f", $node->distance), "  0.003");

$node = $tree->get(1);
is ($node->left, 9);
is ($node->right, 12);
is (sprintf("%7.3f", $node->distance), "  0.029");

$node = $tree->get(2);
is ($node->left, 11);
is ($node->right, -2);
is (sprintf("%7.3f", $node->distance), "  0.033");

$node = $tree->get(3);
is ($node->left, 1);
is ($node->right, 2);
is (sprintf("%7.3f", $node->distance), "  0.061");

$node = $tree->get(4);
is ($node->left, 10);
is ($node->right, -3);
is (sprintf("%7.3f", $node->distance), "  0.077");

$node = $tree->get(5);
is ($node->left, 7);
is ($node->right, -5);
is (sprintf("%7.3f", $node->distance), "  0.092");

$node = $tree->get(6);
is ($node->left, 0);
is ($node->right, -4);
is (sprintf("%7.3f", $node->distance), "  0.242");

$node = $tree->get(7);
is ($node->left, -7);
is ($node->right, -1);
is (sprintf("%7.3f", $node->distance), "  0.246");

$node = $tree->get(8);
is ($node->left, 3);
is ($node->right, -8);
is (sprintf("%7.3f", $node->distance), "  0.287");

$node = $tree->get(9);
is ($node->left, -9);
is ($node->right, 8);
is (sprintf("%7.3f", $node->distance), "  1.936");

$node = $tree->get(10);
is ($node->left, -10);
is ($node->right, -6);
is (sprintf("%7.3f", $node->distance), "  3.432");

$node = $tree->get(11);
is ($node->left, 6);
is ($node->right, -11);
is (sprintf("%7.3f", $node->distance), "  3.535");


#--------------[PCLcluster]-------
$params{method} = 'c';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is (scalar(@$data2) - 1, $tree->length );



$node = $tree->get(0);
is ($node->left, 4);
is ($node->right, 5);
is (sprintf("%7.3f", $node->distance), "  0.003");

$node = $tree->get(1);
is ($node->left, 12);
is ($node->right, 9);
is (sprintf("%7.3f", $node->distance), "  0.029");

$node = $tree->get(2);
is ($node->left, 1);
is ($node->right, 2);
is (sprintf("%7.3f", $node->distance), "  0.061");

$node = $tree->get(3);
is ($node->left, -2);
is ($node->right, 11);
is (sprintf("%7.3f", $node->distance), "  0.063");

$node = $tree->get(4);
is ($node->left, 10);
is ($node->right, -4);
is (sprintf("%7.3f", $node->distance), "  0.109");

$node = $tree->get(5);
is ($node->left, -5);
is ($node->right, 7);
is (sprintf("%7.3f", $node->distance), "  0.189");

$node = $tree->get(6);
is ($node->left, 0);
is ($node->right, -3);
is (sprintf("%7.3f", $node->distance), "  0.239");

$node = $tree->get(7);
is ($node->left, 3);
is ($node->right, -1);
is (sprintf("%7.3f", $node->distance), "  0.390");

$node = $tree->get(8);
is ($node->left, -7);
is ($node->right, -8);
is (sprintf("%7.3f", $node->distance), "  0.382");

$node = $tree->get(9);
is ($node->left, -9);
is ($node->right, 8);
is (sprintf("%7.3f", $node->distance), "  3.063");

$node = $tree->get(10);
is ($node->left, 6);
is ($node->right, -6);
is (sprintf("%7.3f", $node->distance), "  4.578");

$node = $tree->get(11);
is ($node->left, -10);
is ($node->right, -11);
is (sprintf("%7.3f", $node->distance), " 11.536");


#--------------[PMLcluster]-------
$params{method} = 'm';

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is ( scalar(@$data2) - 1, $tree->length );


$node = $tree->get(0);
is ($node->left, 5);
is ($node->right, 4);
is (sprintf("%7.3f", $node->distance), "  0.003");

$node = $tree->get(1);
is ($node->left, 9);
is ($node->right, 12);
is (sprintf("%7.3f", $node->distance), "  0.029");

$node = $tree->get(2);
is ($node->left, 2);
is ($node->right, 1);
is (sprintf("%7.3f", $node->distance), "  0.061");

$node = $tree->get(3);
is ($node->left, 11);
is ($node->right, 10);
is (sprintf("%7.3f", $node->distance), "  0.077");

$node = $tree->get(4);
is ($node->left, -2);
is ($node->right, -4);
is (sprintf("%7.3f", $node->distance), "  0.216");

$node = $tree->get(5);
is ($node->left, -3);
is ($node->right, 0);
is (sprintf("%7.3f", $node->distance), "  0.266");

$node = $tree->get(6);
is ($node->left, -5);
is ($node->right, 7);
is (sprintf("%7.3f", $node->distance), "  0.302");

$node = $tree->get(7);
is ($node->left, -1);
is ($node->right, 3);
is (sprintf("%7.3f", $node->distance), "  0.425");

$node = $tree->get(8);
is ($node->left, -8);
is ($node->right, -6);
is (sprintf("%7.3f", $node->distance), "  0.968");

$node = $tree->get(9);
is ($node->left, 8);
is ($node->right, 6);
is (sprintf("%7.3f", $node->distance), "  3.975");

$node = $tree->get(10);
is ($node->left, -10);
is ($node->right, -7);
is (sprintf("%7.3f", $node->distance), "  5.755");

$node = $tree->get(11);
is ($node->left, -11);
is ($node->right, -9);
is (sprintf("%7.3f", $node->distance), " 22.734");


#-------[treecluster on a distance matrix]------------

my $matrix   =  [
        [],
        [ 3.4],
        [ 4.3, 10.1],
        [ 3.7, 11.5,  1.0],
        [ 1.6,  4.1,  3.4,  3.4],
        [10.1, 20.5,  2.5,  2.7,  9.8],
        [ 2.5,  3.7,  3.1,  3.6,  1.1, 10.1],
        [ 3.4,  2.2,  8.8,  8.7,  3.3, 16.6,  2.7],
        [ 2.1,  7.7,  2.7,  1.9,  1.8,  5.7,  3.4,  5.2],
        [ 1.4,  1.7,  9.2,  8.7,  3.4, 16.8,  4.2,  1.3,  5.0],
        [ 2.7,  3.7,  5.5,  5.5,  1.9, 11.5,  2.0,  1.5,  2.1,  3.1],
        [10.0, 19.3,  2.2,  3.7,  9.1,  1.2,  9.3, 15.7,  6.3, 16.0, 11.5]
];

%params = (
    method     =>       's',
    data       =>   $matrix,
);

$tree = Algorithm::Cluster::treecluster(%params);

# Make sure that @clusters and @centroids are the right length
is ( scalar(@$matrix) - 1, $tree->length );

$node = $tree->get(0);
is ($node->left, 2);
is ($node->right, 3);
is (sprintf("%7.3f", $node->distance), "  1.000");

$node = $tree->get(1);
is ($node->left, 4);
is ($node->right, 6);
is (sprintf("%7.3f", $node->distance), "  1.100");

$node = $tree->get(2);
is ($node->left, 5);
is ($node->right, 11);
is (sprintf("%7.3f", $node->distance), "  1.200");

$node = $tree->get(3);
is ($node->left, 7);
is ($node->right, 9);
is (sprintf("%7.3f", $node->distance), "  1.300");

$node = $tree->get(4);
is ($node->left, 0);
is ($node->right, -4);
is (sprintf("%7.3f", $node->distance), "  1.400");

$node = $tree->get(5);
is ($node->left, -5);
is ($node->right, 10);
is (sprintf("%7.3f", $node->distance), "  1.500");

$node = $tree->get(6);
is ($node->left, -2);
is ($node->right, -6);
is (sprintf("%7.3f", $node->distance), "  1.600");

$node = $tree->get(7);
is ($node->left, 1);
is ($node->right, -7);
is (sprintf("%7.3f", $node->distance), "  1.700");

$node = $tree->get(8);
is ($node->left, 8);
is ($node->right, -8);
is (sprintf("%7.3f", $node->distance), "  1.800");

$node = $tree->get(9);
is ($node->left, -1);
is ($node->right, -9);
is (sprintf("%7.3f", $node->distance), "  1.900");

$node = $tree->get(10);
is ($node->left, -10);
is ($node->right, -3);
is (sprintf("%7.3f", $node->distance), "  2.200");
