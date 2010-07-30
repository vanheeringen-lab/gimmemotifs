use Test::More tests => 74;

use lib '../blib/lib','../blib/arch';

use_ok ("Algorithm::Cluster");
require_ok ("Algorithm::Cluster");

use warnings;


#########################


#------------------------------------------------------
# Tests
# 
my ($mean, $coordinates, $pc, $eigenvalues);

my $data1   =  [
    [ 3.1, 1.2 ],
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

my $data2 = [
    [ 2.3, 4.5, 1.2, 6.7, 5.3, 7.1],
    [ 1.3, 6.5, 2.2, 5.7, 6.2, 9.1],
    [ 3.2, 7.2, 3.2, 7.4, 7.3, 8.9],
    [ 4.2, 5.2, 9.2, 4.4, 6.3, 7.2],
];


($mean, $coordinates, $pc, $eigenvalues) =  Algorithm::Cluster::pca($data1);

is (sprintf("%7.4f", $mean->[0]), " 3.5462");
is (sprintf("%7.4f", $mean->[1]), " 3.5308");
is (sprintf("%7.4f", $coordinates->[0][0]), " 2.0323");
is (sprintf("%7.4f", $coordinates->[0][1]), " 1.2252");
is (sprintf("%7.4f", $coordinates->[1][0]), " 3.0937");
is (sprintf("%7.4f", $coordinates->[1][1]), "-0.1065");
is (sprintf("%7.4f", $coordinates->[2][0]), " 3.1453");
is (sprintf("%7.4f", $coordinates->[2][1]), "-0.4633");
is (sprintf("%7.4f", $coordinates->[3][0]), " 2.5440");
is (sprintf("%7.4f", $coordinates->[3][1]), " 0.2063");
is (sprintf("%7.4f", $coordinates->[4][0]), " 2.4468");
is (sprintf("%7.4f", $coordinates->[4][1]), "-0.2841");
is (sprintf("%7.4f", $coordinates->[5][0]), " 2.4468");
is (sprintf("%7.4f", $coordinates->[5][1]), "-0.2841");
is (sprintf("%7.4f", $coordinates->[6][0]), "-3.2019");
is (sprintf("%7.4f", $coordinates->[6][1]), " 0.0197");
is (sprintf("%7.4f", $coordinates->[7][0]), "-3.2019");
is (sprintf("%7.4f", $coordinates->[7][1]), " 0.0197");
is (sprintf("%7.4f", $coordinates->[8][0]), " 0.4698");
is (sprintf("%7.4f", $coordinates->[8][1]), "-0.1778");
is (sprintf("%7.4f", $coordinates->[9][0]), "-2.5550");
is (sprintf("%7.4f", $coordinates->[9][1]), " 0.1973");
is (sprintf("%7.4f", $coordinates->[10][0]), "-2.5034");
is (sprintf("%7.4f", $coordinates->[10][1]), "-0.1595");
is (sprintf("%7.4f", $coordinates->[11][0]), "-2.4366");
is (sprintf("%7.4f", $coordinates->[11][1]), "-0.2339");
is (sprintf("%7.4f", $coordinates->[12][0]), "-2.2802");
is (sprintf("%7.4f", $coordinates->[12][1]), " 0.0409");
is (sprintf("%7.4f", $pc->[0][0]), "-0.6681");
is (sprintf("%7.4f", $pc->[0][1]), "-0.7441");
is (sprintf("%7.4f", $pc->[1][0]), " 0.7441");
is (sprintf("%7.4f", $pc->[1][1]), "-0.6681");
is (sprintf("%7.4f", $eigenvalues->[0]), " 9.3110");
is (sprintf("%7.4f", $eigenvalues->[1]), " 1.4437");


($mean, $coordinates, $pc, $eigenvalues) =  Algorithm::Cluster::pca($data2);

is (sprintf("%7.4f", $mean->[0]), " 2.7500");
is (sprintf("%7.4f", $mean->[1]), " 5.8500");
is (sprintf("%7.4f", $mean->[2]), " 3.9500");
is (sprintf("%7.4f", $mean->[3]), " 6.0500");
is (sprintf("%7.4f", $mean->[4]), " 6.2750");
is (sprintf("%7.4f", $mean->[5]), " 8.0750");
is (sprintf("%7.4f", $coordinates->[0][0]), " 2.6461");
is (sprintf("%7.4f", $coordinates->[0][1]), "-2.1422");
is (sprintf("%7.4f", $coordinates->[0][2]), "-0.5662");
is (sprintf("%7.4f", abs($coordinates->[0][3])), " 0.0000");
is (sprintf("%7.4f", $coordinates->[1][0]), " 2.0644");
is (sprintf("%7.4f", $coordinates->[1][1]), " 0.5554");
is (sprintf("%7.4f", $coordinates->[1][2]), " 1.4819");
is (sprintf("%7.4f", abs($coordinates->[1][3])), " 0.0000");
is (sprintf("%7.4f", $coordinates->[2][0]), " 1.0687");
is (sprintf("%7.4f", $coordinates->[2][1]), " 1.9994");
is (sprintf("%7.4f", $coordinates->[2][2]), "-1.0007");
is (sprintf("%7.4f", abs($coordinates->[2][3])), " 0.0000");
is (sprintf("%7.4f", $coordinates->[3][0]), "-5.7792");
is (sprintf("%7.4f", $coordinates->[3][1]), "-0.4127");
is (sprintf("%7.4f", $coordinates->[3][2]), " 0.0851");
is (sprintf("%7.4f", abs($coordinates->[3][3])), " 0.0000");
is (sprintf("%7.4f", $pc->[0][0]), "-0.2638");
is (sprintf("%7.4f", $pc->[0][1]), " 0.0648");
is (sprintf("%7.4f", $pc->[0][2]), "-0.9176");
is (sprintf("%7.4f", $pc->[0][3]), " 0.2615");
is (sprintf("%7.4f", $pc->[1][0]), " 0.0507");
is (sprintf("%7.4f", $pc->[1][1]), " 0.6862");
is (sprintf("%7.4f", $pc->[1][2]), " 0.1382");
is (sprintf("%7.4f", $pc->[1][3]), " 0.1978");
is (sprintf("%7.4f", $pc->[2][0]), "-0.6300");
is (sprintf("%7.4f", $pc->[2][1]), " 0.0912");
is (sprintf("%7.4f", $pc->[2][2]), " 0.0456");
is (sprintf("%7.4f", $pc->[2][3]), "-0.6746");
# The last eigenvalue is zero. The corresponding eigenvector is strongly
# affected by roundoff error, so we don't test it. For PCA, it doesn't matter
# since the coordinates along this eigenvector are zero anyway.
is (sprintf("%7.4f", $eigenvalues->[0]), " 6.7679");
is (sprintf("%7.4f", $eigenvalues->[1]), " 3.0109");
is (sprintf("%7.4f", $eigenvalues->[2]), " 1.8776");
is (sprintf("%7.4f", abs($eigenvalues->[3])), " 0.0000");

