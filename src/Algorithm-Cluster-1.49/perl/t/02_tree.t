use Test::More tests => 15;

use lib '../blib/lib','../blib/arch';

use_ok ("Algorithm::Cluster");
require_ok ("Algorithm::Cluster");


#########################

#------------------------------------------------------
# Tests
#

my $node;

my $node1 = Algorithm::Cluster::Node->new(1,2,3.1);
my $node2 = Algorithm::Cluster::Node->new(-1,3,5.3);
my $node3 = Algorithm::Cluster::Node->new(4,0,5.9);
my $node4 = Algorithm::Cluster::Node->new(-2,-3,7.8);
my @nodes = [$node1,$node2,$node3,$node4];

my $tree = Algorithm::Cluster::Tree->new(@nodes);
is ($tree->length, 4);

$node = $tree->get(0);
is ($node->left, 1);
is ($node->right, 2);
is (sprintf ("%7.4f", $node->distance), ' 3.1000');

$node = $tree->get(1);
is ($node->left, -1);
is ($node->right, 3);
is (sprintf ("%7.4f", $node->distance), ' 5.3000');

$node = $tree->get(2);
is ($node->left, 4);
is ($node->right, 0);
is (sprintf ("%7.4f", $node->distance), ' 5.9000');

$node = $tree->get(3);
is ($node->left, -2);
is ($node->right, -3);
is (sprintf ("%7.4f", $node->distance), ' 7.8000');
