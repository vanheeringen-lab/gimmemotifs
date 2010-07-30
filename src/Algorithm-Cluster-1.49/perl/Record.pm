package Algorithm::Cluster::Record;
use strict;

use Algorithm::Cluster;

sub new {
    my $class = shift;
    my $self = {};
    $self->{data} = undef;
    $self->{mask} = undef;
    $self->{geneid} = undef;
    $self->{genename} = undef;
    $self->{gweight} = undef;
    $self->{gorder} = undef;
    $self->{expid} = undef;
    $self->{eweight} = undef;
    $self->{eorder} = undef;
    $self->{uniqid} = undef;
    bless($self, $class);
    return $self;
}

sub read {
    my $self = shift;
    my $handle = shift;
    my $line = <$handle>;
    chomp($line);
    my @words = split(/\t/, $line);
    my $n = scalar @words;
    $self->{uniqid} = $words[0];
    $self->{expid} = [];
    my %cols = (0 => 'GENEID');
    my $i;
    for ($i = 1; $i < $n; $i++) {
        my $word = $words[$i];
        if ($word eq 'NAME') {
            $cols{$i} = $word;
            $self->{genename} = ();
        }
        elsif ($word eq 'GWEIGHT') {
            $cols{$i} = $word;
            $self->{gweight} = ();
        }
        elsif ($word eq 'GORDER') {
            $cols{$i} = $word;
            $self->{gorder} = ();
        }
        else {
            push(@{$self->{expid}}, $word);
        }
    }
    $self->{geneid} = [];
    $self->{data} = [];
    $self->{mask} = [];
    my $needmask = 0;
    while ($line = <$handle>) {
	my $count = ($line =~ tr/\t//);
        @words = split(/\t/, $line);
	chomp @words;
        scalar @words == $n or die "Line with " . scalar @words . " columns found (expected $n): $!";
        my $start = 0;
        for my $key (keys %cols) {
            if ($key > $start) {
                $start = $key;
            }
        }
        if ($words[0] eq 'EWEIGHT') {
            @{$self->{eweight}} = @words[$start+1..$n-1];
        }
        elsif ($words[0] eq 'EORDER') {
            @{$self->{eorder}} = @words[$start+1..$n-1];
        }
        else {
            my @rowdata = ();
            my @rowmask = ();
            for ($i = 0; $i < $n; $i++) {
                my $word = $words[$i];
                if (defined $cols{$i}) {
                    if ($cols{$i} eq 'GENEID') {
                        push(@{$self->{geneid}}, $word);
                    }
                    elsif ($cols{$i} eq 'NAME') {
                        push(@{$self->{genename}}, $word);
                    }
                    elsif ($cols{$i} eq 'GWEIGHT') {
                        push(@{$self->{gweight}}, $word);
                    }
                    elsif ($cols{$i} eq 'GORDER') {
                        push(@{$self->{gorder}}, $word);
                    }
                }
                else {
                    if ($word) {
                        push(@rowdata, $word);
                        push(@rowmask, 1);
                    }
                    else {
                        push(@rowdata, 0.0);
                        push(@rowmask, 0);
                        $needmask = 1;
                    }
                }
            }
            push(@{$self->{data}}, [@rowdata]);
            push(@{$self->{mask}}, [@rowmask]);
        }
    }
    if (not $needmask) {
        $self->{mask} = undef;
    }
}


sub treecluster {
    my ($self, %args) = @_;
    my %default = (
        transpose  =>     0,
        dist       =>   'e',
        method     =>   'm',
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    if ($param{transpose}==0) {
        $param{weight} = $self->{eweight};
    }
    else {
        $param{weight} = $self->{gweight};
    }
    return Algorithm::Cluster::treecluster(%param);
}

sub kcluster {
    my ($self, %args) = @_;
    my %default = (
        nclusters  =>     2,
        transpose  =>     0,
        npass      =>     1,
        method     =>   'a',
        dist       =>   'e',
        initidalid => undef,
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    if ($param{transpose}==0) {
        $param{weight} = $self->{eweight};
    }
    else {
        $param{weight} = $self->{gweight};
    }
    return Algorithm::Cluster::kcluster(%param);
}

sub somcluster {
    my ($self, %args) = @_;
    my %default = (
        transpose  =>     0,
        nxgrid     =>     2,
        nygrid     =>     1,
        inittau    =>  0.02,
        niter      =>     1,
        dist       =>   'e',
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    if ($param{transpose}==0) {
        $param{weight} = $self->{eweight};
    }
    else {
        $param{weight} = $self->{gweight};
    }
    return Algorithm::Cluster::somcluster(%param);
}

sub clustercentroids {
    my ($self, %args) = @_;
    my %default = (
        clusterid  => undef,
        method     =>   'a',
        transpose  =>     0,
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    my @data = @{$self->{data}};
    return Algorithm::Cluster::clustercentroids(%param);
}

sub clusterdistance {
    my ($self, %args) = @_;
    my %default = (
        cluster1   =>   [0],
        cluster2   =>   [0],
        method     =>   'a',
        dist       =>   'e',
        transpose  =>     0,
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    if ($param{transpose}==0) {
        $param{weight} = $self->{eweight};
    }
    else {
        $param{weight} = $self->{gweight};
    }
    return Algorithm::Cluster::clusterdistance(%param);
}

sub distancematrix {
    my ($self, %args) = @_;
    my %default = (
        dist       =>   'e',
        transpose  =>     0,
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    if (defined $self->{mask}) {
        $param{mask} = $self->{mask};
    }
    if ($param{transpose}==0) {
        $param{weight} = $self->{eweight};
    }
    else {
        $param{weight} = $self->{gweight};
    }
    return Algorithm::Cluster::distancematrix(%param);
}

sub save {
    my ($self, %args) = @_;
    my %default = (
        geneclusters => undef,
        expclusters  => undef,
    );
    my %param = (%default, %args);
    $param{data} = $self->{data};
    my $ngenes = scalar @{$self->{geneid}};
    my $nexps = scalar @{$self->{expid}};
    my $jobname = $param{jobname};
    defined ($jobname) or die 'jobname undefined';
    my $geneclusters;
    my $expclusters;
    my $gene_cluster_type;
    my $exp_cluster_type;
    if (defined $param{geneclusters}) {
        $geneclusters = $param{geneclusters};
        if (ref($geneclusters) eq "ARRAY") {
            if (scalar @{$geneclusters} != $ngenes) {
                die "k-means solution found, but its size does not agree with the number of genes";
            }
            $gene_cluster_type = 'k'; # k-means clustering result
        }
        elsif (ref($geneclusters) eq "Algorithm::Cluster::Tree") {
            $gene_cluster_type = 'h'; # hierarchical clustering result
            my $n = $geneclusters->length;
            if ($n != $ngenes - 1) {
                die "Size of the hierarchical clustering tree ($n) should be equal to the number of genes ($ngenes) minus one";
            }
        }
        else {
            die "Cannot understand gene clustering result! $!";
        }
    }
    if (defined $param{expclusters}) {
        $expclusters = $param{expclusters};
        if (ref($expclusters) eq "ARRAY") {
            if (scalar @$expclusters != $nexps) {
                die "k-means solution found, but its size does not agree with the number of experiments";
            }
            $exp_cluster_type = 'k'; # k-means clustering result
        }
        elsif (ref($expclusters) eq "Algorithm::Cluster::Tree") {
            $exp_cluster_type = 'h'; # hierarchical clustering result
            my $n = $expclusters->length;
            if ($n != $nexps - 1) {
                die "Size of the hierarchical clustering tree ($n) should be equal to the number of experiments ($nexps) minus one";
            }
        }
        else {
            die "Cannot understand experiment clustering result! $!";
        }
    }
    my @gorder;
    if (defined $self->{gorder}) {
        @gorder = $self->{gorder};
    }
    else {
        @gorder = (0..$ngenes-1);
    }
    my @eorder;
    if (defined $self->{eorder}) {
        @eorder = $self->{eorder};
    }
    else {
        @eorder = (0..$nexps-1);
    }
    if (defined $gene_cluster_type and defined $exp_cluster_type) {
        if ($gene_cluster_type ne $exp_cluster_type) {
            die 'found one k-means and one hierarchical clustering solution in geneclusters and expclusters';
        }
    }
    my $gid = 0;
    my $aid = 0;
    my $filename = $jobname;
    my $postfix = '';
    my @geneindex;
    my @expindex;
    if ($gene_cluster_type eq 'h') {
        # Hierarchical clustering result
        @geneindex = _savetree(jobname   => $jobname,
                               tree      => $geneclusters,
                               order     => \@gorder,
                               transpose =>  0);
        $gid = 1;
    }
    elsif ($gene_cluster_type eq 'k') {
        # k-means clustering result
        $filename = $jobname . '_K';
        my $k = -1;
        foreach (@$geneclusters) {
            if ($_ > $k) {
                $k = $_;
            }
        }
        $k++;
        my $kggfilename = $jobname . "_K_G$k.kgg";
        @geneindex = $self->_savekmeans(filename => $kggfilename,
                                        clusterids => \@$geneclusters,
                                        order => \@gorder,
                                        transpose => 0);
        $postfix = "_G$k";
    }
    else {
        @geneindex = sort { $gorder[$a] <=> $gorder[$b] } (0..$ngenes-1);
    }
    if ($exp_cluster_type eq 'h') {
        # Hierarchical clustering result
        @expindex = _savetree(jobname   => $jobname,
                              tree      => $expclusters,
                              order     => \@eorder,
                              transpose =>  1);
        $aid = 1;
    }
    elsif ($exp_cluster_type eq 'k') {
        # k-means clustering result
        $filename = $jobname . '_K';
        my $k = -1;
        foreach (@$expclusters) {
            if ($_ > $k) {
                $k = $_;
            }
        }
        $k++;
        my $kagfilename = $jobname . "_K_A$k.kag";
        @expindex = $self->_savekmeans(filename => $kagfilename,
                                       clusterids => \@$expclusters,
                                       order => \@eorder,
                                       transpose => 1);
        $postfix = $postfix . "_A$k";
    }
    else {
        @expindex = sort { $eorder[$a] <=> $eorder[$b] } (0..$nexps-1);
    }
    $filename = $filename . $postfix;
    $self->_savedata(jobname    => $filename,
                     gid        => $gid,
                     aid        => $aid,
                     geneindex  => \@geneindex,
                     expindex   => \@expindex);
}

sub _treesort {
    my %param = @_;
    my @order = @{$param{order}};
    my @nodeorder = @{$param{nodeorder}};
    my @nodecounts = @{$param{nodecounts}};
    my $tree = $param{tree};
    my $nNodes = $tree->length;
    my $nElements = $nNodes + 1;
    my @neworder = (0.0) x $nElements;
    my @clusterids = (0..$nElements-1);
    for (my $i = 0; $i < $nNodes; $i++) {
        my $i1 = $tree->get($i)->left;
        my $i2 = $tree->get($i)->right;
        my ($order1, $order2, $count1, $count2);
        if ($i1 < 0) {
            $order1 = $nodeorder[-$i1-1];
            $count1 = $nodecounts[-$i1-1];
        }
        else {
            $order1 = $order[$i1];
            $count1 = 1;
        }
        if ($i2 < 0) {
            $order2 = $nodeorder[-$i2-1];
            $count2 = $nodecounts[-$i2-1];
        }
        else {
            $order2 = $order[$i2];
            $count2 = 1;
        }
        # If order1 and order2 are equal, their order is determined by the order in which they were clustered
        my $increase;
        if ($i1 < $i2) {
            if ($order1 < $order2) {
                $increase = $count1;
            }
            else {
                $increase = $count2;
            }
            for (my $j = 0; $j < $nElements; $j++)
            {
                my $clusterid = $clusterids[$j];
                if ($clusterid==$i1 and $order1>=$order2) {
                    $neworder[$j] += $increase;
                }
                if ($clusterid==$i2 and $order1<$order2) {
                    $neworder[$j] += $increase;
                }
                if ($clusterid==$i1 or $clusterid==$i2) {
                    $clusterids[$j] = -$i-1;
                }
            }
        }
        else {
            if ($order1<=$order2) {
                $increase = $count1;
            }
            else {
                $increase = $count2;
            }
            for (my $j = 0; $j < $nElements; $j++) {
                my $clusterid = $clusterids[$j];
                if ($clusterid==$i1 and $order1>$order2) {
                    $neworder[$j] += $increase;
                }
                if ($clusterid==$i2 and $order1<=$order2) {
                    $neworder[$j] += $increase;
                }
                if ($clusterid==$i1 or $clusterid==$i2) {
                    $clusterids[$j] = -$i-1;
                }
            }
        }
    }
    my @result = sort { $neworder[$a] <=> $neworder[$b] } (0..$nElements-1);
    return @result;
}

sub _savetree {
    my %param = @_;
    my $jobname = $param{jobname};
    my $tree = $param{tree};
    my @order = @{$param{order}};
    my $transpose = $param{transpose};
    my ($extension, $keyword);
    if ($transpose==0) {
        $extension = 'gtr';
        $keyword = 'GENE';
    }
    else {
        $extension = 'atr';
        $keyword = 'ARRY';
    }
    my $nnodes = $tree->length;
    open OUTPUT, ">$jobname.$extension" or die 'Error: Unable to open output file';
    my @nodeID = ('') x $nnodes;
    my @nodecounts = (0) x $nnodes;
    my @nodeorder = (0.0) x $nnodes;
    my @nodedist;
    my $i;
    for ($i = 0; $i < $nnodes; $i++) {
        my $node = $tree->get($i);
        push (@nodedist, $node->distance);
    }
    for (my $nodeindex = 0; $nodeindex < $nnodes; $nodeindex++) {
        my $min1 = $tree->get($nodeindex)->left;
        my $min2 = $tree->get($nodeindex)->right;
        my $order1;
        my $order2;
        my $counts1;
        my $counts2;
        $nodeID[$nodeindex] = "NODE" . ($nodeindex+1) . "X";
        print OUTPUT $nodeID[$nodeindex];
        print OUTPUT "\t";
        if ($min1 < 0) {
            my $index1 = -$min1-1;
            $order1 = $nodeorder[$index1];
            $counts1 = $nodecounts[$index1];
            print OUTPUT $nodeID[$index1];
            print OUTPUT "\t";
            if ($nodedist[$index1] > $nodedist[$nodeindex]) {
            	$nodedist[$nodeindex] = $nodedist[$index1];
            }
        }
        else {
            $order1 = $order[$min1];
            $counts1 = 1;
            print OUTPUT $keyword . $min1 . "X\t";
        }
        if ($min2 < 0) {
            my $index2 = -$min2-1;
            $order2 = $nodeorder[$index2];
            $counts2 = $nodecounts[$index2];
            print OUTPUT $nodeID[$index2];
            print OUTPUT "\t";
            if ($nodedist[$index2] > $nodedist[$nodeindex]) {
            	$nodedist[$nodeindex] = $nodedist[$index2];
            }
        }
        else {
            $order2 = $order[$min2];
            $counts2 = 1;
            print OUTPUT $keyword . $min2 . "X\t";
        }
        print OUTPUT 1.0-$nodedist[$nodeindex];
        print OUTPUT "\n";
        $nodecounts[$nodeindex] = $counts1 + $counts2;
        $nodeorder[$nodeindex] = ($counts1*$order1+$counts2*$order2) / ($counts1+$counts2);
    }
    close(OUTPUT);
    # Now set up order based on the tree structure
    return _treesort(order      => \@order,
                     nodeorder  => \@nodeorder,
                     nodecounts => \@nodecounts,
                     tree       => $tree);
}

sub _savekmeans {
    my ($self, %param) = @_;
    my $filename = $param{filename};
    my @clusterids = @{$param{clusterids}};
    my @order = @{$param{order}};
    my $transpose = $param{transpose};
    my $label;
    my @names;
    if ($transpose == 0) {
        $label = $self->{uniqid};
        @names = @{$self->{geneid}};
    }
    else {
        $label = 'ARRAY';
        @names = @{$self->{expid}};
    }
    open OUTPUT, ">$filename" or die 'Error: Unable to open output file';
    print OUTPUT "$label\tGROUP\n";
    my $n = scalar @names;
    my @result = sort { $order[$a] <=> $order[$b] } (0..$n-1);
    my @sortedindex;
    my $cluster = 0;
    while (scalar @sortedindex < $n) {
        foreach (@result) {
            my $j = $_;
            my $cid = $clusterids[$j];
            if ($clusterids[$j]==$cluster) {
                print OUTPUT $names[$j] . "\t$cluster\n";
                push (@sortedindex, $j);
            }
        }
        $cluster++;
    }
    close(OUTPUT);
    return @sortedindex;
}

sub _savedata {
    my ($self, %param) = @_;
    my $jobname = $param{jobname};
    my $gid = $param{gid};
    my $aid = $param{aid};
    my @geneindex = @{$param{geneindex}};
    my @expindex = @{$param{expindex}};
    my @genename;
    if (defined $self->{genename}) {
        @genename = @{$self->{genename}};
    }
    else {
        @genename = @{$self->{geneid}};
    }
    my $ngenes = scalar @{$self->{geneid}};
    my $nexps = scalar @{$self->{expid}};
    open OUTPUT, ">$jobname.cdt" or die 'Error: Unable to open output file';
    my @mask;
    if (defined $self->{mask}) {
        @mask = @{$self->{mask}};
    }
    else {
        @mask = ([(1) x $nexps]) x $ngenes;
        # Each row contains identical shallow copies of the same vector;
        # modifying one row would affect the other rows.
    }
    my @gweight;
    if (defined $self->{gweight}) {
        @gweight = @{$self->{gweight}};
    }
    else {
        @gweight = (1) x $ngenes;
    }
    my @eweight;
    if (defined $self->{eweight}) {
        @eweight = @{$self->{eweight}};
    }
    else {
        @eweight = (1) x $nexps;
    }
    if ($gid) {
    	print OUTPUT "GID\t";
    }
    print OUTPUT $self->{uniqid};
    print OUTPUT "\tNAME\tGWEIGHT";
    # Now add headers for data columns
    foreach (@expindex) {
        print OUTPUT "\t" . $self->{expid}[$_];
    }
    print OUTPUT "\n";
    if ($aid) {
        print OUTPUT "AID";
        if ($gid) {
            print OUTPUT "\t";
        }
        print OUTPUT "\t\t";
        foreach (@expindex) {
            print OUTPUT "\tARRY" . $_ . 'X';
        }
        print OUTPUT "\n";
    }
    print OUTPUT "EWEIGHT";
    if ($gid) {
        print OUTPUT "\t";
    }
    print OUTPUT "\t\t";
    foreach (@expindex) {
        print OUTPUT "\t" . $eweight[$_];
    }
    print OUTPUT "\n";
    foreach (@geneindex) {
        my $i = $_;
        if ($gid) {
            print OUTPUT "GENE" . $i . "X\t";
        }
        print OUTPUT $self->{geneid}[$i] . "\t" . $genename[$i] . "\t" . $gweight[$i];
        foreach (@expindex) {
            my $j = $_;
            print OUTPUT "\t";
            if ($mask[$i][$j]) {
                print OUTPUT $self->{data}[$i][$j];
            }
        }
        print OUTPUT "\n";
    }
    close(OUTPUT);
}

1;
