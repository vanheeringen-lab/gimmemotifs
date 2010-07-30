#---------------------------------------------------------------------------

package Algorithm::Cluster;

#---------------------------------------------------------------------------
# Copyright (c) 2003 John Nolan. All rights reserved.
# This program is free software.  You may modify and/or
# distribute it under the same terms as Perl itself.
# This copyright notice must remain attached to the file.
#
# Algorithm::Cluster is a set of Perl wrappers around the
# C Clustering library.
#
#---------------------------------------------------------------------------
# The C clustering library for cDNA microarray data.
# Copyright (C) 2002 Michiel Jan Laurens de Hoon.
#
# This library was written at the Laboratory of DNA Information Analysis,
# Human Genome Center, Institute of Medical Science, University of Tokyo,
# 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
# Contact: mdehoon 'AT' gsc.riken.jp
# 
# The Algorithm::Cluster module for Perl was released under the same terms
# as the Perl Artistic license. See the file artistic.txt for details.
#---------------------------------------------------------------------------


use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS @EXPORT);
use vars qw($DEBUG);
use strict;
use DynaLoader;

require Exporter;

$VERSION     = '1.49';
$DEBUG       = 1;
@ISA         = qw(DynaLoader Exporter);

@EXPORT_OK = qw(
    mean 
    median 
    kcluster 
    kmedoids 
    somcluster 
    treecluster
    clusterdistance 
    clustercentroids 
    distancematrix 
    pca
);

use warnings::register;

bootstrap Algorithm::Cluster $VERSION;


#-------------------------------------------------------------
# Debugging functions
#
sub version {
    return _version();
}


#-------------------------------------------------------------
# Wrapper for printing warnings
#
sub module_warn {
    return unless warnings::enabled();
    warnings::warn("Algorithm::Cluster", join '', @_);
}

#-------------------------------------------------------------
# Make sure that the first parameter is a reference-to-array,
# whose first member is itself a reference-to-array, 
# and that that array has at least one member.
#
sub data_is_valid_matrix {
    unless (ref($_[0]) eq 'ARRAY') {
        module_warn("Wanted array reference, but got a reference to ",
                    ref($_[0]), ". Cannot parse matrix");
        return;
    }
    my $nrows = scalar @{ $_[0] };
    unless ($nrows > 0) {
        module_warn("Matrix has zero rows.  Cannot parse matrix");
        return;
    }
    my $firstrow =  $_[0]->[0];
    unless (defined $firstrow) {
        module_warn("First row in matrix is undef scalar (?). Cannot parse matrix",);
        return;
    }
    unless (ref($firstrow) eq 'ARRAY') {
        module_warn("Wanted array reference, but got a reference to ",
                     ref($firstrow), ". Cannot parse matrix");
        return;
    }
    my $ncols = scalar @{ $_[0]->[0] };
    unless ($ncols > 0) {
        module_warn("Row has zero columns. Cannot parse matrix");
        return;
    }
    unless (defined($_[0]->[0]->[0])) {
        module_warn("Cell [0,0] is undefined. Cannot parse matrix");
        return;
    }
    return 1;
}


#-------------------------------------------------------------
# Wrapper for the mean() function
#
sub mean {
    if(ref $_[0] eq 'ARRAY') {
        return _mean($_[0]);
    } else {
        return _mean([@_]);
    }
}

#-------------------------------------------------------------
# Wrapper for the median() function
#
sub median {
    if(ref $_[0] eq 'ARRAY') {
        return _median($_[0]);
    } else {
        return _median([@_]);
    }
}


#------------------------------------------------------
# This function is called by the wrappers for library functions.
# It checks the dimensions of the data, mask and weight parameters.
#
# Return false if any errors are found in the data matrix. 
#
# Detect the dimension (nrows x ncols) of the data matrix,
# and set values in the parameter hash. 
#
# Also check the mask matrix and weight arrays, and set
# the parameters to default values if we find any errors, 
# however, we still return true if we find errors.
#
sub check_matrix_dimensions {
    my ($param, $default) = @_;
    #----------------------------------
    # Check the data matrix
    #
    return unless data_is_valid_matrix($param->{data});
    #----------------------------------
    # Remember the dimensions of the weight array
    #
    $param->{nrows}   = scalar @{ $param->{data}      };
    $param->{ncols}   = scalar @{ $param->{data}->[0] };
    #----------------------------------
    # Check the mask matrix
    #
    unless (data_is_valid_matrix($param->{mask})) {
        module_warn("Parameter 'mask' is not a valid matrix, ignoring it.");
        $param->{mask}      = $default->{mask}     
    } else {
        my $mask_nrows    = scalar @{ $param->{mask}      };
        my $mask_ncols    = scalar @{ $param->{mask}->[0] };
        unless ($param->{nrows} == $mask_nrows and $param->{ncols} == $mask_ncols) {
            module_warn("Data matrix is $param->{nrows}x$param->{ncols}, but mask matrix is ${mask_nrows}x${mask_ncols}.\nIgnoring the mask.");
            $param->{mask}      = $default->{mask};
        }
    }
    #----------------------------------
    # Check the weight array
    #
    unless(ref $param->{weight} eq 'ARRAY') {
            module_warn("Parameter 'weight' does not point to an array, ignoring it.");
            $param->{weight} = $default->{weight};
    } else {
        my $weight_length    = scalar @{ $param->{weight} };
        if ($param->{transpose} eq 0) {
            unless ($param->{ncols} == $weight_length) {
                module_warn("Data matrix has $param->{ncols} columns, but weight array has $weight_length items.\nIgnoring the weight array.");
                $param->{weight}      = $default->{weight}     
            }
        }
        else {
            unless ($param->{nrows} == $weight_length) {
                module_warn("Data matrix has $param->{nrows} rows, but weight array has $weight_length items.\nIgnoring the weight array.");
                $param->{weight}      = $default->{weight}     
            }
        }
    }
    return 1;
}


sub check_distance_matrix {
    my $distances = $_[0];
    my $i;
    my $row;
    my $column;
    #----------------------------------
    # Check the data matrix
    #
    my $reference = ref($distances);
    if (!$reference) {
        return "Wanted array reference but did not receive a reference";
    }
    elsif ($reference ne 'ARRAY') {
        return "Wanted array reference, but got a $reference";
    }
    my $nobjects = scalar @{ $distances };
    unless ($nobjects > 0) {
        return "Distance matrix has zero rows";
    }
    $i = 0;
    foreach $row (@{ $distances}) {
        unless (defined $row) {
            return "Row $i is undefined";
        }
        unless (ref($row) eq 'ARRAY') {
            return "Row $i is not an array";
        }
        unless (@{$row} == $i) {
            return "Row $i has incorrect columns";
        }
        foreach $column (@{$row}) {
            unless (defined($column)) {
                return "Row $i contains undefined columns";
            }
        }
        $i++;
    }
    return "OK";
}

sub check_initialid {
    my ($param, $default, $nobjects) = @_;
    my $i;
    my @counter = {};
    #----------------------------------
    # Check the initial clustering solution, if specified
    #
    if(ref $param->{initialid} ne 'ARRAY') {
        module_warn("Optional parameter 'initialid' should be an array");
        return;
    }
    if (@{ $param->{initialid}} == 0) {
        # no initial clustering solution specified
        if ($param->{nclusters}==-1) {
            $param->{nclusters} = 2; # default value
        }
        if ($param->{nclusters} > $nobjects) {
            module_warn("More clusters requested than elements available");
            return;
        }
        unless($param->{npass} =~ /^\d+$/ and $param->{npass} > 0) {
            module_warn("Parameter 'npass' must be a positive integer (got '$param->{npass}')");
            return;
        }
        return 1;
    }
    if (@{ $param->{initialid}} != $nobjects) {
        module_warn("Optional parameter 'initialid' should contain $nobjects elements");
        return;
    }
    foreach $i (@{ $param->{initialid}}) {
        unless($i =~ /^\d+$/ and $i >= 0) {
            module_warn("Optional parameter 'initialid' should only contain non-negative integers");
            return;
        }
    }
    if ($param->{nclusters} == -1) {
        # number of clusters was not specified. Infer it from initialid
        foreach $i (@{ $param->{initialid}}) {
            if ($i > $param->{nclusters}) {
                $param->{nclusters} = $i;
            }
        }
        $param->{nclusters}++;
    } else {
        # check if initialid is consistent with number of clusters
        foreach $i (@{ $param->{initialid}}) {
            if ($i >= $param->{nclusters}) {
                module_warn("Optional parameter 'initialid' inconsistent with nclusters");
                return;
            }
        }
    }
    # Check that none of the clusters are empty
    for ($i = 0; $i < $param->{nclusters}; $i++) {
        push(@counter, 0);
    }
    foreach $i (@{ $param->{initialid}}) {
        $counter[$i]++;
    }
    for ($i = 0; $i < $param->{nclusters}; $i++) {
        if ($counter[$i]==0) {
            module_warn("Optional parameter 'initialid' contains empty clusters");
            return;
        }
    }
    # No errors detected
    $param->{npass} = 0;
    return 1;
}

#-------------------------------------------------------------
# Wrapper for the kcluster() function
#
sub kcluster {
    #----------------------------------
    # Define default parameters
    #
    my %default = (
        nclusters =>    -1,
        data      =>  [[]],
        mask      =>    '',
        weight    =>    '',
        transpose =>     0,
        npass     =>     1,
        method    =>   'a',
        dist      =>   'e',
        initialid =>    [],
    );
    #----------------------------------
    # Local variable
    #
    my $nobjects = 0;
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    my @data = @{$param{data}};
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    return unless check_matrix_dimensions(\%param, \%default);
    #----------------------------------
    # Check the transpose parameter
    #
    if ($param{transpose} == 0) {
        $nobjects = $param{nrows};
    } elsif ($param{transpose} == 1) {
        $nobjects = $param{ncols};
    } else {
        module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
        return;
    }
    #----------------------------------
    # Check the initial clustering, if specified, and npass
    #
    return unless check_initialid(\%param, \%default, $nobjects);
    #----------------------------------
    # Check the other parameters
    #
    unless($param{method}    =~ /^[am]$/) {
        module_warn("Parameter 'method' must be either 'a' or 'm' (got '$param{method}')");
        return;
    }
    unless($param{dist}      =~ /^[cauxskeb]$/) {
        module_warn("Parameter 'dist' must be one of: [cauxskeb] (got '$param{dist}')");
        return;
    }
    #----------------------------------
    # Invoke the library function
    #
    return _kcluster(@param{
        qw/nclusters nrows ncols data mask weight transpose npass method dist initialid/
    });
}

#-------------------------------------------------------------
# Wrapper for the kmedoids() function
#

sub kmedoids {
    #----------------------------------
    # Define default parameters
    #
    my %default = (

        nclusters =>     2,
        distances =>  [[]],
        npass     =>     1,
        initialid =>    [],
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the distance matrix
    #
    my $message = check_distance_matrix($param{distances});
    unless ($message eq "OK") {
        module_warn($message); 
        return;
    }
    $param{nobjects} = scalar @{ $param{distances} }; 
    #----------------------------------
    # Check the initial clustering, if specified, and npass
    #
    return unless check_initialid(\%param, \%default, $param{nobjects});
    #----------------------------------
    # Invoke the library function
    #
    return _kmedoids(@param{
        qw/nclusters nobjects distances npass initialid/
    });
}

#-------------------------------------------------------------
# treecluster(): Wrapper for the treecluster function
#
sub treecluster {
    #----------------------------------
    # Define default parameters
    #
    my %default = (
        data       =>  [[]],
        mask       =>    '',
        weight     =>    '',
        transpose  =>     0,
        dist       =>   'e',
        method     =>   's',
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    my $message = check_distance_matrix($param{data});
    if ($message eq "OK") {
        $param{nrows}     = scalar @{ $param{data} }; 
        $param{ncols}     = scalar @{ $param{data} }; 
        $param{mask}      = $default{mask};
        $param{weight}    = $default{weight};
        $param{transpose} = $default{transpose};
        $param{dist}      = $default{dist};
        #----------------------------------
        # Check the clustering method
        #
        unless($param{method}    =~ /^[sma]$/) {
            module_warn("Parameter 'method' must be one of [sma] (got '$param{method}')");
            return;
        }
    } else {
        return unless check_matrix_dimensions(\%param, \%default);
        unless($param{transpose} =~ /^[01]$/) {
            module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
            return;
        }
        unless($param{dist}      =~ /^[cauxskeb]$/) {
            module_warn("Parameter 'dist' must be one of: [cauxskeb] (got '$param{dist}')");
            return;
        }
        unless($param{method}    =~ /^[smca]$/) {
            module_warn("Parameter 'method' must be one of [smca] (got '$param{method}')");
            return;
        }
    }
    #----------------------------------
    # Invoke the library function
    #
    return _treecluster(@param{
        qw/nrows ncols data mask weight transpose dist method/
    });
}



#-------------------------------------------------------------
# Wrapper for the clusterdistance() function
#
sub clusterdistance {
    #----------------------------------
    # Define default parameters
    #
    my %default = (
        data      =>  [[]],
        mask      =>    '',
        weight    =>    '',
        cluster1  =>    [],
        cluster2  =>    [],
        dist      =>   'e',
        method    =>   'a',
        transpose =>     0,
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the cluster1 and cluster2 arrays
    #
    if($param{cluster1} =~ /^\d+$/) {
        $param{cluster1} = [int($param{cluster1})];
    } elsif(ref $param{cluster1} ne 'ARRAY') {
        module_warn("Parameter 'cluster1' does not point to an array. Cannot compute distance.");
        return;
    } elsif(@{ $param{cluster1}} <= 0) {
        module_warn("Parameter 'cluster1' points to an empty array. Cannot compute distance.");
        return;
    }
    if($param{cluster2} =~ /^\d+$/) {
        $param{cluster2} = [int($param{cluster2})];
    } elsif(ref $param{cluster2} ne 'ARRAY') {
        module_warn("Parameter 'cluster2' does not point to an array. Cannot compute distance.");
        return;
    } elsif(@{ $param{cluster2}} <= 0) {
        module_warn("Parameter 'cluster2' points to an empty array. Cannot compute distance.");
        return;
    }
    $param{cluster1_len} = @{ $param{cluster1}};
    $param{cluster2_len} = @{ $param{cluster2}};
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    return unless check_matrix_dimensions(\%param, \%default);
    #----------------------------------
    # Check the other parameters
    #
    unless($param{transpose} =~ /^[01]$/) {
        module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
        return;
    }
    unless($param{method}    =~ /^[amsxv]$/) {
        module_warn("Parameter 'method' must be 'a', 'm', 's', 'x', or 'v' (got '$param{method}')");
        return;
    }
    unless($param{dist}      =~ /^[cauxskeb]$/) {
        module_warn("Parameter 'dist' must be one of: [cauxskeb] (got '$param{dist}')");
        return;
    }
    #----------------------------------
    # Invoke the library function
    #
    return _clusterdistance(@param{
        qw/nrows ncols data mask weight cluster1_len cluster2_len 
        cluster1 cluster2 dist method transpose/
    });
}


#-------------------------------------------------------------
# Wrapper for the clustercentroids() function
#
sub clustercentroids {
    #----------------------------------
    # Define default parameters
    #
    my %default = (
        data      =>  [[]],
        mask      =>    '',
        clusterid =>    [],
        method    =>   'a',
        transpose =>     0,
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    return unless check_matrix_dimensions(\%param, \%default);
    #----------------------------------
    # Check the other parameters
    #
    unless($param{transpose} =~ /^[01]$/) {
        module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
        return;
    }
    unless($param{method}    =~ /^[am]$/) {
        module_warn("Parameter 'method' must be 'a' or 'm' (got '$param{method}')");
        return;
    }
    #----------------------------------
    # Check the clusterid arrays
    #
    if($param{clusterid} =~ /^\d+$/) {
        $param{clusterid} = [int($param{clusterid})];
    } elsif(ref $param{clusterid} ne 'ARRAY') {
        module_warn("Parameter 'clusterid' does not point to an array. Cannot compute distance.");
        return;
    } elsif(@{ $param{clusterid}} <= 0) {
        module_warn("Parameter 'clusterid' points to an empty array. Cannot compute distance.");
        return;
    }
    my $clusterid_len = @{ $param{clusterid}};
    my $nrows = $param{nrows};
    my $ncols = $param{ncols};
    if ($param{transpose}==0 and $clusterid_len != $nrows) {
        die "Parameter 'clusterid' should have a size of $nrows; found $clusterid_len";
    }
    elsif ($param{transpose}==1 and $clusterid_len != $ncols) {
        die "Parameter 'clusterid' should have a size of $ncols; found $clusterid_len";
    }
    my $nclusters = -1;
    foreach (@{$param{clusterid}}) {
        if ($_ > $nclusters) {
            $nclusters = $_;
        }
    }
    $param{nclusters} = $nclusters + 1;
    #----------------------------------
    # Invoke the library function
    #
    return _clustercentroids(@param{
        qw/nclusters nrows ncols data mask clusterid transpose method/
    });
}


#-------------------------------------------------------------
# Wrapper for the distancematrix() function
#
sub distancematrix {
    #----------------------------------
    # Define default parameters
    #
    my %default = (

        data      =>  [[]],
        mask      =>    '',
        weight    =>    '',
        dist      =>   'e',
        transpose =>     0,
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    return unless check_matrix_dimensions(\%param, \%default);
    #----------------------------------
    # Check the transpose parameter
    #
    unless($param{transpose} =~ /^[01]$/) {
        module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
        return;
    }
    #----------------------------------
    # Check the other parameters
    #
    unless($param{dist}      =~ /^[cauxskeb]$/) {
        module_warn("Parameter 'dist' must be one of: [cauxskeb] (got '$param{dist}')");
        return;
    }
    #----------------------------------
    # Invoke the library function
    #
    return _distancematrix(@param{
        qw/nrows ncols data mask weight transpose dist/
    });
}


#-------------------------------------------------------------
# Wrapper for the somcluster() function
#
sub somcluster {
    #----------------------------------
    # Define default parameters
    #
    my %default = (
        data      =>  [[]],
        mask      =>    '',
        weight    =>    '',
        transpose =>     0,
        nxgrid    =>    10,
        nygrid    =>    10,
        inittau   =>  0.02,
        niter     =>   100,
        dist      =>   'e',
    );
    #----------------------------------
    # Accept parameters from caller
    #
    my %param = (%default, @_);
    #----------------------------------
    # Check the data, matrix and weight parameters
    #
    return unless check_matrix_dimensions(\%param, \%default);
    #----------------------------------
    # Check the other parameters
    #
    unless($param{transpose} =~ /^[01]$/) {
        module_warn("Parameter 'transpose' must be either 0 or 1 (got '$param{transpose}')");
        return;
    }
    unless($param{nxgrid}     =~ /^\d+$/ and $param{nxgrid} > 0) {
        module_warn("Parameter 'nxgrid' must be a positive integer (got '$param{nxgrid}')");
        return;
    }

    unless($param{nygrid}     =~ /^\d+$/ and $param{nygrid} > 0) {
        module_warn("Parameter 'nygrid' must be a positive integer (got '$param{nygrid}')");
        return;
    }
    unless($param{inittau}     =~ /^\d+.\d+$/ and $param{inittau} >= 0.0) {
        module_warn("Parameter 'inittau' must be a non-negative number (got '$param{inittau}')");
        return;
    }
    unless($param{niter}     =~ /^\d+$/ and $param{niter} > 0) {
        module_warn("Parameter 'niter' must be a positive integer (got '$param{niter}')");
        return;
    }
    unless($param{dist}      =~ /^[cauxskeb]$/) {
        module_warn("Parameter 'dist' must be one of: [cauxskeb] (got '$param{dist}')");
        return;
    }
    #----------------------------------
    # Invoke the library function
    #
    return _somcluster(@param{
        qw/nrows ncols data mask weight transpose nxgrid nygrid inittau niter dist/
    });
}


#-------------------------------------------------------------
# Wrapper for the pca() function
#
sub pca {
    #----------------------------------
    # Accept parameters from caller
    #
    my $data = shift;
    #----------------------------------
    # Check the data matrix
    #
    return unless data_is_valid_matrix($data);
    #----------------------------------
    # Remember the dimensions of the data array
    #
    my $nrows   = scalar @{$data};
    my $ncols   = scalar @{$data->[0]};

    #----------------------------------
    # Invoke the library function
    return _pca($nrows, $ncols, $data);
}


1;

__END__


=head1 NAME

Algorithm::Cluster - Perl interface to the C Clustering Library.


=head1 DESCRIPTION

This module is an interface to the C Clustering Library,
a general purpose library implementing functions for hierarchical 
clustering (pairwise simple, complete, average, and centroid linkage), 
along with k-means and k-medians clustering, and 2D self-organizing 
maps.  This library was developed at the Human Genome Center of the
University of Tokyo. The C Clustering Library is distributed along 
with Cluster 3.0, an enhanced version of the famous 
Cluster program originally written by Michael Eisen 
while at Stanford University.

=head1 EXAMPLES

See the scripts in the examples subdirectory of the package.

=head1 CHANGES

=over 4

=item * C Clustering Library version 1.49 (2010.04.05)

=head1 TO DO

=over

=head1 THANKS

Thanks to Michael Eisen, for creating the software packages
Cluster and TreeView. 

=head1 AUTHOR

John Nolan jpnolan@sonic.net 2003.  
Michiel de Hoon mdehoon "AT" gsc.riken.jp 2003-2010.
Seiya Imoto imoto "AT" ims.u-tokyo.ac.jp 2003-2010.
Satoru Miyano 2003-2010.
A copyright statement is contained in the source code itself. 

This module is a Perl wrapper for the C clustering library for 
cDNA microarray data, Copyright (C) 2002 Michiel Jan Laurens de Hoon.

See the source of Cluster.pm for a full copyright statement. 

=cut

1;
