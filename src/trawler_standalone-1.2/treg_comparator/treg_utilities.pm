#!/bin/perl -w

=head1 NAME

treg_utilities.pl

=head1 SYNOPSIS

=cut


package treg_utilities;

use strict;
use POSIX;

require Exporter;
use vars       qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = 0.90;
@ISA         = qw(Exporter);

@EXPORT      = qw(&testsub &get_reg_pmat_from_cmat &revcom_matrix &_prob_to_iupac &get_iupac_from_pmat);

%EXPORT_TAGS = ();     # eg: TAG => [ qw!name1 name2! ],
@EXPORT_OK   = ();
######################################################################

sub testsub  {
	my $s = "Hellow ow"; 
	return $s;
}

=pod

=item I<get_reg_pmat_from_cmat>

        Input is a reference to the count matrix (array). 
	Output is the regularized array.

=cut

sub get_reg_pmat_from_cmat   {   ### based on _calc_reg_p_mat_sophisticated
    my $ref = shift;  ## ref to count matrix
    my $len = scalar(@$ref)/4;   
    my @p = _calc_mat_dist($ref);
    my @c_mat = @{$ref};
    my ($i, $j);
    my ($delta, $deltaw);
    my ($v, $w);
    my $prec = 0.0001;
    my $cutoff = 1.5;
    my @reg_p_mat = ();

    # calculate raw_p_mat
    my @raw_p_mat = ();
    my $rs;
    for ($i = 0; $i < $len; $i++) {
	$rs = $c_mat[$i*4] + $c_mat[$i*4 + 1] + $c_mat[$i*4 + 2] + $c_mat[$i*4 + 3];
	if ($rs == 0) {
	    for ($j = 0; $j < 4; $j++) {
		$raw_p_mat[$i*4+$j] = $p[$j];
	    }
	} else {
	    for ($j = 0; $j < 4; $j++) {
		$raw_p_mat[$i*4+$j] = $c_mat[$i*4 + $j]/$rs;
	    }
	}
    }


    for ($i = 0; $i < $len; $i++) {
	$delta = 0;
	$rs = 0;
	for ($j = 0; $j < 4; $j++) {
	    if ($c_mat[$i*4 + $j] > 0) {
		$delta += $raw_p_mat[$i*4 + $j] * log($raw_p_mat[$i*4 + $j]/$p[$j]);
	    }
	    $rs += $c_mat[$i*4 + $j];
	}
	$delta *= 2 * $rs;
	
	if ($delta <= $cutoff) {
	    @reg_p_mat = (@reg_p_mat, @p);
	} else {
	    $deltaw = $delta;
	    $w = 0.5;
	    $v = 0.25;
	    while (abs($delta - $cutoff - $deltaw) > $prec) {
		$deltaw = 0;
		for ($j = 0; $j < 4; $j++) {
		    $deltaw += ((1 - $w)*$raw_p_mat[$i*4 + $j] + $w*$p[$j]) * log(((1 - $w)*$raw_p_mat[$i*4 + $j] + $w*$p[$j])/$p[$j]);
		}
		$deltaw *= 2 * $rs;
		if ($deltaw >= $delta - 1.5) {
		    $w += $v;
		} else {
		    $w -= $v;
		}
		$v *= 0.5;
	    }
	    for ($j = 0; $j < 4; $j++) {
		push(@reg_p_mat, (1 - $w)*$raw_p_mat[$i*4 + $j] + $w*$p[$j]);
	    }
	}
    }

    return @reg_p_mat;
}


sub _calc_mat_dist {
    my $ref = shift; 
    my $len = scalar(@$ref)/4;   
    my @c = (0,0,0,0);
    my @p = (0,0,0,0);
    my $n = 0;
    my $add = 0;
    for (my $i = 0; $i < $len; $i++) {
	for (my $j = 0; $j < 4; $j++) {
	    $c[$j] += $ref->[$i*4+$j];
	}
    }
    $n = $c[0] + $c[1] + $c[2] + $c[3];

    if ($c[0] == 0 || $c[1] == 0 || $c[2] == 0 || $c[3] == 0) {
	$add = 1;
    }
    for (my $j = 0; $j < 4; $j++) {
	$p[$j] = ( $c[$j] + $add/$len ) / ( $n + $add );
    }
    return @p;
}


sub revcom_matrix  {
    my $mref = shift;
    my @m = @{$mref};
    my $len = scalar(@m)/4;

    my @rvm = @m;
    for (my $i = 0; $i < $len; $i++) {
	$rvm[4*$i]   = $m[4*($len-$i-1)+3];
	$rvm[4*$i+1] = $m[4*($len-$i-1)+2];
	$rvm[4*$i+2] = $m[4*($len-$i-1)+1];
	$rvm[4*$i+3] = $m[4*($len-$i-1)];
    }

    return @rvm;
}


=item I<get_iupac>

        $iupac_string = $pssm->get_iupac();

On taking a pssm this returns the iupac sequence as a string based on the regularized p values.

=cut

#sub get_iupac  {
#    	my $ref_reg_p_mat=shift;
#	my @reg_p_mat = @{$ref_reg_p_mat};
#	my $pwm_length = scalar(@reg_p_mat) / 4;
#	my $iupac_seq;
#	my $count;
#	my $x=0;
#	for(my $i=0;$i<$pwm_length;$i++) {
#	    $x=0;
#	    my @c=($reg_p_mat[4*$i],$reg_p_mat[4*$i+1],$reg_p_mat[4*$i+2],$reg_p_mat[4*$i+3]);
#	    #print "@c";
#	    $iupac_seq .= _prob_to_iupac(@c);
#	}
#	return($iupac_seq);
#}
#

sub _prob_to_iupac {
    my @prob_dis=@_;
    my @iupac=();
    @iupac[0..3] = ("A", "C", "G", "T");
    $iupac[10] = 'M';   # [CA]
    $iupac[20] = 'R';   # [GA]
    $iupac[30] = 'W';   # [TA]
    $iupac[21] = 'S';   # [GC]
    $iupac[31] = 'Y';   # [TC]
    $iupac[32] = 'K';   # [TG]
    for(my $j=0;$j<4;$j++) {
	if($prob_dis[$j] > 0.5) {
	    return $iupac[$j];
	}
    }
    for(my $j=0; $j<3; $j++) {
	for(my $k=$j+1; $k<4; $k++) {
	    if(($prob_dis[$j]+$prob_dis[$k]) >= 0.75) {
		return $iupac[10*$k + $j];
	    }
	}
    }
    return 'N';
}


sub get_iupac_from_pmat  {
    	my $ref=shift;
    	my @reg_p_mat=@{$ref};
    	my $len = scalar(@reg_p_mat)/4;
	my $iupac_seq;
	my $count;
	my $x=0;
	for(my $i=0;$i<$len;$i++) {
	    $x=0;
	    my @c=($reg_p_mat[4*$i],$reg_p_mat[4*$i+1],$reg_p_mat[4*$i+2],$reg_p_mat[4*$i+3]);
	    #print "@c";
	    $iupac_seq .= _prob_to_iupac(@c);
	}
	return $iupac_seq;
}


1;
