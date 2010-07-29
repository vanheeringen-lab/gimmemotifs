# $Id: ETimer.pm,v 1.1 2008/08/19 08:04:05 haudry Exp $

=head1 NAME

Trawler::ETimer - Utility package to compute elapsed time.

=head1 DESCRIPTION

This module helps in calculating elapsed time.

=head1 CONTACT

Yannick Haudry - EMBL (haudry@embl.de)

=cut


#  Usage => see Ensembl doc
#  use Trawler::ETimer;
#  $timer = new Trawler::ETimer;
#  sleep 2;
#  print 'Time elapsed: ', $timer->elapsed(), " seconds.\n";

package Trawler::ETimer;

use strict;
use Carp;

our $VERSION = '1.00';

###############################################################################

# Constructor
sub new {
    my $class_name = shift; # was shift

    my $self = {
        _start => time,
        _last  => time,
    };

    bless $self, $class_name;
    return $self;
}

# accessor method for ETimer last call
sub last_call {
    my ($self, $last) = @_;
    $self->{_last} = $last if defined($last);
    return $self->{_last};
}

# get elapsed time (start = object create)
sub elapsed {
    my $self = shift;
    return _get_etime( time - $self->{_start} );
}

# get interval time (since last call)
sub interval {
    my $self = shift;
    my $current_time = time;

    my $interval = _get_etime( $current_time - $self->{_last} );
    $self->{_last} = $current_time;

    return $interval;
}

# process seconds to return in human readable format
sub _get_etime {

    my $time_diff = shift;

    my @int = (
        [ 'second', 1                ],
        [ 'minute', 60               ],
        [ 'hour',   60*60            ],
        [ 'day',    60*60*24         ],
        [ 'week',   60*60*24*7       ],
        [ 'month',  60*60*24*30.5    ],
        [ 'year',   60*60*24*30.5*12 ],
    );
    my $i = $#int;
    my @r;
    while ( ($i >= 0) && ($time_diff) ) {
        if ($time_diff / $int[$i] -> [1] >= 1) {
            push @r, sprintf "%d %s%s",
                         $time_diff / $int[$i] -> [1],
                         $int[$i]->[0],
                         ( sprintf "%d", $time_diff / $int[$i] -> [1] ) > 1
                             ? 's'
                             : '';
        }
        $time_diff %= $int[$i] -> [1];
        $i--;
    }

    my $etime;
    if (@r) {
      $etime = join ", ", @r;
    } else {
      $etime = "NA";
    }
    #warn sprintf "ETIME %d\n", $etime;
    return $etime;
}

#------------------------------------------------------------------------------

#warn "Trawler::ETimer successfully loaded!\n";
1;
__END__
