# $Id: FileUtils.pm,v 1.3 2008/08/15 08:23:01 haudry Exp $

=head1 NAME

Trawler::FileUtils - Common Filesystem Utility methods

=head1 DESCRIPTION

This module contains generic methods to manipulate files

=head1 CONTACT

Yannick Haudry - EMBL (haudry@embl.de)

=cut

package Trawler::FileUtils;

use strict;
use Carp;

require Exporter;
our $VERSION = '1.00';
our @ISA = qw(Exporter);

our @EXPORT = (); # Symbols to autoexport (:DEFAULT tag)
our @EXPORT_OK = qw(_get_tstamp); # Symbols to export on request()
our %EXPORT_TAGS = (); # Define names for sets of symbols - eg: TAG => [ qw(!name1 name2) ]

#==============================================================================
# Filesystem subroutines
#==============================================================================

# Trawler::FileUtils::get_tstamp(time) [PUBLIC]
# Returns human readable current date
# output: YYYY-MM-DD_HH.mm
sub _get_tstamp {
    my $string = shift;
    my ($second, $minute, $hour, $day, $month, $year, $wday, $yday, $isdst);

    # Prefill date fields with current date
    if($string) {
        ($second, $minute, $hour, $day, $month, $year, $wday, $yday, $isdst) = localtime($string);
    } else {
        ($second, $minute, $hour, $day, $month, $year, $wday, $yday, $isdst) = localtime(time);
    }

    # Adjust month and year (Human readable date)
    $month++;       # month numbering starts from 0 !
    $year += 1900;  # year is the number of years since 1900 !

    #print "date is $day $month $year\n";
    # eg. 2008-02-17_14h01
    my $date_string = sprintf "%02s-%02s-%02s_%02sh%02s_%02s", $year, $month, $day, $hour, $minute, $second;

    return $date_string;

} # end get_tstamp()

#------------------------------------------------------------------------------

#sub _add_dir_separator {
#
#    my $directory = shift;
#
#    if ($directory !~ /\/$/) {
#        $directory = $directory."/";
#    }
#
#    return $directory;
#}

#------------------------------------------------------------------------------

#warn "Trawler::FileUtils successfully loaded!\n";
1;
__END__
