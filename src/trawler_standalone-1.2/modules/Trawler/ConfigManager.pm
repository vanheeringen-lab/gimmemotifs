# $Id: ConfigManager.pm,v 1.4 2008/08/19 08:04:05 haudry Exp $

=head1 NAME

Trawler::ConfigManager - Manage property file

=head1 DESCRIPTION

This module read configuration data from a property file (key=value pairs)

=head1 CONTACT

Yannick Haudry - EMBL (haudry@embl.de)

=cut

package Trawler::ConfigManager;

use strict;
use Carp;

our $VERSION = '1.00';

###############################################################################

# Constructor
sub new {
    my ($class_name) = @_;

    my ($self) = {};

    bless ($self, $class_name);

    $self->{'_created'} = 1;
    return $self;
}

# read property file
sub read {
    my ($self, $file) = @_;

    open(CONFIGFILE, $file) or croak "Can't read config file '$file': $!";

    # We'll set a special property that tells what filename we just read.
    $self->{'_filename'} = $file;

    my $section;
    while(<CONFIGFILE>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?

        if ($_ =~ /^\[(.*)\]/) {
            $section = $1;
        } else {
            my ($config_name, $config_val) = split(/\s*=\s*/, $_, 2);
            if ($section) {
                $self->{"$section.$config_name"} = $config_val;
            } else {
                $self->{$config_name} = $config_val;
            }
        }

    } # end while

    close(CONFIGFILE) or croak "Can't close config file '$file': $!";

    return 1;

} # end read()

# Getter
sub getProperty {
    my ($self, $key) = @_;

    if (!defined $self->{$key}) {
        croak "property '$key' is not defined!\n";
    }
    return $self->{$key};
}

# Setter
sub setProperty {
    my ($self, $key, $value) = @_;

    $self->{$key} = $value;
}

#------------------------------------------------------------------------------

#warn "Trawler::ConfigManager successfully loaded!\n";
1;
__END__
