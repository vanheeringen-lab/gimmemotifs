#!/usr/bin/perl

# $Id: overrepresentation.pl,v 1.5 2009/04/22 18:44:41 haudry Exp $

=head1 NAME

  overrepresentation.pl

=head1 SYNOPSIS

    this program takes 2 files [1] the sample file (in FASTA format) and [2] the background file (FASTA format)
    and will [1] run trawler and cluster all the motifs.
    the parameters used are described in the OPTIONS.

=head1 DESCRIPTION

This program is the main program runing Trawler.

=head1 OPTIONS

    -sample (FASTA format) better to be repeat-masked.
    -background (FASTA format)

    OPTIONAL PARAMETERS
    ===================

    [MOTIF DISCOVERY]
    -occurrence (optional) minimum occurrence in the sample sequences. [DEFAULT = 10]
    -wildcard (optional) number of wild card in motifs. [DEFAULT = 2]

    [VISUALIZATION]
    -directory (optional) output directory. [DEFAULT = "TRAWLER_HOME/myResults"]

=head1 CONTACT

Contact laurence Ettwiller (EMBL) ettwille@embl.de

=head1 DEPENDENCY

    perl / java
    trawler (from pecan)


=cut

use strict;
use Carp;
use File::Basename;
use Getopt::Long;

# Locate Trawler modules
use FindBin ();
use lib "$FindBin::RealBin/../modules";
my $script_name = $FindBin::RealScript;

# Internal Modules
use Trawler::Constants 1.0 qw(trawler_usage _read_config _tcst);

# START processing
print "\n## Running $script_name\n";

#==============================================================================
# Read and Set properties
#==============================================================================

# Read config file
_read_config($FindBin::RealBin);
my %tcst = _tcst();

# Logging Levels
my $DEBUG = $tcst{DEBUG};
my $INFO  = $tcst{INFO};


# TRAWLER input ===============================================================
# Required
my $DNA_file; # that's the sample sequence file in FASTA format
my $background_DNA_file; # that the background FASTA file for trawler
# Motif discovery options
my $occurrence = 3; # minimum occurrence of the motif
my $WC = 2; # number of wild cards in motifs
# Visualization options
my $DIR = undef;
# excluded for now
my $K_cutoff = 0; # SD from trawler default 0 (no filtering based on the standard deviation)
#=================================================================

GetOptions(
  'sample:s'     => \$DNA_file,
  'background:s' => \$background_DNA_file,
  'occurrence:s' => \$occurrence, # default 3
  'wildcard:s'   => \$WC, # default 2
  'directory:s'  => \$DIR,
  'k:s'          => \$K_cutoff,
);

# Whack out some help if lost...
if(!defined $DNA_file || !defined $background_DNA_file) {
    print STDERR "===================\n";
    print STDERR "\nError: Invalid command line arguments\n";
    trawler_usage();
}

# if no directory provided => exit [initialized in trawler.pl]
unless($DIR) {
    print STDERR "\nERROR: no directory to run\n";
    exit(1);
}
# extract tmp_xxx
my $tmp_dir_name = fileparse($DIR);
# default result directory is like $TRAWLER_HOME/tmp_YYYY-MM-DD_HHhmm:ss/result
my $tmp_result_dir = File::Spec->catdir( $DIR, $tcst{RES_DIR_NAME} );


# .trawler file
my $trawler_file = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{TRAWLER_FILE_EXT});

#==============================================================================
#now that everthing is read, start to run trawler

#background and sample file and output file
run_trawler($background_DNA_file, $DNA_file, $trawler_file, $occurrence, $WC);

#==============================================================================
# Subroutines
#==============================================================================

sub run_trawler {

    my ($background, $sample, $outputfile, $occurrence, $WC) = @_;

    #my $command = "java -Xmx500m -cp pecan.jar bp.trawler.Trawler -E $sample -F $background -N $outputfile -J $occurrence -G ACGTMRWSYKN -I $WC";
    my $command = "java " . $tcst{JAVA_OPTS} . " -cp '" . $tcst{LIB_PECAN} . "' bp.trawler.Trawler -E '$sample' -F '$background' -N '$outputfile' -J $occurrence -G ACGTMRWSYKN -I $WC";
    print "trawler cmd: $command\n" if $INFO;

    #run the commande now
    system($command) == 0 or croak "FAILURE command $command:\n $!\n with return code $?\n";

}

1;
