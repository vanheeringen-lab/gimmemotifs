# $Id: Constants.pm,v 1.15 2009/07/15 13:33:55 ramialis Exp $

=head1 NAME

Trawler::Constants - definition of Constants for Trawler scripts

=head1 DESCRIPTION

This module loads Trawler constants from a property file and put them in a hash

=head1 CONTACT

Yannick Haudry - EMBL (haudry@embl.de)

=cut

package Trawler::Constants;

use strict;
use Carp;
use File::Basename;
use File::Spec::Functions qw[catfile catdir];

use Trawler::ConfigManager 1.0;
#use Trawler::FileUtils 1.0 qw(_add_dir_separator);

require Exporter;
our $VERSION = '1.00';
our @ISA = qw(Exporter);

our @EXPORT = (); # Symbols to autoexport (:DEFAULT tag)
our @EXPORT_OK = qw(trawler_usage trawler_version _read_config _tcst); # Symbols to export on request()
our %EXPORT_TAGS = (); # Define names for sets of symbols - eg: TAG => [ qw(!name1 name2) ]

###############################################################################

# DIR AND FILES
use constant CONFIG_FILE  =>  'trawler.cfg';
use constant CONFIG_DIR   =>  'conf';
use constant BIN_DIR      =>  'bin';
use constant RESULTS_DIR  =>  'myResults';

# TRAWLER SCRIPTS
use constant SCRIPT_OVEREP  => 'overrepresentation.pl';
use constant SCRIPT_P1SOM   => 'pipeline_trawler_01_som.pl';
use constant SCRIPT_P2      => 'pipeline_trawler_02.pl';
use constant SCRIPT_P2NORTH => 'pipeline_trawler_02_no_orthologs.pl';
use constant SCRIPT_P3WEB   => 'pipeline_trawler_03.pl';

use constant SCRIPT_P1SOM_SSTR   => 'pipeline_trawler_01_som_single_strand.pl';
use constant SCRIPT_P2_SSTR      => 'pipeline_trawler_02_single_strand.pl';
use constant SCRIPT_P2NORTH_SSTR => 'pipeline_trawler_02_no_orthologs_single_strand.pl';

# TRAWLER LIBS
use constant LIB_SEQLOGO     => 'weblogo/seqlogo';
use constant LIB_TREG_COMP   => 'treg_comparator/treg_comparator.pl';
use constant LIB_TREG_MATRIX => 'treg_comparator/matrix_set_all';
use constant LIB_PECAN       => 'lib/pecan.jar';
use constant LIB_JALVIEW     => 'lib/jalviewApplet.jar';

###############################################################################

our %tcst = ();

sub _tcst {
    return %tcst;
}

#------------------------------------------------------------------------------

sub _read_config {

    my $bin_directory = shift;

    ### Configuration files to look up
    # Default configuration
    # First check with the env variable
    # NOT NEEDED
    #my $config_file = File::Spec->catfile( ( $ENV{TRAWLER_HOME}, CONFIG_DIR ), CONFIG_FILE);
    #my ($config_dir, $tr_path, $bin_dir);
    #($config_dir, $tr_path) = fileparse($config_file);

    # Then check with relative path from bin directory [Default]
    my ($config_file, $bin_dir, $tr_path);
    #if (!check_file_er($config_file)) {
    ($bin_dir, $tr_path) = fileparse($bin_directory);
    $config_file = File::Spec->catfile( ( $tr_path, CONFIG_DIR ), CONFIG_FILE);
    #}

    # Finally, If we can't resolve the config file path, throw error
    if (!check_file_er($config_file)) {
        croak "FAILURE configuration file $config_file not found or not readbale";
    }

    # Read config file (properties) [overrride properties]
    my $config = new Trawler::ConfigManager;
    $config->read($config_file) or die "Couldn't read default config file $config_file: $!";

    # Load user config (~/.trawlerc) if any !
    my $user_config = File::Spec->catfile( $ENV{HOME}, $config->getProperty('user.trawlerc') );
    if ( check_file_er($user_config) ) { # skip loading if not found
        print "loading user configuration: $user_config\n";
        $config->read($user_config) or die "Couldn't read user config file $user_config: $!";
    }

    ### Trawler Path setting [default: $TRAWLER_HOME environment variable]
    # Base path (trawler install directory) from config file
    my $base_path = $config->getProperty('user.basepath');
    unless($base_path) {
        $base_path = "$ENV{TRAWLER_HOME}";
    }
    # Try to guess base_path even if not configured using relative path
    unless($base_path) {
        $base_path = $tr_path;
    }
    # Finally check if we found the Trawler install dir
    unless(-e $base_path && -d $base_path) {
       die "Can not locate Trawler install directory: $base_path does not exists \n";
    }

    ### Now we can process the constants
    # bin PATH
    my $bin_path = File::Spec->catfile( $base_path, BIN_DIR );

    # Where do we store results ? [default: $TRAWLER_HOME]
    my $res_path = $config->getProperty('user.resultpath');
    unless($res_path) {
        $res_path = File::Spec->catdir( $base_path, RESULTS_DIR );
    }

    ### put constants in this hash
    %tcst = (
        # Paths
        BASE_PATH => $base_path,
        RES_PATH  => $res_path,

        # Log level
        DEBUG => $config->getProperty('logger.debug'), # 0:disabled 1:enabled
        INFO  => $config->getProperty('logger.info'),  # 0:disabled 1:enabled

        # trawler scripts
        overrepresentation                    => File::Spec->catfile( $bin_path, SCRIPT_OVEREP ),
        pipeline_trawler_01_som               => File::Spec->catfile( $bin_path, SCRIPT_P1SOM ),
        pipeline_trawler_02                   => File::Spec->catfile( $bin_path, SCRIPT_P2 ),
        pipeline_trawler_02_no_orthologs      => File::Spec->catfile( $bin_path, SCRIPT_P2NORTH ),
        pipeline_trawler_01_som_sstr          => File::Spec->catfile( $bin_path, SCRIPT_P1SOM_SSTR ),
        pipeline_trawler_02_sstr              => File::Spec->catfile( $bin_path, SCRIPT_P2_SSTR ),
        pipeline_trawler_02_no_orthologs_sstr => File::Spec->catfile( $bin_path, SCRIPT_P2NORTH_SSTR ),
        pipeline_trawler_03                   => File::Spec->catfile( $bin_path, SCRIPT_P3WEB ),

        # pipeline_trawler_01_som.pl constants (use absolute paths)
        SEQLOGO     => File::Spec->catfile( $base_path, LIB_SEQLOGO),
        TREG_COMP   => File::Spec->catfile( $base_path, LIB_TREG_COMP ),
        TREG_MATRIX => File::Spec->catfile( $base_path, LIB_TREG_MATRIX ),
        LIB_PECAN   => File::Spec->catfile( $base_path, LIB_PECAN ),
        LIB_JALVIEW => File::Spec->catfile( $base_path, LIB_JALVIEW ),
        JAVA_OPTS   => $config->getProperty('lib.JAVA_OPTS'),

        # internal
        trawler_version => $config->getProperty('internal.version'),

        # HTML links
        #group_url => $config->getProperty('html.group_url'),

        # Dir handling
        RES_DIR_PREFIX    => "tmp_",
        RES_DIR_NAME      => "result",
        FASTA_DIR_NAME    => "fasta",
        FEATURES_DIR_NAME => "features",
        INPUT_DIR_NAME    => "input",

        # HTML DIR
        HTML_CSS      => "css",
        HTML_DOWNLOAD => "downloads",
        HTML_IMG      => "images",
        HTML_JS       => "js",
        HTML_LIB      => "lib",
        HTML_INPUT    => "input",

        # INPUT FILES
        INPUT_FILE => "input.txt",
        README_FILE => "README.txt",
        LICENSE_FILE => "LICENSE.txt",

        # Files handling
        CULSTER_FILE_EXT       => ".cluster",
        GRAPH_FILE_EXT         => ".graph",
        TRAWLER_FILE_EXT       => ".trawler",
        TRAWLER_SHORT_FILE_EXT => ".trawler_short",
        COMPARE_FILE_EXT       => "_compare",
        PWM_FILE_EXT           => ".pwm",
        FASTA_FILE_EXT         => ".fasta",

        STAT_FILE_NAME   => "STAT.tmp",

        MOTIF_PNG_EXT => "_all_motif.png",
        HTML_EXT      => ".html",

        );

}


sub check_file_er {

   my $file_path = shift;

   if (-e $file_path && -r $file_path) {
     return 1;
   }
   else {
      return;
   }

}

sub trawler_usage {
my $usage = <<"USAGE";

    Usage: trawler.pl -sample [file containing the enriched sequences] -background [file containing the background sequences]

    -sample (FASTA format) better to be repeat-masked.
    -background (FASTA format)

    OPTIONAL PARAMETERS
    ===================

    [MOTIF DISCOVERY]
    -occurrence (optional) minimum occurrence in the sample sequences. [DEFAULT = 10]
    -mlength (optional) minimum motif length. [DEFAULT = 6]
    -wildcard (optional) number of wild card in motifs. [DEFAULT = 2]
    -strand (optional) single or double [DEFAULT = double]

    [CLUSTERING]
    -overlap (optional) in percentage. [DEFAULT = 70]
    -motif_number (optional) total number of motifs to be clustered. [DEFAULT = 200]
    -nb_of_cluster (optional) fixed number of cluster. if this option is set to an integer (1 and above), the k-mean clustering algorithm with fixed k will be used instead of the strongly connected component (SCC). [DEFAULT = NULL]

    [VISUALIZATION]
    -directory (optional) output directory. [DEFAULT = "TRAWLER_HOME/myResults"]
    -dir_id (optional) gives an id to the results directory. [DEFAULT = NULL]
    -xtralen (optional) add bases around the motifs for the logo. [DEFAULT = 0]
    -alignments (optional) file containing the list of files containing the aligned sequences (see README file for more info) [DEFAULT = NULL]
    -ref_species (optional) name of the reference species [DEFAULT = NULL]
    -clustering (optional) if 1 the program clusters the instances, if 0 no clustering. [DEFAULT = 1]
    -web (optional) if 1 the output will be a web page with all the information [DEFAULT = 1]

USAGE

    print STDERR $usage."\n";

    exit(1);
}

sub trawler_version {
    my $version = "trawler standalone \"$tcst{trawler_version}\"\n";
    print $version;
    exit(1);
}

#------------------------------------------------------------------------------

#warn "Trawler::Constants successfully loaded!\n";
1;
__END__
