#!/usr/bin/perl

# $Id: trawler.pl,v 1.30 2009/04/24 12:01:11 haudry Exp $

=head1 NAME

    trawler.pl

=head1 SYNOPSIS


=head1 DESCRIPTION

    file that does everything.

=head1 OPTIONS

    Usage: trawler.pl -sample [file containing the enriched sequences] -background [file containing the background sequences]

    -sample (FASTA format) better to be repeat-masked.
    -background (FASTA format)

    OPTIONAL PARAMETERS
    ===================

    [MOTIF DISCOVERY]
    -occurrence (optional) minimum occurrence in the sample sequences. [DEFAULT = 10]
    -mlength (optional) minimum motif length. [DEFAULT = 6]
    -wildcard (optional) number of wild card in motifs. [DEFAULT = 2]
    -strand (optional) single or double. [DEFAULT = double]

    [CLUSTERING]
    -overlap (optional) in percentage. [DEFAULT = 70]
    -motif_number (optional) total number of motifs to be clustered. [DEFAULT = 200]
    -nb_of_cluster (optional) fixed number of cluster. if this option is used, the k-mean clustering algorithm with fixed k will be used instead of the strongly connected component (SCC). Alternatively this option can be set to 'som' and in this case, the self organizing map (SOM) will be use. [DEFAULT = NULL]

    [VISUALIZATION]
    -directory (optional) output directory. [DEFAULT = "TRAWLER_HOME/myResults"]
    -dir_id (optional) gives an id to the results directory. [DEFAULT = NULL]
    -xtralen (optional) add bases around the motifs for the logo. [DEFAULT = 0]
    -alignments (optional) file containing the list of files containing the aligned sequences (see README file for more info) [DEFAULT = NULL]
    -ref_species (optional) name of the reference species [DEFAULT = NULL]
    -clustering (optional) if 1 the program clusters the instances, if 0 no clustering. [DEFAULT = 1]
    -web (optional) if 1 the output will be a web page with all the information [DEFAULT = 1]

=head1 CONTACT

Contact laurence Ettwiller (EMBL) ettwille@embl.de

=head1 DEPENDENCY

    perl / java
    weblogo (http://weblogo.berkeley.edu/)
    treg_comparator.pl

=cut

use 5.6.0;
use strict;
use Carp;
use File::Basename;
use File::Copy;
use File::Path;
use File::Spec::Functions qw[catfile catdir];
use Getopt::Long;
use Time::Local;

# Locate Trawler modules
use FindBin ();
use lib "$FindBin::RealBin/../modules";
my $script_name = $FindBin::RealScript;

# Internal Modules
use Trawler::Constants 1.0 qw(trawler_usage trawler_version _read_config _tcst);
use Trawler::FileUtils 1.0 qw(_get_tstamp);
use Trawler::ETimer 1.0;

# START processing
#print "\n## Running $script_name\n";

#==============================================================================
# Read and Set properties
#==============================================================================

# Read config file
_read_config($FindBin::RealBin);
my %tcst = _tcst();

# Logging Levels
my $DEBUG = $tcst{DEBUG};
my $INFO  = $tcst{INFO};

#==============================================================================
# Initialize options

# Help
my $help;    # print USAGE
my $version; # print VERSION

# Required
my $sample = undef;
my $background = undef;
# Required [with orthologs]
my $alignments; # list of file containing the aligned sequences explained in README
my $ref_species;

# Motif discovery
my $occurrence = 10;
my $mlength = 6;
my $wildcard = 2;
my $strand = 'double';

# Clustering
my $overlap = 70;
my $motif_number = 200;
my $nb_of_cluster = undef;

# Visualization
my $directory = undef;
my $dir_id = undef; # give meaningful name to the result dir
my $xtralen = 0;
my $clustering = 1;
my $web = 1;

my $k = 0;
my $seqlogo = 1;
my $avoid;
my $low_content = 0;
my $treg_threshold = 0.9;


#==============================================================================
# Get Options
my $trawler_input = join(" ", @ARGV);
############### TEST process input params
my @trawler_argv = @ARGV;

###############

GetOptions(
  'help'           => \$help,
  'version'        => \$version,
  # Required
  'sample=s'       => \$sample,
  'background=s'   => \$background,
  # Motif discovery
  'occurrence:s'   => \$occurrence,
  'mlength:s'      => \$mlength,
  'wildcard:s'     => \$wildcard,
  'strand:s'       => \$strand,
  # Clustering
  'overlap:s'      => \$overlap,
  'motif_number:s' => \$motif_number, # default 200
  'nb_of_cluster:s'=> \$nb_of_cluster, # default undef.
  # Visualization
  'directory=s'    => \$directory,
  'dir_id=s'       => \$dir_id,
  'xtralen:s'      => \$xtralen, # default 0
  'alignments:s'   => \$alignments,
  'ref_species:s'  => \$ref_species,
  'clustering:s'   => \$clustering, # default 1
  'web:s'          => \$web, # default 1
  # excluuded for now
  'k:s'            => \$k,
  'seqlogo:s'      => \$seqlogo, # default 1
  'avoid:s'        => \$avoid,
  'low_content:s'  => \$low_content,
  'treg_threshold' => \$treg_threshold, # default 0.9
);

# Unprocessed options
#print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
#foreach (@ARGV) {
#  print "$_\n";
#}

#==============================================================================
# Process options

# if help is called, then usage
if ($help) { trawler_usage(); }

# if version is called, print it
if ($version) { trawler_version(); }

# start Running Trawler
msg_config();

# if no sample and background, then usage
unless($sample && $background) {
    print "\nERROR: must specify a sample and background with -sample and -background\n";
    trawler_usage();
}

# Apply SOM by default
if (!defined $nb_of_cluster) {
    $nb_of_cluster = 0; # enables SOM
}

#==============================================================================
# START
my $timer = new Trawler::ETimer;
my $timer_stack;

#==============================================================================
# Directories handling

# create a default directory where everything will be stored (to change use option -directory).
unless($directory) {
    # set the base result dir path (Default: $TRAWER_HOME/myResults)
    $directory = $tcst{RES_PATH};
    print "default base result directory: $directory \n" if $DEBUG;
}
# handle tmp directory to store current results (can be customized with option -dir_id).
my $tmp_dir_name;
if ($dir_id) {
  # <id>_<localetime>/
  $tmp_dir_name = $dir_id . "_" . _get_tstamp(); # get timestamp to create resut dir
} else {
  # tmp_<localetime>/
  $tmp_dir_name = $tcst{RES_DIR_PREFIX} . _get_tstamp(); # get timestamp to create resut dir
}
# default tmp store directory is like $TRAWLER_HOME/tmp_YYYY-MM-DD_HHhmm:ss
$directory = File::Spec->catfile($directory, $tmp_dir_name);

# default procesing results directory is like $TRAWLER_HOME/tmp_YYYY-MM-DD_HHhmm:ss/result
my $tmp_result_dir = File::Spec->catdir( $directory, $tcst{RES_DIR_NAME} );
my $tmp_fasta_dir = File::Spec->catdir( $tmp_result_dir, $tcst{FASTA_DIR_NAME} );
my $tmp_features_dir = File::Spec->catdir( $tmp_result_dir, $tcst{FEATURES_DIR_NAME} );

my $input_file = create_output_arch();

# $TRAWLER_RESULT/result/file_motif.cluster
my $file_motif = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{CULSTER_FILE_EXT});


#==============================================================================
# Run Trawler
#==============================================================================

# double strand
if ($strand eq 'double') {

if ($clustering == 0) {
    my $c = "perl \"" . $tcst{overrepresentation} . "\" -sample \"$sample\" -background \"$background\" -directory \"$directory\" -k $k -occurrence $occurrence -wildcard $wildcard";
    print "LOG[INFO]: trawler will be run only for the over-representation function and the output will be stored in $directory\n the commande is:\n $c\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c':\n $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with program $c\n" if $DEBUG;
    $timer_stack .= "overrepresentation.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif ($clustering == 1 && $web == 0) {
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    print "LOG[INFO]: trawler will be run for the over-representation function and and the clustering function (no web) and the output will be stored in $directory\n the command is:\n $c\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif ($alignments && $clustering == 1 && $web == 1) { # then the first schema needs to be run.
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    my $c1 = "perl \"" . $tcst{pipeline_trawler_02} . "\" -directory \"$directory\" -list \"$alignments\" -ref_species \"$ref_species\" -motif \"$file_motif\"";
    my $c2 = "perl \"" . $tcst{pipeline_trawler_03} . "\" -directory \"$directory\" -conservation 1";
    print "LOG[INFO]: trawler will be run and the output will be stored in $directory\n the commands are :\n $c\n AND\n $c1\n AND\n $c2\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with first program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som.pl => " . $timer->interval() . "\n";

    system($c1) == 0 or croak "FAILURE command '$c1': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with second program $c1\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_02.pl => " . $timer->interval() . "\n";

    system($c2) == 0 or croak "FAILURE command '$c2': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with third program $2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_03.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif (!$alignments && $clustering == 1 && $web == 1) { # then the first schema needs to be run.
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    my $c1 = "perl \"" . $tcst{pipeline_trawler_02_no_orthologs} . "\" -directory \"$directory\" -sequences \"$sample\" -motif \"$file_motif\"";
    my $c2 = "perl \"" . $tcst{pipeline_trawler_03} . "\" -directory \"$directory\"";
    print "LOG[INFO]: trawler will be run (but no orthologous information provided) and the output will be stored in $directory\n the commands are :\n $c\n AND\n $c1\n AND\n $c2\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with first program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som.pl => " . $timer->interval() . "\n";

    my $t1 = time;
    system($c1) == 0 or croak "FAILURE command '$c1': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with second program $c2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_02_no_orthologs.pl => " . $timer->interval() . "\n";

    my $t2 = time;
    system($c2) == 0 or croak "FAILURE command '$c2': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with third program $c2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_03 => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
else {
    trawler_usage();
}

# single strand
} elsif ($strand eq 'single') {

if ($clustering == 0) {
    my $c = "perl \"" . $tcst{overrepresentation} . "\" -sample \"$sample\" -background \"$background\" -directory \"$directory\" -k $k -occurrence $occurrence -wildcard $wildcard";
    print "LOG[INFO]: trawler will be run only for the over-representation function and the output will be stored in $directory\n the commande is:\n $c\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c':\n $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with program $c\n" if $DEBUG;
    $timer_stack .= "overrepresentation.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif ($clustering == 1 && $web == 0) {
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som_sstr} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    print "LOG[INFO]: trawler will be run for the over-representation function and and the clustering function (no web) and the output will be stored in $directory\n the command is:\n $c\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som_single_strand.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif ($alignments && $clustering == 1 && $web == 1) { # then the first schema needs to be run.
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som_sstr} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    my $c1 = "perl \"" . $tcst{pipeline_trawler_02_sstr} . "\" -directory \"$directory\" -list \"$alignments\" -ref_species \"$ref_species\" -motif \"$file_motif\"";
    my $c2 = "perl \"" . $tcst{pipeline_trawler_03} . "\" -directory \"$directory\" -conservation 1";
    print "LOG[INFO]: trawler will be run and the output will be stored in $directory\n the commands are :\n $c\n AND\n $c1\n AND\n $c2\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with first program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som_single_strand.pl => " . $timer->interval() . "\n";

    system($c1) == 0 or croak "FAILURE command '$c1': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with second program $c1\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_02_single_strand.pl => " . $timer->interval() . "\n";

    system($c2) == 0 or croak "FAILURE command '$c2': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with third program $2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_03.pl => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
elsif (!$alignments && $clustering == 1 && $web == 1) { # then the first schema needs to be run.
    my $c = "perl \"" . $tcst{pipeline_trawler_01_som_sstr} . "\" -sample \"$sample\" -background \"$background\" -k $k -overlap $overlap -seqlogo $seqlogo -directory \"$directory\" -xtralen $xtralen -avoid $avoid -mlength $mlength -motif_number $motif_number -low_content $low_content -occurrence $occurrence -wildcard $wildcard -treg_threshold $treg_threshold -nb_of_cluster $nb_of_cluster";
    my $c1 = "perl \"" . $tcst{pipeline_trawler_02_no_orthologs_sstr} . "\" -directory \"$directory\" -sequences \"$sample\" -motif \"$file_motif\"";
    my $c2 = "perl \"" . $tcst{pipeline_trawler_03} . "\" -directory \"$directory\"";
    print "LOG[INFO]: trawler will be run (but no orthologous information provided) and the output will be stored in $directory\n the commands are :\n $c\n AND\n $c1\n AND\n $c2\n" if $INFO;
    $timer_stack .= "trawler.pl => " . $timer->interval() . "\n";

    system($c) == 0 or croak "FAILURE command '$c': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with first program $c\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_01_som_sstr.pl => " . $timer->interval() . "\n";

    my $t1 = time;
    system($c1) == 0 or croak "FAILURE command '$c1': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with second program $c2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_02_no_orthologs_sstr.pl => " . $timer->interval() . "\n";

    my $t2 = time;
    system($c2) == 0 or croak "FAILURE command '$c2': $!\n with return code $?\n";
    print "LOG[DEBUG]: finished with third program $c2\n" if $DEBUG;
    $timer_stack .= "pipeline_trawler_03 => " . $timer->interval() . "\n";

    $timer_stack .= "Elapsed time => " . $timer->elapsed() . "\n";
    msg_output();
}
else {
    trawler_usage();
}

}

create_input_file($input_file);


#==============================================================================
# Subroutines
#==============================================================================

sub create_output_arch {

  # FIXME[YH]: create full directory result structure here !
  # create Trawler results directories (result fasta features)
  eval { mkpath([$tmp_result_dir, $tmp_fasta_dir, $tmp_features_dir]) };
  if ($@) {
    print STDERR "Cannot create $tmp_result_dir: $@";
    exit(1);
  }

  # TODO[YH]: Refactor with prepare_html_output
  # TODO[YH]: use constants for file/dir naming
  # create HTML output directories
  my $html_css_dir = File::Spec->catdir( $directory, $tcst{HTML_CSS} );
  my $html_download_dir = File::Spec->catdir( $directory, $tcst{HTML_DOWNLOAD} );
  my $html_imgs_dir = File::Spec->catdir( $directory, $tcst{HTML_IMG} );
  my $html_js_dir = File::Spec->catdir( $directory, $tcst{HTML_JS} );
  my $html_lib_dir = File::Spec->catdir( $directory, $tcst{HTML_LIB} );
  my $html_input_dir = File::Spec->catdir( $directory, $tcst{HTML_INPUT} );
  eval { mkpath( [ $html_css_dir, $html_js_dir, $html_imgs_dir, $html_download_dir, $html_lib_dir, $html_input_dir ] ) };
  if ($@) {
    print STDERR "Cannot create HTML ouput directories: $@";
    exit(1);
  }

  # TODO[YH]: get arrays from Config (hash: key=>srcDir, value=>file array)
  # copy JS files
  my @js_files = ('base.js', 'excanvas.js');
  foreach (@js_files) {
    copy( File::Spec->catfile( ($tcst{BASE_PATH}, "jsscripts"), $_),
        File::Spec->catfile( $html_js_dir, $_) )
        or die "Copy failed: $!";
  }
  # copy CSS files
  my @css_files = ('style.css', 'asc.gif', 'bg.gif', 'desc.gif', 'loading.gif', 'tab.png', 'hover.png');
  foreach (@css_files) {
    copy( File::Spec->catfile( ($tcst{BASE_PATH}, "css"), $_),
          File::Spec->catfile( $html_css_dir, $_) )
          or die "Copy failed: $!";
  }
  # copy LIB files
  copy( File::Spec->catfile( ($tcst{BASE_PATH}, "lib"), "jalviewApplet.jar"),
        File::Spec->catfile( $html_lib_dir, "jalviewApplet.jar") )
        or die "Copy failed: $!";

  # copy README file
  copy( File::Spec->catfile( $tcst{BASE_PATH}, $tcst{README_FILE} ),
        File::Spec->catfile( $html_input_dir, $tcst{README_FILE} ) )
        or die "Copy failed: $!";

  # copy LICENSE file
  copy( File::Spec->catfile( $tcst{BASE_PATH}, $tcst{LICENSE_FILE} ),
        File::Spec->catfile( $html_input_dir, $tcst{LICENSE_FILE} ) )
        or die "Copy failed: $!";

  my $input_f = File::Spec->catfile( $html_input_dir, $tcst{INPUT_FILE} );
  return $input_f;

}

# trawler options
sub get_trawler_opts {

    my %tparams = (
        # Required
        'sample' => '',
        'background' => '',
        # Motif discovery
        'occurrence' => '',
        'mlength' => '',
        'wildcard' => '',
        'strand' => '',
        # Clustering
        'overlap' => '',
        'motif_number' => '',
        'nb_of_cluster' => '', # new
        # Visualization
        'directory' => '',
        'dir_id' => '',
        'xtralen' => '',
        'alignments' => '',
        'ref_species' => '',
        'clustering' => '',
        'web' => '',
        #'k' => '',
        #'seqlogo' => '',
        #'avoid' => '',
        #'low_content' => '',
        #'treg_threshold' => '',
    );

    return %tparams;
}

# create input summary file
sub create_input_file {

  my $input_file = shift;

  #==============================================================================
  # Create Input file (HTML output).
  # TODO[YH]: create params.txt like file

  my $input_opts = create_input_sum();

  print "### Runnning with options: " . $trawler_input . "\n" if $DEBUG;
  open(FINPUT, ">$input_file") or croak "Cannot open file $input_file: $!";
  print FINPUT "\$ trawler.pl " . $trawler_input . "\n"; # command line
  print FINPUT $input_opts; # trawler options
  print FINPUT "\npipeline\n";
  print FINPUT "========\n\n";
  print FINPUT $timer_stack; # pipeline and etime
  close(FINPUT) or croak "Can't close file '$input_file': $!";

}

sub create_input_sum {

  my $opts;
  print "trawler_input(ARGV): $trawler_input\n" if $DEBUG;

  my %params = get_trawler_opts();

  # extract command line options
  $opts .= "\ncommand line options\n";
  $opts .= "====================\n\n";
  my %cmd_opts = {};
  my $len = $#trawler_argv;
  my $i = 0;
  while ($i < $len) {
    $trawler_argv[$i] =~ s/\-+//g; # remove option switch
    if(exists $params{$trawler_argv[$i]}) { # get option=value
      $opts .= "$trawler_argv[$i]=$trawler_argv[$i+1]\n";
      $cmd_opts{$trawler_argv[$i]} = $trawler_argv[$i+1];
      $i++;
    }
    $i++;
  }

  # then extract default values
  $opts .= "\ndefault options\n";
  $opts .= "===============\n";
  $opts .= "\n[motif discovery]\n";
  $opts .= "occurrence=" . $occurrence . "\n" if defined($occurrence) && !exists($cmd_opts{"occurrence"});
  $opts .= "mlength=" . $mlength . "\n" if defined($mlength) && !exists($cmd_opts{"mlength"});
  $opts .= "wildcard=" . $wildcard . "\n" if defined($wildcard) && !exists($cmd_opts{"wildcard"});
  $opts .= "strand=" . $strand . "\n" if defined($strand) && !exists($cmd_opts{"strand"});
  $opts .= "\n[clustering]\n";
  $opts .= "overlap=" . $overlap . "\n" if defined($overlap) && !exists($cmd_opts{"overlap"});
  $opts .= "motif_number=" . $motif_number . "\n" if defined($motif_number) && !exists($cmd_opts{"motif_number"});
  $opts .= "nb_of_cluster=" . $nb_of_cluster . "\n" if defined($nb_of_cluster) && !exists($cmd_opts{"nb_of_cluster"});
  $opts .= "\n[visualization]\n";
  $opts .= "directory=" . $directory . "\n" if defined($directory) && !exists($cmd_opts{"directory"});
  $opts .= "dir_id=" . $dir_id . "\n" if defined($dir_id) && !exists($cmd_opts{"dir_id"});
  $opts .= "xtralen=" . $xtralen . "\n" if defined($xtralen) && !exists($cmd_opts{"xtralen"});
  $opts .= "alignments=" . $alignments . "\n" if defined($alignments) && !exists($cmd_opts{"alignments"});
  $opts .= "ref_species=" . $ref_species . "\n" if defined($ref_species) && !exists($cmd_opts{"ref_species"});
  $opts .= "clustering=" . $clustering . "\n" if defined($clustering) && !exists($cmd_opts{"clustering"});
  $opts .= "web=" . $web . "\n" if defined($web) && !exists($cmd_opts{"web"});
  #$opts .= "k=" . $k . "\n" if defined($k) && !exists($cmd_opts{"k"});
  #$opts .= "seqlogo=" . $seqlogo . "\n" if defined($seqlogo) && !exists($cmd_opts{"seqlogo"});
  #$opts .= "avoid=" . $avoid . "\n" if defined($avoid) && !exists($cmd_opts{"xtralen"});
  #$opts .= "low_content=" . $low_content . "\n" if defined($low_content) && !exists($cmd_opts{"low_content"});
  #$opts .= "treg_threshold=" . $treg_threshold . "\n" if defined($treg_threshold) && !exists($cmd_opts{"treg_threshold"});

  return $opts;

}

sub msg_output {
  my $msg_res = <<"OUTPUT";
## Job's done !
  ==========
$timer_stack
  ----------
  Output files have been stored in $directory
  Open the index.html file in your web browser.
  ==========
OUTPUT

  print $msg_res ."\n";
}

sub msg_config {
  my $msg_config = <<"CONFIG";
########################################
  Trawler standalone $tcst{trawler_version}
########################################

## Running $script_name
  ==========
  Reading configuration and pipeline settings...
  ==========
CONFIG

  print $msg_config . "\n";
}

1;
