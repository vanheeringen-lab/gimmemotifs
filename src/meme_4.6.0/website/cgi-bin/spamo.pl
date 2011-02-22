#!@WHICHPERL@

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Cwd qw(getcwd abs_path);
use Fcntl;
use File::Basename qw(fileparse);
use File::Copy qw(cp);
use File::Temp qw(tempdir);
use HTML::Template;
use List::Util qw(sum);

use CatList qw(open_entries_iterator load_entry DB_ROOT DB_NAME);
use MotifUtils qw(seq_to_intern matrix_to_intern intern_to_iupac intern_to_meme motif_id motif_width);
use MemeWebUtils qw(getnum valid_meme_version add_status_msg update_status);
use OpalServices;
use OpalTypes;


$CGI::POST_MAX = 1024 * 50000; # 50Mb upload limit

#defaults
my $SAFE_NAME_CHARACTERS = 'a-zA-Z0-9:_\.-';
my $SAFE_NAME = '^[a-zA-Z0-9:_\.][a-zA-Z0-9:_\.-]*$';
my $MOTIFS_DB_INDEX = '@MEME_DIR@/etc/motif_db.index';
my $VERSION = '@S_VERSION@';
my $PROGRAM = 'SPAMO';
my $HELP_EMAIL = '@contact@';
my $BIN_DIR = '@MEME_DIR@/bin';
my $SPAMO_SERVICE_URL = "@OPAL@/SPAMO_@S_VERSION@";

process_request();
exit(0);

#
# Checks for the existance of the cgi parameter search
# and invokes the approprate response
#
sub process_request {
  my $q = CGI->new;
  my $utils = new MemeWebUtils($PROGRAM, $BIN_DIR);
  #decide what to do
  my $search = $q->param('search'); #ignore taint
  
  if (not defined $search) {
    display_form($q, $utils);
  } else {
    #debug
    my ($v_msg, $address, $sequences_name, $sequences_fasta, $primary_name, $primary_motif, 
      $secondary_name, $secondary_databases, $uploaded_name, $uploaded_database, $margin) 
        = check_params($q, $utils);
    if ($utils->has_errors()) {
      $utils->print_tailers();
    } else {
      submit_to_opal($q, $utils, $v_msg, $address, $sequences_name, $sequences_fasta, 
        $primary_name, $primary_motif, $secondary_name, 
        $secondary_databases, $uploaded_name, $uploaded_database, $margin);
    }
  }
}


#
# display_form
#
# Displays the spamo page.
#
sub display_form {
  my $q = shift;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  my $utils = shift;
  die("Expected MemeWebUtils object") unless ref($utils) eq 'MemeWebUtils';
  
  #open the template
  my $template = HTML::Template->new(filename => 'spamo.tmpl');

  #fill in parameters
  $template->param(version => $VERSION);
  $template->param(help_email => $HELP_EMAIL);

  #create the parameter for the database selection
  my $iter = open_entries_iterator($MOTIFS_DB_INDEX, 0); #load the first category
  my @dbs_param = ();
  while ($iter->load_next()) { #auto closes
    my %row_data = (db_id => $iter->get_index(), db_name => ($iter->get_entry())[DB_NAME]);
    push(@dbs_param, \%row_data);
  }
  $template->param(dbs => \@dbs_param);
  
  if ($utils->has_errors()) {
    #finish off the errors
    $utils->print_tailers();
  } else {
    #output the html
    print $q->header, $template->output;
  }
}


#
# check_params
#
# Reads and checks the parameters in preperation for calling spamo
#
sub check_params {
  my $q = shift;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  my $utils = shift;
  die("Expected MemeWebUtils object") unless ref($utils) eq 'MemeWebUtils';
  
  #get the template for the verification message
  my $v_msg = HTML::Template->new(filename => 'spamo_verify.tmpl');
  # get the address
  my $address = $q->param('email');
  $utils->whine("An email address must be specified.") unless $address;
  my $address_confirm = $q->param('email_confirm');
  $utils->whine("An email confirmation address must be specified.") unless $address_confirm;
  if ($address && $address_confirm && $address ne $address_confirm) {
    $utils->whine("The email address does not match the confirmation address");
  }
  # get the sequence file
  my ($sequences_name, $sequences_fasta) = param_file($q, $utils, 'sequences', 'sequences.fasta'); 
  # get the primary motif file
  my ($primary_name, $primary_motif) = param_file($q, $utils, 'primary_motif', 'primary_motif.meme');
  # get the secondary motif file
  my ($secondary_name, $secondary_databases, $uploaded_name, $uploaded_database) = param_secondaries($q, $utils);
  # get the margin

  #my $margin = $q->param('margin');
  #if ($margin < 0) {
  #  $margin = 0;
  #  $utils->whine("Margin must be positive.");
  #}
  my $margin = 150;

  return ($v_msg->output(), $address, $sequences_name, $sequences_fasta, $primary_name, $primary_motif, 
    $secondary_name, $secondary_databases, $uploaded_name, $uploaded_database, $margin);
}

sub safe_name {
  my $name = shift;
  $name =~ tr/ /_/;
  $name =~ s/[^$SAFE_NAME_CHARACTERS]//g;
  #remove leading dashes
  $name =~ s/^-+//g;
  #this is primarily for satisfying taint mode
  if ($name =~ /^([$SAFE_NAME_CHARACTERS]+)$/) { 
    $name = $1;
  } else { #should never be executed
    die("Name contains invalid characters");
  }
  return $name;
}

sub param_filename {
  my ($q, $param_name) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';

  my $upload_path = $q->param($param_name);
  my ($filename, $path) = fileparse($upload_path);
  return safe_name($filename);
}

sub param_file {
  my ($q, $utils, $param_name, $substitute_file_name) = @_;
  my ($data, $name) = ('', '');

  my $file = $q->upload($param_name);
  $data = do { local $/;  <$file> };
  close($file);
  $name = param_filename($q, $param_name);
  $name = $substitute_file_name unless length($name) > 0;
  return ($name, $data);
}

sub param_secondaries {
  my ($q, $utils) = @_;
  my ($secondary_name, $secondary_databases, $uploaded_name, $uploaded_database) = ('','','','');
  my $db = $q->param('secondary_db');
  die("Need secondary_db param.") unless defined $db;
  if ($db eq 'upload') {
    ($uploaded_name, $uploaded_database) = param_file($q, $utils, 'secondary_file', 'uploaded_secondaries.meme');
  } else {
    my $index = getnum($db);
    if (defined $index) {
      my @entry = load_entry($MOTIFS_DB_INDEX, $index);
      $secondary_name = $entry[DB_NAME];
      $secondary_databases = $entry[DB_ROOT];
    } else {
      $utils->whine("You must select a database to search.");
    }
  }
  return ($secondary_name, $secondary_databases, $uploaded_name, $uploaded_database);
}

#
# submit_to_opal
#
# Submits the job to a queue which will probably run on a cluster.
#
sub submit_to_opal {
  my ($q, $utils, $v_message, $address, $sequences_name, $sequences_fasta, $primary_name, 
        $primary_motif, $secondary_name, $secondary_databases, 
        $uploaded_name, $uploaded_database, $margin) = @_;

  # get a handle to the webservice
  my $service = OpalServices->new(service_url => $SPAMO_SERVICE_URL);

  # create OPAL requst
  my $req = JobInputType->new();
  $req->setArgs(($margin ? "--margin $margin " : "").($uploaded_name ? "--uploaded $uploaded_name " : "")."$sequences_name $primary_name $secondary_databases");

  # create list of OPAL file objects
  my @infilelist = ();
  # 1) Sequences file
  push(@infilelist, InputFileType->new($sequences_name, $sequences_fasta));
  # 2) primary motif file
  push(@infilelist, InputFileType->new($primary_name, $primary_motif));
  # 2.5) Uploaded database
  push(@infilelist, InputFileType->new($uploaded_name, $uploaded_database)) if $uploaded_name;
  # 3) Email address file (for logging purposes only) - Disabled as could be scraped by web crawler.
  #push(@infilelist, InputFileType->new("address_file", $address));
  # 4) Submit time file (for logging purposes only)
  push(@infilelist, InputFileType->new("submit_time_file", `date -u '+%d/%m/%y %H:%M:%S'`));
  # Add file objects to request
  $req->setInputFile(@infilelist);

  # output HTTP header
  print $q->header('text/html');

  # Submit the request to OPAL
  my $result = $service->launchJob($req);

  # Give user the verification form and email message
  $utils->verify_opal_job($result, $address, $HELP_EMAIL, $v_message);
  $utils->print_tailers();
}
