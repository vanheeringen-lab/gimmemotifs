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

$CGI::POST_MAX = 1024 * 5000; # 5Mb upload limit

#defaults
my $SAFE_NAME_CHARACTERS = 'a-zA-Z0-9:_\.-';
my $SAFE_NAME = '^[a-zA-Z0-9:_\.][a-zA-Z0-9:_\.-]*$';
my $MOTIFS_DB_INDEX = '@MEME_DIR@/etc/motif_db.index';
my $VERSION = '@S_VERSION@';
my $PROGRAM = 'TOMTOM';
my $HELP_EMAIL = '@contact@';
my $BIN_DIR = '@MEME_DIR@/bin';
my $TOMTOM_OUTPUT_DIR = "../output"; # directory for results
my $TOMTOM_SERVICE_URL = "@OPAL@/TOMTOM_@S_VERSION@";

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
    my ($is_queue_job, $query_name, $query_motif, $query_filter, $target_databases, 
      $uploaded_name, $uploaded_database, $dist, $ttype, $threshold, $address, $v_msg) = check_params($q, $utils);
    if ($utils->has_errors()) {
      $utils->print_tailers();
    } else {
      if ($is_queue_job) {
        submit_to_opal($q, $utils, $query_name, $query_motif, $query_filter, $target_databases,
          $uploaded_name, $uploaded_database, $dist, $ttype, $threshold, $address, $v_msg);
      } else {
        run_on_wserver($q, $utils, $query_name, $query_motif, $query_filter, $target_databases, 
          $uploaded_name, $uploaded_database, $dist, $ttype, $threshold);
      }
    }
  }
}

#
# display_form
#
# Displays the tomtom page.
#
sub display_form {
  my $q = shift;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  my $utils = shift;
  die("Expected MemeWebUtils object") unless ref($utils) eq 'MemeWebUtils';
  
  #there may be inline motifs
  my $inline_motifs = $q->param('inline_motifs');
  my $inline_name = $q->param('inline_name');
  my $inline_count = 0;

  if ($inline_motifs) {
    my ($found, $alpha, $matricies, $columns) = $utils->check_meme_motifs('meme-io', $inline_motifs);
    $inline_count = $matricies;
  }

  #open the template
  my $template = HTML::Template->new(filename => 'tomtom.tmpl');

  #fill in parameters
  $template->param(version => $VERSION);
  $template->param(help_email => $HELP_EMAIL);
  $template->param(inline_motifs => $inline_motifs);
  $template->param(inline_name => $inline_name);
  $template->param(inline_count => $inline_count);

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
# Reads and checks the parameters in preperation for calling tomtom.
#
sub check_params {
  my $q = shift;
  die("CGI required") unless defined $q;
  my $utils = shift;
  die("Expected MemeWebUtils object") unless ref($utils) eq 'MemeWebUtils';
  
  #get the template for the verification message
  my $v_msg = HTML::Template->new(filename => 'tomtom_verify.tmpl');

  #read parameters
  #is this a queued job
  my $is_queue_job = 0; #default to no
  #lookup the meme files for the selected database
  my ($target_name, $target_databases, $uploaded_name, $uploaded_database) = param_database($q, $utils);
  if ($uploaded_name) {
    $is_queue_job = 1;
    $v_msg->param(databases => "Uploaded Database \"$uploaded_name\""); 
  } else {
    $v_msg->param(databases => $target_name); 
  }
  # get the comparison function
  my $dist = param_comparison_function($q, $utils);
  my $dist_name = ($dist eq 'pearson' ? "Pearson's correlation coefficient" : ($dist eq 'ed' ? "Euclidean distance" : "Sandelin-Wasserman similarity"));
  $v_msg->param(distance => $dist_name);
  #get the threshold
  my $threshold = param_num($q, $utils, 'thresh', 0.5, "You must enter a number for the threshold.");
  $utils->whine("The threshold must be larger than zero.") if ($threshold < 0);
  my $ttype = param_thresh_type($q, $utils);
  $v_msg->param(threshold => $threshold, ethresh => ($ttype eq 'evalue'));
  #check the threshold range based on the type
  $utils->whine("The threshold must be smaller than 1 for the q-value") unless ($ttype eq 'evalue' || $threshold < 1);

  #read the motif source
  my $motif_src = $q->param('motif_src'); #leave tainted
  die("Need motif_src param.") unless (defined $motif_src);
  $v_msg->param($motif_src => 1);
  #read the motif into query_motif and the file name into query_name
  my $query_motif;
  my $query_name;
  my $query_filter = '1';
  if ($motif_src eq 'consensus') { 
    # motif is created from a sequence fragement
    my $sequence = $q->param('consensus_sequence'); # taint is removed by conversion
    my $sites = param_num($q, $utils, 'consensus_sites', 20);
    my $pseudo = param_num($q, $utils, 'consensus_pseudocount', 1);
    my $bg = param_background($q, $utils, 'consensus_background_');
    # convert into internal format ignoring parsing errors
    my ($intern) = seq_to_intern($bg, $sequence, $sites, $pseudo);
    # check that the motif has a non-zero length
    if (!$intern) {
      $utils->whine("The IUPAC motif you specified did not contain any IUPAC letter codes.");
    }
    # get the file name from the parsed sequence
    $query_name = motif_id($intern) . ".meme";
    $v_msg->param(sequence => motif_id($intern));
    # create the motif
    my $create_pssm = 0;
    my $create_pspm = 1;
    my $create_header = 1;
    $query_motif = intern_to_meme($intern, $create_pssm, $create_pspm, $create_header);
  } elsif ($motif_src eq 'matrix') {
    # motif is created from a matrix
    my $matrix = $q->param('matrix_data'); # taint is removed by conversion
    die("Matrix required.") unless defined $matrix;
    my $orient = $q->param('matrix_orientation');
    die("Orientation required.") unless defined $orient;
    die("Orientation must be 'row', 'col' or 'auto'.") unless $orient =~ m/^(row|col|auto)$/;
    my $sites = param_num($q, $utils, 'matrix_sites', 20);
    my $pseudo = param_num($q, $utils, 'matrix_pseudocount', 1);
    my $bg = param_background($q, $utils, 'matrix_background_');
    # convert into internal format
    my ($intern, $errors) = matrix_to_intern($bg, $matrix, $orient, $sites, $pseudo);
    if ($intern) { # make the motif
      $v_msg->param(sequence => intern_to_iupac($intern));
      $query_motif = intern_to_meme($intern, 0, 1, 1);
      $query_name = motif_id($intern) . ".meme";
    } else { # display errors
      $utils->whine(@$errors);
    }
  } elsif ($motif_src eq 'meme') {
    # motif is from an uploaded meme file
    my $file = $q->upload('meme_file');
    # copy the file into the query_motif variable
    $query_motif = do { local $/;  <$file> };
    close($file);
    # parse the motif
    my ($found, $alpha, $matricies, $columns) = $utils->check_meme_motifs('meme-io', $query_motif);
    $v_msg->param(count => $matricies, columns => $columns);
    # check the motif
    $utils->whine("Tomtom could not find motifs in query file.") unless $found;
    $utils->whine("Tomtom does not support protein motifs.") if $alpha eq 'PROTEIN';
    # extract the file name and make it safe
    $query_name = param_filename($q, 'meme_file');
    $query_name = "query.meme" unless length($query_name) > 0; 
    $v_msg->param(file => $query_name);
    # get the motif filter
    if ($matricies > 1) {
      my $no_filter = $q->param('meme_ignore_filter');
      if ($no_filter) {
        $query_filter = '';
        $is_queue_job = 1;
        $v_msg->param(filter => 'All motifs');
      } else {
        $query_filter = $q->param('meme_filter');
        die("Need query_filter param.") unless defined $query_filter;
        #check the filter entries
        $query_filter = '1' if $query_filter =~ m/^\s*$/; #make non-empty
        my @entries = split(/\s+/, $query_filter);
        foreach my $entry (@entries) {
          if ($entry !~ m/$SAFE_NAME/) {
            $utils->whine("Filter contains unsupported motif names. Motif names",
                " may not begin with a '-' and may only contain alphanumeric ",
                "characters plus ':', '_', '.' and '-'.");
            last;
          }
        }
        $is_queue_job = 1 if scalar(@entries) > 1;
        $v_msg->param(filter => 'Some motifs (' . join(', ', @entries) . ')');
      }
    }
  } elsif ($motif_src eq 'inline') {
    # motif is from an inline variable
    $query_motif = $q->param('inline_motifs');
    # parse the motif
    my ($found, $alpha, $matricies, $columns) = $utils->check_meme_motifs('meme-io', $query_motif);
    $v_msg->param(count => $matricies, columns => $columns);
    # check the motif
    $utils->whine("Tomtom could not find motifs in inline motifs.") unless $found;
    $utils->whine("Tomtom does not support protein motifs.") if $alpha eq 'PROTEIN';
    # get the file name and make it safe
    $query_name = param_filename($q, 'inline_name');
    $query_name = "query.meme" unless length($query_name) > 0; 
    # remove the filter if more than one motif
    if ($matricies > 1) {
      $query_filter = '';
      $is_queue_job = 1;
    }
  } else {
    die("Illegal value for motif_src param.");
  }
  # get the email address to send the results (if queuing)
  my $address = $q->param('email');
  $utils->check_address($address, $HELP_EMAIL) if ($is_queue_job);
  # add the description to the message
  my $description = $q->param('job_description');
  $description =~ s/[^a-zA-Z0-9:;\-_"\(\)<>%]/ /g;
  $v_msg->param(description => $description);

  return ($is_queue_job, $query_name, $query_motif, $query_filter, $target_databases, 
    $uploaded_name, $uploaded_database, $dist, $ttype, $threshold, $address, $v_msg->output());
}

#
# submit_to_opal
#
# Submits the job to a queue which will probably run on a cluster.
#
sub submit_to_opal {
  my ($q, $utils, $query_name, $query_motif, $query_filter, $target_databases, 
    $uploaded_name, $uploaded_database, $dist, $ttype, $thresh, $address, $v_message) = @_;

  # get a handle to the webservice
  my $service = OpalServices->new(service_url => $TOMTOM_SERVICE_URL);

  # create OPAL requst
  my $req = JobInputType->new();
  my $filter = tag_filter($query_filter);
  my $upload = ($uploaded_name ? "--uploaded $uploaded_name" : '');
  $req->setArgs("$filter $upload $dist $ttype $thresh $query_name $target_databases");

  # create list of OPAL file objects
  my @infilelist = ();
  # 1) Meme file
  push(@infilelist, InputFileType->new($query_name, $query_motif));
  # 1.5) Uploaded database
  push(@infilelist, InputFileType->new($uploaded_name, $uploaded_database)) if $uploaded_name;
  # 2) Email address file (for logging purposes only) - Disabled as could be scraped by web crawler.
  #push(@infilelist, InputFileType->new("address_file", $address));
  # 3) Submit time file (for logging purposes only)
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

#
# run_on_wserver
#
# Runs tomtom on the webserver.
#
sub run_on_wserver {
  my ($q, $utils, $query_name, $query_motif, $query_filter, $target_databases, 
    $uploaded_name, $uploaded_database, $dist, $ttype, $thresh) = @_;

  # create a directory for the results
  my $dir_rel = tempdir('tomtom_XXXXXXXX', DIR => $TOMTOM_OUTPUT_DIR);
  my $dir = abs_path($dir_rel);
  chmod 0777, $dir; # I'm suspicious that the 0777 may be too permissive (but that's what we had before)
  chdir($dir);
  # create a file to view
  my $msg_list = [];
  update_status('index.html', 'Tomtom', 3, [], add_status_msg('Please wait', $msg_list));
  # create the motif file to be searched
  my $query_fh;
  sysopen($query_fh, $query_name, O_CREAT | O_WRONLY, 0777); #too permissive?
  print $query_fh $query_motif;
  close($query_fh);
  # create the database file to be searched
  my $uploaded = "";
  if ($uploaded_name) {
    my $db_fh;
    sysopen($db_fh, $uploaded_name, O_CREAT | O_WRONLY, 0777); #too permissive?
    print $db_fh $uploaded_database;
    close($db_fh);
    $uploaded = "--uploaded $uploaded_name";
  }
  # prepare the pattern to be passed (otherwise it would be auto expanded)
  $target_databases =~ s/\*/\\*/g;
  $target_databases =~ s/\?/\\?/g;
  # prepare the filter
  my $filter = tag_filter($query_filter);
  # send redirect
  send_redirect($q, $dir_rel);
  # call the webservice script directly
  my $out = 'stdout.txt';
  my $err = 'stderr.txt';
  my $id = fork();
  if ($id == 0) {# child process
    # close file descriptors
    close(STDOUT);
    close(STDIN);
    close(STDERR);
    # run webservice
    my $result = system("$BIN_DIR/tomtom_webservice --niced $uploaded $filter $dist $ttype $thresh $query_name $target_databases 1>> $out 2>> $err");
    # check the result for errors
    if ($result != 0) {
      # make a copy of the index in case it had useful messages
      cp('index.html', 'index_old.html');
      # create a new index, describing the problem
      my $file_list = [{file => 'index_old.html', desc => 'Old index'},
          {file => $out, desc => 'Processing messages'},
          {file => $err, desc => 'Error messages'}];
      if ($result == -1) {
        add_status_msg('Tomtom webservice failed to run', $msg_list);
      } elsif ($result & 127) {
        add_status_msg("Tomtom webservice process died with signal ".($result & 127).", ".(($result & 128) ? 'with' : 'without')." coredump", $msg_list);
      } else {
        add_status_msg("Tomtom webservice exited with error code ".($result >> 8), $msg_list);
      }
      update_status('index.html', 'Tomtom', 0, $file_list, $msg_list);
    }
  }
}

#
# sends a redirect request
#
sub send_redirect {
  my ($q, $dir_rel) = @_;
  my $namelen = length($q->url(-relative=>1));
  #turn off buffering
  $| = 1;
  print $q->redirect(substr($q->url(-path_info=>1),-$namelen,-$namelen) .$dir_rel.'/index.html');
}

#
# takes a filter and tags it so it can be passed
#
sub tag_filter {
  my $filter = shift;
  return '' if $filter eq '';
  my $out = '';
  foreach my $part (split(/\s+/, $filter)) {
    if ($part =~ m/^\d+$/) {
      $out .= ' --mi ' . $part;
    } else {
      $out .= ' --m ' . $part;
    }
  }
  return $out;
}

# Getters that pass taint validation

sub param_database {
  my ($q, $utils) = @_;
  my ($target_name, $target_databases, $uploaded_name, $uploaded_database) = ('','','','');
  my $db = $q->param('target_db');
  die("Need target_db param.") unless defined $db;
  if ($db eq 'upload') {
    my $file = $q->upload('db_file');
    #copy the file into the uploaded_database variable
    $uploaded_database = do { local $/;  <$file> };
    close($file);
    $uploaded_name = param_filename($q, 'db_file');
    $uploaded_name = "targets.meme" unless length($uploaded_name) > 0; 
  } else {
    my $index = getnum($db);
    if (defined $index) {
      my @entry = load_entry($MOTIFS_DB_INDEX, $index);
      $target_name = $entry[DB_NAME];
      $target_databases = $entry[DB_ROOT];
    } else {
      $utils->whine("You must select a database to search.");
    }
  }
  return ($target_name, $target_databases, $uploaded_name, $uploaded_database);
}

sub param_background {
  my ($q, $utils, $base) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  die("Expected Utils object") unless ref($utils) eq 'MemeWebUtils';

  my %bg = (A => param_num($q, $utils, $base. 'A', 0.25), C => param_num($q, $utils, $base.'C', 0.25), 
      G => param_num($q, $utils, $base . 'G', 0.25), T => param_num($q, $utils, $base . 'T', 0.25));
  #normalise the background
  my $bg_sum = $bg{'A'} + $bg{'C'} + $bg{'G'} + $bg{'T'};
  $bg{'A'} /= $bg_sum;
  $bg{'C'} /= $bg_sum;
  $bg{'G'} /= $bg_sum;
  $bg{'T'} /= $bg_sum;
  $bg{'dna'} = 1;
  $bg{'source'} = "web form";
  return \%bg;
}

sub param_thresh_type {
  my ($q, $utils) = @_;
  #get the threshold type
  my $thresh_type = $q->param('thresh_type');
  die("Need thresh_type param.") unless defined $thresh_type;
  if ($thresh_type =~ /^(evalue|qvalue)$/) {
    $thresh_type = $1;
  } else {
    die("Thresh type must be qvalue or evalue");
  }
  return $thresh_type;
}

sub param_comparison_function {
  my ($q, $utils) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  die("Expected Utils object") unless ref($utils) eq 'MemeWebUtils';

  #get the comparison function
  my $fun = $q->param('comparison_function');
  die("Missing parameter comparison_function.") unless defined $fun;

  if ($fun =~ /^(pearson|ed|sandelin)$/) {#allowed options
    $fun = $1; #detaint
  } else {
    die("Distance function not an accepted option.");
  }
  return $fun;
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

sub param_num {
  my ($q, $utils, $param_name, $default, $message) = @_;
  die("Expected CGI object") unless ref($q) eq 'CGI';
  die("Expected Utils object") unless ref($utils) eq 'MemeWebUtils';

  my $val = $q->param($param_name);
  if (not defined $val) {
    #only die if given a message (hence this is a post request we're processing)
    die("Error: param_num failed to retrieve an expected parameter '$param_name'.") if $message; 
    return $default;
  }
  my $num = getnum($val);
  if (not defined $num) {
    $num = $default;
    #only whine if given a message
    $utils->whine($message) if $message;
  }
  return $num;
}
