#!@WHICHPERL@ -Tw
# FILE: mcast-submit.cgi
# AUTHOR: Timothy Bailey and William Stafford Noble
# CREATE DATE: 1/2002
# PROJECT: MCAST
# DESCRIPTION: Main script for MCAST server.

use lib qw(@PERLLIBDIR@);
use MetaGlobals;
use Submit;
use UID;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use CGIutil;
use strict;
no strict 'subs'; # let us use EAGAIN
$ENV{PATH} = "";
$|++;
$CGI::POST_MAX = 2000000;
main();
exit(0);

####################################################
# MAIN
####################################################
sub main {
  my $query = new CGI;
  my $backto = $query->referer();
  my $redirectto = "$MAIN_SCRIPT";

  # process the form
  if ($query->param()) { 
    my $error = $query->cgi_error;
    if ($error) {
			die "Fatal CGI error: $error";
			#error($error, $query, "", $backto);
    }

    # create user id and directory for their stuff.
    my $uid = getUID($query);

    # clobber any existing log for this id (bug: this can happen...).
    my $log = new Log($uid, 1); 

    # Start logging.
    $log->log("User id is $uid");
    $log->debug("--------------------- $uid ------------------------");
    doCancel($query, $uid, $log, $redirectto); # have we canceled?

    # get the settings.
    my %settings;
    getSettings($query, \%settings, $log); # load our MCAST settings
    $log->debug("Settings done");

    # get and validate the files. We can do this by getting a
    # filehandle using $query->upload, and then reading/writing the
    uploadFiles($query, \%settings, $log, $uid);

    # Send the job off to the analysis script. 
    startJob($query, \%settings, $log, $uid, \&run_mcast);
  } else {
    printform($query); # just print a blank form.
  }
} # sub main

#####################################################################
# getSettings  -  get the parameters from the form
#####################################################################
sub getSettings {
  my ($query, $settings, $log) = @_;

  my @keys = qw(p
		g
		e
		b);

	my(%param_text) = ();
	$param_text{p} = "<em>p</em>-value threshold for motif hits";
	$param_text{g} = "maximum allowed distance between adjacent hits";
	$param_text{e} = "<em>E</em>-values printing threshold";
	$param_text{b} = "pseudocount weight",

  $log->log("Settings:");
  foreach (@keys) {
    my $param = validateParam($query, $_, 0, $SITE_MANAGER);
		if (not is_number($param)) {
    	error("$param is not a valid value for the $param_text{$_},", $query, $log);
		}
    $settings->{$_} = $param;
  }

  # now write the final choices to the log.
  foreach (@keys) {
    $log->log("$_ : $settings->{$_}");
  }

  # Various tests.
  if ($settings->{'p'} < $MINP || $settings->{'p'} > $MAXP) {
    error("<em>p</em>-value threshold must be between $MINP and $MAXP",
	  $query, $log);
  }
  if ($settings->{'g'} < 0 ) {
    error("$param_text{g} must be positive",
	  $query, $log);
  }
  if ($settings->{'e'} < 0 ) {
    error("$param_text{e} must be positive",
	  $query, $log);
  }
  if ($settings->{'b'} < 0 ) {
    error("$param_text{b} must be positive",
	  $query, $log);
  }

} #getSettings

#####################################################################
# uploadFiles
#####################################################################
sub uploadFiles {
  my ($query, $settings, $log, $uid) = @_;

  if (validateParam($query, 'usetest', 0, $SITE_MANAGER)) {
    $log->debug("Using demo files"); # and no upload is needed.
    if (!-e $demoquery || !-e $demodatabase) {
      error("Could not locate the demo files ($demoquery, $demodatabase).", 
	    $query, $log);
    }
    $settings->{qfile} = $demoquery;
    $settings->{datafile} = $demodatabase;
  } else {
    my ($qfile, $datafile);
    $log->debug("Using user's files");
    $qfile = validateParam($query, 'qfile', 1, $SITE_MANAGER);
    $datafile = validateParam($query, 'datafile', 1, $SITE_MANAGER);

    my $fh;
    my ($file, $outfile, $n);
    my @input_files = ($qfile, $datafile);
    my @input_file_types = qw(qfile datafile);

    foreach $file (@input_files) {
      my $type = shift @input_file_types;
      $log->debug("Processing $type $file");
      if (!($fh = $query->upload($type))) {
        error("Could not upload $file", $query, $log);
      }
      $outfile = "$OUTPUT_DIR/${uid}/${uid}.${type}";
      $n = checkfile($fh, $query, $outfile, $log, 1);
      $log->debug("File contains $n lines");
      close $fh;
      $settings->{$type} = $outfile;
    }
  }
  $log->debug("Files ok.");
  return 0;
} # uploadFiles

#####################################################################
# checkfile: check that the file format is valid, and copy it to the
# needed location.
#####################################################################
sub checkfile {
  my ($file, $query, $outfile, $log, $isclassfile) = @_;

  open(OUT, ">$outfile") or error("Could not open $outfile for writing",
				  $query, $log);

  my $line;
  my $m = 0;
  while ($line = <$file>) {
    print OUT $line;
    $m++;
  }
  close OUT;
  chmod(0666, "$outfile");
  return $m;
}

#####################################################################
# Run MCAST
#####################################################################
sub run_mcast {
  my ($settings, $uid, $runlog, $log) = @_;

  # deal with options.
  my $params = "mcast -p-thresh " . $settings->{'p'};
  $params .= " -max-gap " . $settings->{'g'};
  $params .= " -e-thresh " . $settings->{'e'};
  $params .= " -bg-weight " . $settings->{'b'};
  $params .= " -scratch $OUTPUT_DIR/$uid";
  $params .= " " . $settings->{'qfile'};
  $params .= " " . $settings->{'datafile'};

  #Calculate location of platform specific binaries
  my $arch = `/usr/bin/uname -ms|/bin/tr -d ' '`;
  chomp($arch);
  $arch =~ /(\w+)/; # detaint: should only be alphanumeric
  $arch = $1;

  my $mcasthtml = "$OUTPUT_DIR/$uid/${uid}.mcast.html";

  # run MCAST
  $params = "PATH=/bin:/usr/bin:$BINDIR:\$PATH; export PATH;". $params;
  my $status = runProgram($mcasthtml, $runlog, $log, $params);

  return $status;
} #run_mcast

#########################################################################
# Print the form.
#########################################################################
sub printform {
  my ($query) = @_;

    # print the form.
    print($query->header(-nph=>1));

    # title and preamble.
    print($query->start_html("-title"=>"MCAST - Motif Cluster Alignment Search Tool $version - Data submission",
                           -target=>'_self',
                           -bgcolor=>"White",
			     -onload=>"disableFields()"));

    print("<BLOCKQUOTE>",
	  qq(<IMG alt="" src="$titleimage"/><br />
	     <H3>Data Submission Form</H3> <p>
	     MCAST searches a set of sequences (the "database") 
	     using a set of motifs (the "query") for
	     clusters of "hits" to the motifs.
	     This page allows you to submit a query and a database to MCAST.
	     </p>\n));

    print("<p>There is more documentation <A HREF=\"../mcast-intro.html\">available</a>.</p>\n<HR />\n");

    &submit_default($query);
    &submit_advanced($query);

    print("<HR /><p>Send reports of problems to ", 
	  $query->a({href=>"mailto:$SITE_MANAGER"}, "$SITE_MANAGER"), 
	  "</p>\n</BLOCKQUOTE>");
    print($query->end_html);

}

#########################################################################
# Default Options
#########################################################################
sub submit_default  {
  my($query) = @_;

  # Trick Netscape into thinking it's loading a new script:
  my($s) = $query->path_info=~/(\d+)/; # get sequence
  $s++;                                # bump it up

  print($query->start_multipart_form(-name=>'form1',
				     -action=>$query->script_name . "/$s",
				     -onSubmit=>"return validateForm()"),

	"<H3>Inputs</H3>
	 <p>Instead of uploading your own files, you can run MCAST using
            a demonstration query and database by checking this box and
            clicking the 'submit' button: ",
        $query->checkbox(-name=>'usetest', -label=>''),
	"</p>\n",
	"\n

	<OL>
	<li><font size=\"+1\">Query</font><BR />\n",
	$query->filefield(-name=>'qfile', -size=>40, -maxlength=>80),
	"\n<BR />This file should contain a set of motifs.  ",
	$query->a({href=>"../sample.query"}, "Here"), 
	" is an example query.\n<br /><br /></li>

        <li><font size=\"+1\">Database</font><BR />\n",
	$query->filefield(-name=>'datafile',
			  -size=>40,
			  -maxlength=>80),
	"\n<BR />This file should contain the FASTA sequences.  ",
	$query->a({href=>"../sample.database"}, "Here"),
	" is an example.<br /><br /></li>\n
        </OL>");

} # submit_default


#########################################################################
# Advanced Options
#########################################################################
sub submit_advanced {
  my($query) = @_;

  print("<p>",
	"You may click ",
	$query->submit(-label=>'Submit'), 
	" now to use the default options.</p>",

	"<H3>Search options</H3>\n",
	"<p>These options control the search.</p>",
	"<UL>",

	"<p>",
	$query->textfield(-name=>'p', -default=>0.0005, -size=>8, 
			  onBlur=>'validateInput(this)'),
	" <em>p</em>-value threshold for motif hits</li>",

	"<p>",
	$query->textfield(-name=>'g', -default=>200, -size=>5, 
			  onBlur=>'validateInput(this)'),
	" maximum allowed distance between adjacent hits</li>",

	"<p>",
	$query->popup_menu(-name=>'e', "-values"=>[qw/1 10 100 1000/],
			   -default=>'10'), 
	" print matches with <em>E</em>-values less than this value</li>",
		       
	"<p>",
	$query->textfield(-name=>'b', -default=>4, -size=>5, 
			  onBlur=>'validateInput(this)'),
	" pseudocount weight</li>",
	"</UL>",

	$query->submit(-label=>'Submit'),
	$query->end_multipart_form);
}

sub is_number() {
  my $var = $_[0];
	if ($var =~ /^([+-]?)(\d+\.|\.\d­+|\d+)\d*([Ee]([+-]?­\d+))?$/ ) { 
		return(1);
	} else { 
		return(0) 
	}
}
