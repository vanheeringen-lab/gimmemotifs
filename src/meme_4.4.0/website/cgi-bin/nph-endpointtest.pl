#!@WHICHPERL@ -Tw
# Paul Pavlidis

use strict;
use lib qw(@PERLLIBDIR@);
use CGI qw(:standard);
use CGIutil;
use CGI::Carp qw(fatalsToBrowser);

use lib qw(.);
use MetaGlobals;
use UID qw(validateUID);
use Process;
use Log;

$ENV{PATH} = "";
$|++;

&main();
#----------------------
# Main
#----------------------
sub main(;) {
  my $query = new CGI();
  my $backto = $query->referer() || $MAIN_SCRIPT; 	# go back to the referral page
  my $redirectto = "$SITE_URL/cgi-bin/nph-placeholder.cgi";

  if ($query->param()) { # GET method
    my $error = $query->cgi_error;
    if ($error) {
      outputError($error, $query, $SITE_MANAGER, $backto);
    }

    # validate that our user is requesting this from the right place
    my $uid = validateParam($query, $UID_FIELD, 1, $SITE_MANAGER, $backto);
    validateUID($query, $uid, $SITE_MANAGER);

    # check to see if the process we're waiting for is still alive. If
    # it is dead, it could be because it is done, or because it
    # suffered an error.
    my $process = validateParam($query, 'PROCESS', 0, $SITE_MANAGER);
    my $sema = readSemaphore($uid);
    if ($sema eq "") {
      my $log = new Log($uid, 0);
      error("No signal $UPLOAD_DIR/$uid/${uid}.$SEMAPHORE_SUFFIX !", $query, $log);
    }
    my $status = isDone($process, $uid, $sema);

    # Do what needs to be done based on the status.
    if ($status eq 'finished') {
      print $query->header(-nph=>1);
      print "<html><head><title>Finished</title></head><body><h1 style=\"background-color:#55AA55\">Finished!</h1>\n";

      #my $results = `/usr/bin/cat $OUTPUT_DIR/$uid/${uid}.mcast.html`;
      #print "<pre>$results</pre>\n";
      print "<h3><a href=\"../output/$uid/${uid}.mcast.html\">View results</h3></a>\n";
      print "<br><a href=\"$backto\">Return to the input form to do another run.</a><br>\n";

      dumplogs($uid);
    } elsif ($status eq 'died' || $status eq 'aborted') {
      print redirect(-location=>"$redirectto?$UID_FIELD=$uid&CALLER=nph-endpointtest.cgi&REFERRER=$backto&PROGRESS=-1&PROCESS=$process&STATUS=$status", -nph=>1);
    } elsif ($status eq 'running' || $status eq 'waiting') {
      print redirect(-location=>"$redirectto?$UID_FIELD=$uid&CALLER=nph-endpointtest.cgi&REFERRER=$backto&PROGRESS=1&PROCESS=$process&STATUS=$status", -nph=>1);
    } else {
      outputError("Insane Status $status returned", $query, $SITE_MANAGER);
    }
  } # endif param
  else {
      outputError("No parameters", $query, $SITE_MANAGER);
  }
} # main
