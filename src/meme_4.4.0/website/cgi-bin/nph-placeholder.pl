#!@WHICHPERL@ -Tw
# Paul Pavlidis

# A (fairly) generic page which users see while waiting for analysis. This page
# refreshes back to the referring script, which returns here if the
# analysis is not done yet. Another required field is the 'referrer',
# which is used to make a 'cancel' link, which returns the user back
# to the form.

use strict;
use lib qw(@PERLLIBDIR@);
use CGI qw(:standard);
use CGIutil;
use CGI::Carp qw(fatalsToBrowser);

use lib qw(.);
use MetaGlobals;
use UID qw(validateUID);
use Log;
$ENV{PATH} = "";
$|++;

main();

sub main {
  my $query = new CGI;
  my $backto = $ENV{HTTP_REFERER}; # go back to the referral page
  if ($query->param()) {
    my $error = $query->cgi_error;
    if ($error) {
      error($error, $query, 0, $backto);
    }
    printPlaceHolder($query, $backto);
  } else {
    error("No parameters", $query, 0, $backto);
  }
}


sub printPlaceHolder {
  my ($query, $backto) = @_;

  # validate that our user is requesting this from the right place
  my $uid = validateParam($query, $UID_FIELD, 1, $SITE_MANAGER, $backto);
  validateUID($query, $uid, $SITE_MANAGER);

  my $caller = validateParam($query, 'CALLER', 1, $SITE_MANAGER, $backto);
  my $referrer = validateParam($query, 'REFERRER', 1, $SITE_MANAGER, $backto);
  my $progress = validateParam($query, 'PROGRESS', 0, $SITE_MANAGER);
  my $process = validateParam($query, 'PROCESS', 0, $SITE_MANAGER);
  my $status = validateParam($query, 'STATUS', 0, $SITE_MANAGER);

  # were we sent here with an error condition?
  if ($progress < 0) {
    print $query->header(-nph=>1);
    print $query->start_html();
    print $query->h2({-style=>"background-color:#EE9999"}, "There was an error during the analysis (status was \"$status\")");
    print "<p>You may <a href=\"${referrer}?$UID_FIELD=$uid&$RECALL_FIELD=FALSE&CANCEL=$process\">return to the form</a> and try again.</p>";
    dumplogs($uid);
    print $query->end_html();

  } else { # still waiting.
    print header(-refresh=>"$REFRESH; URL=$caller?$UID_FIELD=$uid&PROCESS=$process", -nph=>1);
    print start_html();
    print h2("Your analysis is underway");

    $status = "Checking..." unless $status;
    print $query->h2("Progress: $status");
    print $query->p("This page will be updated every $REFRESH seconds until your results are ready.");
    print "<p>You may <a href=\"${referrer}?$UID_FIELD=$uid&$RECALL_FIELD=FALSE&CANCEL=$process\">cancel</a> the analysis and return to the form.</p>";
    print "<p>You may also bookmark this page and return later to check on the status of your analysis.</p>";

    dumplogs($uid);
    print end_html();
  }
}
