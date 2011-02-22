#!@WHICHPERL@
##
## $Id: querystatus.pl 4911 2010-09-03 07:10:33Z james_johnson $
##
## $Log$
## $Rev: 4911 $ $Date: 2010-09-03 00:10:33 -0700 (Fri, 03 Sep 2010) $ $Author: james_johnson $
## uses opal to query job status
##
# Author: Timothy Bailey

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);
use Globals;
use CGI qw/:standard/;          # use the CGI package
use SOAP::Lite;
use OpalServices;
use OpalTypes;

my $service_url = "@OPAL@";
my $service_version = "@S_VERSION@";
my $refresh = 60;
my $activity_url;
if (substr $service_url, -1, 1 eq '/') {
  $activity_url = "$service_url../dashboard?command=statistics";
}
else {
  $activity_url = "$service_url/../dashboard?command=statistics";
}
my $jobid = param('jobid');
my $service = param('service');

$service_url=$service_url."/".$service."_$service_version";

my $meme = OpalServices->new(service_url => $service_url);
my $status = $meme->queryStatus($jobid);

if (eval {$status->fault}) {
  print "<h2>Web Server Error!</h2>";
  print "<p>The following information may help your system administrator ";
  print "solve the problem</p>";
  print "String:".$status->faultstring."<br>";
  print "Code:".$status->faultcode."<br>";
} else { # our soap call didn't fail
  # find out what our status is
  my $resp_code = $status->getCode();
  my $resp_msg = $status->getMessage();
  my $out_url = $status->getBaseURL();
  if ($resp_code==8) { # GramJob.STATUS_DONE
    # Versions of tomcat before 5.5 have a bug that
    # prevents redirection use <meta refresh tag instead.
    #print redirect(-uri=>$out_url,-status=>302);
    printredirect($out_url);
    exit(0);
  } elsif ($resp_code==4) { # GramJob.STATUS_FAILED
    &printheaders;
    print "<h2>Job Failed</h2><br><hr><br>";
    print "The output from your job will be found here: ";
    print "<a href=\"$out_url\">$out_url</a><br>";
    print "<p>The following information may help your system administrator ";
    print "solve the problem</p>";
    print "Message: $resp_msg<br>";
    print "Output code: $resp_code<br>";
  } else {
    &printheaders;
    print "<h2>$resp_msg</h2><br><hr><br>";
    #print "When your job finishes, the output will be found here: ";
    #print "<a href=\"$out_url\">$out_url</a><br>";
    print "This page will contain your job output when it is done.<br>\n";
    print "<br>You can bookmark it for later reference.<br>";
    print "<br>The status of your job will be checked again in $refresh seconds.<br>";
    print "<br><a href=\"$activity_url\">View server activity</a>.<br>";
    #print "<p>The following information may help your system administrator ";
    #print "solve any problems that occur:</p>";
    #print "Message: $resp_msg<br>";
    #print "Output code: $resp_code<br>";
    #print "Message: $resp_msg<br>";
    #print "Base Output URL: <a href=\"$out_url\">$out_url</a><br>";
    #print "Page will refresh in $refresh seconds.<br>";
  }
}

print "
<HR>
</BODY>
</HTML>
";

exit(0);

#
# start the response form
#
sub printheaders {
print <<END; 
Content-type: text/html

<HTML>
<HEAD>
<META HTTP-EQUIV="Refresh" CONTENT="$refresh">
<TITLE> MEME Suite - Query Job Status </TITLE>
</HEAD>
<BODY BACKGROUND=\"../images/bkg.jpg\">
<HR>
END
}

sub printredirect {
my $out_url = shift;
print <<END; 
Content-type: text/html

<HTML>
<HEAD>
<meta http-equiv=\"refresh\" content=\"0;url=$out_url\" />
<TITLE> MEME Suite - Query Job Status </TITLE>
</HEAD>
<BODY BACKGROUND=\"../images/bkg.jpg\">
</BODY>
</HTML>
END
}
