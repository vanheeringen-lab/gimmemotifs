# Paul Pavlidis 4/2/2000
# DO NOT import microarray specific functions into this library! It is meant to be generic.

package CGIutil;
use CGI qw(:standard);
use CGI::Carp qw(fatalsToBrowser);

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(cleanName dumpParams my_upload outputError sanity sendEmail trace validateParam validateEmail validateFilename validateParam printToFile checkQueryError);

my $SITE_MANAGER = "cegrant\@u.washington.edu";

use strict; 

sub cleanName($$$$) {
  my ($name, $query, $notify, $backto) = @_;
  unless ($name && $name=~/^[\/\\\w\._\-:]+$/) {
    outputError("Illegal characters in file name $name. Make sure it has no \"funny\" characters such as %, \$ or *", $query, $notify, $backto);
  }
  return $name;
}

# Print out the parameters from a cgi query. (debugging tool)
sub dumpParams($) {
  my ($query) = @_;
  my (%params, $item);
  %params = $query->Vars;
  my $string;
  foreach $item (sort keys %params) {
    $string.= "$item $params{$item}\n";
  }
  return $string;
}

# return the name, incl. the subdir within $uploadPath
# $uploadPath is not revealed to the user
sub my_upload(*$$$$$) {
  my ($query, $fieldName, $maxSize, $uploadPath, $mailto, $backto) = @_;
  my $debug = 0;
  my $bufSize = 16384;
  my ($filepath, $fh, $savedFile);
  if ($filepath = $query->param($fieldName)) {
    $filepath = cleanName($filepath, $query, $mailto, $backto);
    if ($filepath =~ /([^\/\\]+)$/) { # strip path to get basename
      $savedFile = "$1";
    } else {
      $savedFile = "$filepath";
    }
    $savedFile =~ s/\s+//g; # remove whitespace from name
    #$savedFile = "$ENV{REMOTE_ADDR}.$savedFile";  
    if (! ($fh = $query->upload($fieldName))) {
      outputError("Failed to obtained uploaded file handle", 
		  $query, $mailto, $backto);
    }
    unless ($debug) { # do not read file if debugging
      unless (-e "$uploadPath/$savedFile") { # let user use existing files. 
	# We'll have to clear it out by hand.???
	if (!open(OUT, ">$uploadPath/$savedFile")) { # failed to open file 
	  outputError("Error opening file for writing: $!", 
		      $query, $mailto, $backto);
	}
	my $fileSize = 0;
	my ($byteRead, $buf);
	while ($byteRead = read($fh, $buf, $bufSize)) {
	  $fileSize += $byteRead;
	  if ($fileSize > $maxSize) {
	    close(OUT);
	    unlink("$uploadPath/$savedFile");
	    outputError ("That file was too large (max = $maxSize bytes)", 
			 $query, $mailto, $backto);
	  }
	  print OUT $buf;
	}
	close(OUT);
	if ((stat "$uploadPath/$savedFile")[7] <= 0) {
	  unlink("$uploadPath/$savedFile");
	  outputError ('Failed to upload file. Make sure you entered the file name correctly and the file is readable', $query, $mailto, $backto);  
	}
      }
    }
  }
  return $savedFile;
} # sub my_upload


sub outputError($$$;$$$) {
  my ($error, $query, $mailto, $backto, $logCallback, $nph) = @_;
#  my ($error, $query, $mailto, $backto) = @_;

  if ($logCallback && ref($logCallback) eq 'CODE') {
    &$logCallback('Error', $error);
  }
  print $query->head();
  print $query->start_html("-title"=>'Error',
			   -BGCOLOR=>'white');
  print $query->h1("There was an error processing your request.");
  print $query->strong($error);
  print $query->p("Try again. If you continue to have problems, ",
		  "or think you have identified a problem with this web site,",
		  " please notify the", a({-href=>"mailto:$mailto"}, 
					  'maintainer'), "of this site.");
  if ($backto) {
    print $query->p(a({-href=>$backto}, 'Back'));
  }
  print $query->end_html;
#  Apache::exit(1);
  if ($ENV{MOD_PERL}) {
  Apache::exit(1);
  } else {
    exit(1); # Cannot call exit under mod_perl.
  }
} # sub sub outputError


# Use this to make a file go to as 'save' window, not the browser.
sub printToFile {
  my ($filename, $query, $SITE_MANAGER, $backto) = @_;
  open (IN, "<$filename") or outputError ("Could not access file $filename.", 
					  $query, $SITE_MANAGER, $backto);
  print $query->header(-type=>'application/octet-stream', -attachment=>"${filename}");
  while (<IN>) {
    print;
  }
  close IN;
} # printToFile 


# Something like this is recommended in the CGI security
# FAQ. http://www.w3.org/Security/Faq/wwwsf4.html What they actually
# suggest is replacing metacharacters; even better is to 'check that
# you are getting exactly what you want'; but that isn't practical
# here, so we do something simple - if it looks bad, we say NO
sub sanity($;$$$$) {
  my ($string, $query, $mailto, $backto, $suppress) = @_;
  return 0 if !defined $string; # no string is okay. This happens when null parameters are tested.
  my $length;
  
  # check length, and look out for * < >, things that begin with '-'.
  if (($length = $string =~ tr/a-zA-Z0-9_-//) > 1000 || 
      $string =~ /(&;\`\'\\\"\|\*\?\~\<\>\^\(\)\[\]\{\}\$\n\r)/o ||
      $string =~ /\s-/
      ) {
    unless ($suppress || !defined $query) { # suppress error message
      outputError("Illegal string: too long or illegal characters", 
		  $query, $mailto, $backto);
    }
    return 0; # error condition
  }
  return $string; # no error
} # sub sanity


sub sendEmail($$$$) {
#  use Mail::Mailer;
#  use Mail::Sendmail;
  my ($from, $address, $subject, $message) = @_;
  my %mail = (To => $address, From => $from, Message => $message, Subject => $subject, smtp => 'genome3.cpmc.columbia.edu');
#  sendmail (%mail) or die "Mail::Sendmail error: $!\n";
  
#  my $mailer = Mail::Mailer->new("sendmail");
#  $mailer->open({From => $from,
#		 To => $address,
#		 Subject => $subject
#		   }) or die "$!\n" ;
#  print $mailer $message;
#  $mailer->close();
  return 0; # no error.
}

sub trace($$) {
  my ($trace, $query) = @_;
  print $query->header;
  print $query->start_html("-title"=>'Debugging',
			   -BGCOLOR=>'white');
  print $query->strong($trace);
  print $query->end_html;
} # trace


# this is recommended in the CGI security FAQ. http://www.w3.org/Security/Faq/wwwsf4.html
sub validateEmail($$$$) {
  my ($email, $query, $mailto, $backto) = @_;
  sanity($email, $query, $mailto, $backto, 0);
  unless ($email =~ m/^[\w.+-]+\@[\w.+-]+$/) {
    outputError ("You did not enter a valid email address $email", $query, $mailto, $backto);
  }
} # sub validateEmail


sub validateFilename($$$$) {
  my ($query, $filename, $mailto, $backto) = @_;
  if ($filename !~ /\w$/) {
    outputError ("You omitted the file name or used illegal characters", $query, $mailto, $backto);
  } elsif ($filename =~ /\*\.\?\$/) {
    outputError ("File names cannot contain most punctuation marks", $query, $mailto, $backto);
  }
  my $filenametemp;
  ($filenametemp = $filename) =~ s/.+[\\\/]//;
  if ($filenametemp =~ /\s/) {
    outputError ("File names $filename cannot contain white space", $query, $mailto, $backto);
  }
} # sub validateFilename

sub validateParam($$$$;$) {
  my ($query, $fieldName, $required, $mailto, $backto) = @_;
  $backto = "" if !defined $backto;
  if ($required) {
    if (!length($query->param($fieldName))) {
      outputError("You omitted a required field: $fieldName",
		  $query, $mailto, $backto);
    }
  }
  my $param = sanity($query->param($fieldName), $query, $mailto, $backto, 0);
  # everything is okay...
  # untaint the param since we already checked it.
  $param =~ /(.*)/;
  return $1;
} # sub validateParam

sub checkQueryError {
  my ($query, $mailto, $backto) = @_;
  my $error = $query->cgi_error;
  if ($error) {
    outputError($error, $query, $mailto, $backto);
  }
}

1;
#package CGIutil;

