#!@WHICHPERL@ -Tw
# FILE: metameme.cgi
# AUTHOR: William Stafford Noble
# CREATE DATE: 4-17-98
# PROJECT: MHMM
# DESCRIPTION: Submit a Meta-MEME job to the server and direct you to the output page.
# COPYRIGHT: 2001, Columbia University

use lib qw(@PERLLIBDIR@);
use mhmm_globals;
$ENV{"PATH"} = "";

use CGI;
$cgi = new CGI;
$form =  new CGI;

chdir($outputdir) || die("Can't cd to $outputdir");

print $cgi->header;
print $cgi->start_html;

$output_page = &submit_job($cgi); 

# Direct you to the output page
print "<META HTTP-EQUIV=\"Refresh\" CONTENT=\"1; URL=$output_page\">\n";
print "<BODY bgcolor=\"#FFFFFF\">";

print $cgi->end_html;

############################################################################
# Submit the job
############################################################################
sub submit_job {
  my($formdata) = @_;
  my(%form);
  my($jobnum, $linkpage, $command, $pid, $memefile, $seqfile, 
     @stat, $job_log, $job_log_filename);

  # Create a local hash so that the passed-in parameters can be
  # accessed through the hash.  Note that, in this approach,
  # $form{$key} can only store a scalar value.  If $key's value is an
  # array, $form{$key} stores only the first element.  The alternative
  # approach is to use the $formdata->param() method all the way.
  foreach my $key ($form->param) {
    $form{$key} = $form->param($key);
    #print "$key","\n";
  }

  # Extract the filenames.
  $jobnum = $form{'jobnum'}; #print $jobnum,"\n";
  $job_log_filename = $form{'logfile'}; #print $job_log,"\n";
	# Catch shell metacharacters in the file name
	if ($job_log_filename && $job_log_filename =~ /^([\w\.]+)$/o) {
		$job_log_filename = $1;
	}
  $memefile = '';
  if (defined $form{'memefile'}) {
    $memefile = $form{'memefile'};
  }
  $seqfile = '';
  if (defined $form{'seqfile'}) {
    $seqfile = $form{'seqfile'};
  }

  # Create the output filename.
  $linkpage = "$metasite/output/$jobnum.html";

  # Make sure the MEME file is readable.
	# We are only going to read from this file so we can 
	# allow '/' characters in the file name when we untaint it
	if ($memefile =~ /^([\w\.\/-]+)$/o) {
		$memefile = $1;
	} else {
		&return_error("Invalid file name", "$memefile is not a valid file name");
	}	
  chmod(0664, $memefile);

  # Building the model.
  $command = "$bindir/metameme -bin $bindir ";
  $command .= "-meme $memefile -train $seqfile ";
  $command .= "-hmm mhmm.$jobnum ";
  $command .= "-draw mhmm.$jobnum.gif ";
  $command .= "-topology $form{'topology'} ";
  if (defined $form{'maxseqs'} and ($form{'maxseqs'} ne '')
      and ($form{'maxseqs'} > 0)) {
    $command .= "-maxseqs $form{'maxseqs'} ";
  }

  # $form{'motif'} can only return the first element of the array 
  # in this way of implementation.  To access the whole array, use 
  # the param() method to access the original $formdata.
  if (defined $form{'motif'} && $form{'motif'} ne "")  {
    my @selected_motifs = $form->param('motif');
    foreach my $motif (@selected_motifs) {
      $command .= "-motif $motif ";
    }
  }

  # 
  # Add sequence file description 
  #
  if (defined $form{'Description'} and $form{'Description'} ne ""){
      my $descriptiontext = $form->param('Description');
      $command .= "-description \'$descriptiontext\' ";
  }

  $command .= "-nspacer $form{'nspacer'} ";

  # Finding motif occurrences.
  if (defined $form{'occur'} && $form{'occur'} eq "on") {
      #$command .= "-scoring $form{'scoring'} ";
      #$command .= "-scoring local-viterbi ";
      $command .= "-scoretrain occur.$jobnum ";
  }

  # Searching a database.
  if (defined $form{'search'} && $form{'search'} eq "on") {
    $command .= "-test $datadir/$form{'db'} ";
    #$command .= "-scoring $form{'scoring'} ";
    #$command .= "-scoring local-viterbi ";
    $command .= "-scoretest search.$jobnum ";
    $command .= "-threshold $form{'threshold'} ";
    #$command .= "-minprint $form{'minprint'} ";
    #$command .= "-maxprint $form{'maxprint'} ";
  }

  #dna or protein?
  if (defined $form{'BLOSUM'} && $form{'BLOSUM'} ne ""){
      $command .= "-sc $bindir/$form{'BLOSUM'} ";
  }
  #if ($form{'pamvalue'} ne ""){
  #    $command .= "-pam $form{'pamvalue'} ";
  #}

  $command .= "-web $jobnum $linkpage $webmaster $metasite ";
  #@stat   = $formdata->param('statistic');
  #$command .= "-log $job_log_filename $stat[0] $stat[1] $stat[2] $stat[3] ";
  $command .= "-verbose 1 -notime ";
  $command .= "-statusto $job_log_filename ";

  # Now the server begin to output html file and model graph
  $command .= "-html ";
  #$command .= "-graph ";

  # The exact name of user input
  $command .= "-userinputmotif \'$form{'userinputmotif'}\' ";
  $command .= "-userinputseq \'$form{'userinputseq'}\' ";

	# Detaint, looking for dangerous shell metacharacters in the command 
	if ($command =~ /^([-:\s\w\.\@\/\'\~]+)$/o) {
		$command = $1;
	} else {
    &return_error(
      "Error running metameme, invalid characters in command: $command"
    );
	}

	# Detaint, looking for dangerous shell metacharacters in the file name 
	if ($job_log_filename =~ /^([\w\.\/-]+)$/o) {
		$job_log_filename= $1;
	} else {
    &return_error(
      "Error running metameme, invalid characters in file name:"
     . " $job_log_filename"
    );
	}
  # Pipe output to the logfile.
  $command .= ">> $job_log_filename";

  # Store the command in the logfile.
  open($job_log, ">$job_log_filename") 
    || &return_error("Can't open $job_log_filename.\n");
  print($job_log "$command\n\n");
  close($job_log);
  chmod(0664, $job_log_filename);

  #print $command;

  # Prevent duplicate output after fork by flushing.
  $| = 1;

  # Create a child process just like this one.
  $pid = fork;

  # The child runs the command and then exits.
  if ($pid == 0) {
    # http://www.perl.com/CPAN-local/doc/manual/html/pod/perlfunc/fork.html
    # Note that if your forked child inherits system file descriptors like 
    # STDIN and STDOUT that are actually connected by a pipe or socket, even
    # if you exit, the remote server (such as, say, httpd or rsh) won't think
    # you're done. You should reopen those to /dev/null if it's any issue.

    open(STDIN, "/dev/null");
    open(STDOUT, "/dev/null");
    open(STDERR, "/dev/null");
	
		# run the command
    `$command`;

    exit(0);
  }

  # The parent makes sure that the child started properly.
  elsif (!$pid) {
    &return_error("Error running metameme", "\n$command\n");
  }
  return($linkpage);
}

############################################################################
# Return an error.
############################################################################
sub return_error {
  my($keyword, $message) = @_;

  print $cgi->header;
  print $cgi->start_html(-title=>"Meta-MEME error - $keyword",
                         -target=>'_self');
  print $cgi->h1($keyword),
  $cgi->hr($message),
  $cgi->p("If you have questions, please contact $webmaster.<BR>");
  exit(1);
}

