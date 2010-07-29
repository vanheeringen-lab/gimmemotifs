#!@WHICHPERL@ -Tw
# FILE: submit-verify.cgi
# AUTHOR: William Stafford Noble
# CREATE DATE: 5-26-98
# PROJECT: MHMM
# DESCRIPTION: Submit and verify Meta-MEME inputs and gather advanced options.
# COPYRIGHT: 2001, Columbia University

use lib qw(@PERLLIBDIR@);
use mhmm_globals;
$ENV{"PATH"} = "";
chdir($outputdir) || die("Can't cd to $outputdir");

$hostname = `/bin/hostname`;
chomp($hostname);
$ipaddress = $ENV{REMOTE_ADDR};

# Exact name of user input motif file and sequence file
#$userinput_motif = "";
#$userinput_seq = "";

# Copy of the given meme file.
$new_meme_filename = "$inputdir/memefile.$$";

# Copy of the given sequence file.
$new_seq_filename = "$inputdir/seqfile.$$";

# Job log.
$job_log_filename = "$logdir/metameme.log.$$";

# Use the latest version of CGI.pm package
use CGI;
use constant MAX_FILE_SIZE => 1_048_576;
$CGI::POST_MAX = MAX_FILE_SIZE;
$maincgi = new CGI;

# Print the output header.
print $maincgi->header;
print "<BLOCKQUOTE>";

if (!$maincgi->param){

  &print_intro();

} else {
  my(@statistic, @dbase, @all_Motif, $userinputmotif, $userinputseq);

  # Pass the main CGI object and 3 array references (pointers) to the callee.
  ($userinputmotif,$userinputseq,$memetype) 
    = &verify($maincgi, $global_log_filename, $job_log_filename, \@statistic, \@dbase, \@all_Motif);
  if ($userinputseq eq ""){
      &display_statistic("nonstat", $memetype,$maincgi, $dbase[0], $#all_Motif+1, @statistic);
  }
  else{
      &display_statistic("stat", $memetype,$maincgi, $dbase[0], $#all_Motif+1, @statistic);
  }
      
  &prompt_default_option($maincgi, $dbase[1], $job_log_filename,
    $new_meme_filename, $new_seq_filename, "$userinputmotif", 
    "$userinputseq", @statistic);

  &prompt_advanced_option($maincgi, @all_Motif);
}

print $maincgi->a({href=>"../metameme-intro.html"}, "Return"),
    " to the Meta-MEME home page.<BR><BR>",
    "<I>Please send comments and questions to ",
    $maincgi->a({href=>"mailto:$webmaster"},
		"$webmaster"),
    "</I></BLOCKQUOTE>",
    $maincgi->end_html;

##############################################################################
# Print intro
##############################################################################
sub print_intro {
  
  print $maincgi->start_html(-title=>"Meta-MEME v3.3 - Data submission",
                             -style=>{'src'=>"../metameme.css"},
                             -target=>'_self',
                             -bgcolor=>"White",
			     ),
			     #-script=>{-language=>'JAVASCRIPT',
           #          -src=>'../validate.js'}),b
  "<IMG ALIGN='LEFT' SRC='../images/meta-meme.gif' Height='81' Width='444'>",
  "<BR><BR><BR><BR><BR>",
  "<BR><H3>Data submission Form</H3>",
  "Use this form to submit DNA or protein sequences to Meta-MEME. ",
  "In addition to your sequences, Meta-MEME requires as input a ",
  $maincgi->a({href=>"http://meme.sdsc.edu"},"MEME"),
  " motif analysis of your data.  Using the output from MEME, Meta-MEME ",
  "will build a motif-based hidden Markov model of your sequences and ",
  "will use the model to produce a multiple alignment of your sequences ",
  "and to search a sequence database for homologous sequences.",
  #"<P>",
  #"Your data will be processed at the ",
  #$maincgi->a({href=>"http://www.sdsc.edu/"}, 
	#      "San Diego Supercomputer Center"),
  #" and you will be directed to a URL that contains the links to ",
  #"your results.",
  "<HR>",
  "<H3>Inputs</H3>",
  "Please enter the following inputs:<BR>";

  &submit($maincgi);
}

##############################################################################
# Submit
##############################################################################
sub submit {
  my($query) = @_;
  
  print $query->start_multipart_form(-name=>'form1',
				     -onSubmit=>"return validateForm()"),
  "<OL>",
  "<LI>",
  " The name of a file containing a set of ",
  $query->a({href=>"http://meme.sdsc.edu"}, "MEME motifs"),
  " from your sequences:<BR>",
  $query->filefield(-name=>'meme',
		    -size=>40,
		    -maxlength=>80),
  "<BR>",
  "<LI>",
  " The name of a file containing a set of ",
  $query->radio_group(-name=>'seqtype',
		      -"values"=>['DNA','protein'],
		      -default=>'protein',
		      -labels=>{'DNA'=>"DNA",
				'protein'=>"protein"}),
  " sequences in an acceptable ",
  $query->a({href=>"/help_format.html"}, "format"),
  "<BR>",
  $query->filefield(-name=>'seq',
		    -size=>40,
		    -maxlength=>80),
  "<BR>",
  "<LI>",
  " <B>[Optional]</b> Description of your sequences: ",
  "<BR>",
  $query->textfield(-name=>'Description',
                    -default=>'',
                    -size=>'40'),
  "<LI>",
  " <B>[Optional]</B>The sequence database you want to search:<BR>",
  $query->popup_menu(-name=>'database',
		     -"values"=>['',
				 		 'PDB pdb.aa pdb.na',
			       'Non-redundant nr.aa nr.na',
			       'Monthly month.aa month.na',
			       'Yeast yeast.aa yeast.na',
			       'C_elegans c_elegans.aa c_elegans.na',
			       'Ecoli ecoli.aa ecoli.na',
			       'SWISS-PROT swissprot none',
			       'Genpept genpept none',
			       'DBEST none est',
			       'DBSTS none sts',
			       'Vector none vector',
			       'EPD none epd'],
		     -labels=>{'PDB pdb.aa pdb.na'=>"Protein Databank",
			       'Non-redundant nr.aa nr.na'=>"NCBI non-redundant database",
			       'Monthly month.aa month.na'=>"New NCBI entries this month",
			       'Yeast yeast.aa yeast.na'=>"Yeast (S. cerevisiae) genome",
			       'C_elegans c_elegans.aa c_elegans.na'=>"C. elegans genome",
			       'Ecoli ecoli.aa ecoli.na'=>"E. coli genome",
			       'SWISS-PROT swissprot none'=>"SWISS-PROT (protein only)",
			       'Genpept genpept none'=>"Genpept (protein only)",
			       'DBEST none est'=>"NCBI non-redundant EST database (DNA only)",
			       'DBSTS none sts'=>"NCBI non-redundant STS database (DNA only)",
			       'Vector none vector'=>"Vector subset of Genbank (DNA only)",
			       'EPD none epd'=>"Eukaryotic Promoter Database (DNA only)"}),
  "</OL><BR><HR>",

  "<H3>Outputs</H3>",
  "Please select the desired output(s).<BR>",
  "<OL>",
  
  #MHMM model is not optional for output
  "<LI> A motif-based hidden Markov model.<BR>",
  "<LI>",
  $query->checkbox(-name=>'occur', -value=>'on',
		   -label=>"An analysis of the motif occurrences in your sequences.",
		   -onClick=>"changeValue(this)"),
  "<LI>",
  $query->checkbox(-name=>'search', -value=>'on',
		   -label=>"A list of homologous sequences from the selected database.",
		   -onClick=>"changeValue(this)"),
  "<LI>",
  'Maximum number of sequences to return (optional)',
  $query->textfield(-name=>'maxseqs', -default=>'', -size=>'10'),
  "</OL><BR><HR>",

  $query->submit(-label=>'Submit'),
  $query->reset(-name=>'Clear Input'),
  $query->end_form,
  "<HR>";
}

##############################################################################
# Verify
##############################################################################
sub verify {
  # Create a dynamic scope object to be available to the callee.
  local ($form) = $_[0]; 

  my($inline_motifs) = '';
  if (defined $form->param('inline_motifs')) {
    $inline_motifs = $form->param('inline_motifs');
  }
  my($motif_summary) = '';
  if (defined $form->param('motif-summary')) {
    $motif_summary = $form->param('motif-summary');
  }
  my($motif_file_name) = '';
  if (defined $form->param('memename')) {
    $motif_file_name = $form->param('memename');
  }
  my($meme) = '';
  if (defined $form->param('meme')) {
    $meme = $form->param('meme');
  }
  my($seq) = '';
  if (defined $form->param('seq')) {
    $seq = $form->param('seq');
  }
  my($occur) = '';
  if (defined $form->param('occur')) {
    $occur = $form->param('occur');
  }
  my($search) = '';
  if (defined $form->param('search')) {
    $search = $form->param('search');
  }
  my($maxseqs) = '';
  if (defined $form->param('maxseqs')) {
    $maxseqs = $form->param('maxseqs');
  }
  my($database) = '';
  if (defined $form->param('database')) {
    $database = $form->param('database');
  }
  my($seqtype) = '';
  if (defined $form->param('seqtype')) {
    $seqtype = $form->param('seqtype');
  }

  my($motif_summary_begin) 
    = "<INPUT TYPE = HIDDEN NAME = motif-summary VALUE = \"";
  my($motif_summary_end) = "\">";
  my($global_log_filename, $job_log_filename, $stat, $dbse, $all_mtf)
    = @_[1..5];

  my($getsize1, $getsize2, $x);
  my(@stat, @db_list, $all_motifs_str);

  # Use these lexical scoping variables to store data that will be returned
  # to the caller.
  my(@db_items, @all_motifs);

  # Make sure we got all the inputs.
  if ($meme eq '') {
    # No MEME file
    &return_error(
      $global_log_filename, 
      $job_log_filename, 
      "No MEME file provided.",
      "You must provide a MEME analysis of your sequences."
    );
    exit (1);
  }
  else {
    $userinput_motif = $meme;
    $userinput_motif =~ s/.*[\/\\]([\w\s]*)/$1/;
  }
  if ($seq eq '') {
    # No sequence file
    $userinput_seq = '';
    if ($occur eq 'on'){
      &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "No sequence file provided.",
        "If you select analysis of the motif occurrences as output"
        .  " you must provide a sequence file."
      );
    }
  }
  else {
    $userinput_seq = $seq;
    # We only what the file name, not the path
    $userinput_seq =~ s/.*[\/\\]([\w\s]*)/$1/;
  }
  if ((($occur eq 'on') or ($search eq 'on')) and (not $maxseqs eq '')) {
    # Check that maxseqs field is a nummber
    if (($maxseqs + 0) <= 0) {
      # maxseqs not a positive number
      &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Maximum number of sequences not a positive number.", 
        "The maximum number of sequences must be a positive number or empty."
      );
    }
  }

  if ($search eq 'on') {
    # If search is turned on the must have selected a db.
    if ($database eq '') {
      &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "No database selected.",
        "If you request a list of homologous sequences from a database as "
        . " output you must select a database."
      );
    }
    # Choose the appropriate database.
    # The actual name of the database file is in positions 1 (peptide)
    # or 2 (nucleotide) of $database.
    # If no such type exists, that name is "none".
    @db_list = split(/ +/, $database);
    $db_items[0] = $db_list[0];
   
    if ($seqtype eq "DNA") {
        $db_items[1] = $db_list[2];
    } elsif ($seqtype eq "protein") {
        $db_items[1] = $db_list[1];
    }
    # Make sure the right flavor of database exists
    if ($db_items[1] && $db_items[1] eq "none") {
      &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Database and sequence type mismatch",
        "There is no $seqtype version of $db_list[0]."
      );
    }
  }

  # Copy the uploaded files to local files to be used by the
  # metameme.cgi script later on because CGI.pm automatically
  # unlinks them after submit-verify.cgi finishes executing. 
  # In the process, remove extra carriage returns.
  my($meme_type);
  open ($new_meme_file, ">$new_meme_filename");
  if ($inline_motifs){
    #my $meme_file = $form->tmpFileName($form->param('inline_motifs'));
    #my $meme = `cat $meme_file`;
    $sub="";
    $inline_motifs =~ s/\r/$sub/g;
    $motif_summary =~ s/\r/$sub/g;
      
    print ($new_meme_file $inline_motifs);
    print ($new_meme_file $motif_summary_begin);
    print ($new_meme_file $motif_summary);
    print ($new_meme_file "$motif_summary_end\n");

    $meme_type = $seqtype;
  }
  else {
    my $meme_file = $form->tmpFileName($meme);
		# We are only going to read from this file so we can 
		# allow '/' characters in the file name when we untaint it
		if ($meme_file =~ /^([\w\.\/]+)$/o) {
			$meme_file = $1;
		} 
    else {
			&return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Invalid file name",
				"$meme_file is not a valid file name"
      );
		}	
    my $meme = `/bin/cat $meme_file`;
    #s/\r\n/\n/g;
    #s/\r/\n/g;
    $sub = "";
    $meme =~ s/\r/$sub/g;
    $meme_type = type_check($meme);
      
    if ($seq ne "") {
	  	if ($meme_type ne $seqtype) {
	      &return_error(
          $global_log_filename, 
          $job_log_filename, 
          "Type does not match",
			    "The sequence file you provided is in $seqtype format"
          . " but your motif file is in $meme_type format.<BR>"
        );
	  	}
    }
    $meme_version = version_check($meme);

    if ($meme_version eq "old"){
	  	$temp_meme_file = "temp_meme_file.$$";
	  	$temp_meme_error = "temp_meme_error.$$";
	  	$status = system("$bindir/oldmeme2meme $meme_file \\
                    1> $temp_meme_file 2>$temp_meme_error");
	  	$error = `/bin/cat $temp_meme_error`;
	  	unlink ($temp_meme_error);
	  	if ($status){
	      &return_error(
          $global_log_filename, 
          $job_log_filename, 
          "oldmeme2meme error",
			    "An error occurred when oldmeme2meme program"
          . " attemped to convert your MEME file to new version."
          . " <BR> oldmeme2meme returned: $error"
        );
	    }
   
	    $meme = `/bin/cat $temp_meme_file`;
	    $meme =~ s/\r/$sub/g;
      unlink ($temp_meme_file);
    }
      
    print ($new_meme_file $meme);
  }
  close($new_meme_file);
  chmod(0664, $new_meme_filename);

  if ($seq ne ""){
    my $seq_file = $form->tmpFileName($seq);
		# We are only going to read from this file so we can 
		# allow '/' characters in the file name when we untaint it
		if ($seq_file =~ /^([\w\.\/]+)$/o) {
			$seq_file = $1;
		} 
    else {
			&return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Invalid file name",
				"$seq_file is not a valid file name"
      );
		}	
    my $sequences = `/bin/cat $seq_file`;
    if (open($new_seq_file, ">$new_seq_filename")) {
	  	$sub = "";
	  	$sequences =~ s/\r/$sub/g;
	  	print($new_seq_file $sequences);
	  	close($new_seq_file);
    }
    chmod(0664, $new_seq_filename);
      
    # Convert the data to FASTA format.
    my $readseq_file = "readseq.$$";
    my $readseq_error = "readseq.error.$$";
    $status = system("$bindir/readseq -a -f=8 $new_seq_filename \\
                    1> $readseq_file 2> $readseq_error");
    $error = `/bin/cat $readseq_error`;
    unlink($readseq_file, $readseq_error);
    if ($status) {
	  	&return_error(
        $global_log_filename, 
        $job_log_filename, 
        "READSEQ error", 
			  "An error occurred when the READSEQ program attempted"
        . "to convert your dataset to FASTA format.<BR>"
        . "READSEQ returned: $error"
      );
    }
      
    # Run the 'getsize' program to get information on converted data.
    my $getsize_file = "getsize.$$";
    my $getsize_error = "getsize.error.$$";
    my $status = system("$bindir/getsize $new_seq_filename \\
                       1> $getsize_file 2> $getsize_error");
    my $error = `/bin/cat $getsize_error`;
    chop($getsize = `/bin/cat $getsize_file`); 
    unlink($getsize_file, $getsize_error);
    if ($error || $status || ($getsize eq "")) {
	    &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Sequence error",
			  "After converting to FASTA format using the READSEQ"
        . "program, the following errors in your dataset were"
        . "detected:<BR>$error.\n"
      );
    }
    ($num, $min, $max, $ave, $total, $letters) = split (' ', $getsize);
    # Check for problem reading the dataset.
    if ($num <= 0) {
	    &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Sequence format error",
			  "Your dataset appears to be in a format that Meta-MEME"
        . " does not recognize.<BR> Please check to be sure that"
        . " your data is <A HREF=\"/help_format.html\">"
        . " formatted</A> properly."
      );
    }
      
    # Make sure there isn't too much data.
    if ($total > $maxdataset) {
	    &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Data set too large.", 
			  "The data you have entered contains more than $maxdataset"
			  . " characters.  Meta-MEME cannot process it at this time. "
			  . "<BR>Please submit a smaller dataset."
      );
    }
      
    # Check for bad sequences.
    if ($num > 0 && $min == 0) {
	    &return_error(
        $global_log_filename, 
        $job_log_filename, 
        "Zero-length sequence.",
			  "Your dataset appears to contain one or more zero-length"
        . " sequences.<BR> Please check to be sure that your data"
        . " is <A HREF=\"/help_format.html\"> formatted</A> properly."
      );
    }
      
    # Check for illegal characters.
		if ($letters) {
    	$x = $letters;
      $x =~ tr/ABCDEFGHIKLMNPQRSTUVWXYZ//cd;
      if (length($x) != length($letters)) {
	      $x = $letters;
	      $x =~ tr/ABCDEFGHIKLMNPQRSTUVWXYZ//d;
	      &return_error(
          $global_log_filename, 
          $job_log_filename, 
          "Unrecognized letters",
			    "Your sequences contain the following unrecognized"
          . " letters: $x<BR> Please convert your sequences to one"
          . " of the sequence <A HREF=/help_alphabet.html>"
          . " alphabets</A> that Meta-MEME recognizes."
        );
			}
    }
  }

  # Find motifs in the MEME file.
  my $motif_ids_file = "motifs.$$";
  my $motif_ids_error = "motifs.error.$$";
  $status = system(
    "$bindir/meme-io --indices --verbosity 1  $new_meme_filename"
    . " 1> $motif_ids_file 2> $motif_ids_error"
  );
  $error = `/bin/cat $motif_ids_error`;
  unlink($motif_ids_error);
  if ( $status) {
    &return_error(
      $global_log_filename, 
      $job_log_filename, 
      "MEME file error",
      "Meta-MEME could not parse the MEME file that you provided."
      . " <BR>$error status = $status"
      . " <BR>Meta-MEME requires the HTML version of the MEME file."
    );
  }

  $all_motifs_str = `/bin/cat $motif_ids_file`;
  unlink($motif_ids_file);
  chop($all_motifs_str);
  @all_motifs = split(' ', $all_motifs_str);

  # Make sure we found some motifs.
  if ($#all_motifs < 0) {
    &return_error(
      $global_log_filename, 
      $job_log_filename, 
      "Motifs not found",
      "Meta-MEME could not find any motifs in the MEME file"
      . " that you provided."
    );
  }

  @{$stat}    = ($num, $min, $max, $ave, $total);
  @{$dbse}    = @db_items;
  @{$all_mtf} = @all_motifs;
  
  return($userinput_motif,$userinput_seq,$meme_type);
}


############################################################################
##### REPORT JOB CHARACTERISTICS
############################################################################

sub display_statistic{
  my($stat,$memetype,$entry, $db_name, $num_motifs, @stat) = @_;
  my($num, $min, $max, $ave, $total) = @stat;

  # The entry parameters have to be copied to local variables
  # because the 'print' statement, for some reason, won't
  # work directly with the $entry->param() method.
  my $seqtype = $entry->param('seqtype');
  my $mult    = $entry->param('mult');
  my $occur   = $entry->param('occur');
  my $search  = $entry->param('search');
  my $maxseqs  = $entry->param('maxseqs');

  print $entry->start_html(-title=>"Meta-MEME v3.3 - Data Verification",
                           -style=>{"Src"=>"../metameme.css"},
                           -target=>'_self',
                           -bgcolor=>"White"),
  "<IMG ALIGN='LEFT' SRC='../images/meta-meme.gif' Height='81' Width='444'>",
  "<BR><BR><BR><BR><BR>";
  #"You have submitted Meta-MEME job \#$$."
 
  if ($stat eq "stat"){
      # Print characteristics of the data set.
      print $entry->p("<H3>Statistics of your dataset:</H3>
     <ul>
     <TABLE BORDER>
      <TR> <TD> Type of sequences <TH ALIGN=RIGHT> $seqtype
      <TR> <TD> Number of sequences <TH ALIGN=RIGHT> $num
      <TR> <TD> Shortest sequence (residues) <TH ALIGN=RIGHT> $min 
      <TR> <TD> Longest sequence (residues) <TH ALIGN=RIGHT> $max 
      <TR> <TD> Average sequence length (residues) <TH ALIGN=RIGHT> $ave 
      <TR> <TD> Total dataset size (residues) <TH ALIGN=RIGHT> $total 
    </TABLE>");
      print "</ul>\n";
  }

  print "The MEME output file you supplied contains $num_motifs motifs.<BR>";

  # Tell the user what files to expect.
  print $entry->p("The following output will be produced:"),
  "<UL>",
  "<LI> A motif-based hidden Markov model.</LI>\n";
  if ($occur && ($occur eq "on")) {
    if (($maxseqs ne '') and ($maxseqs > 0)) {
      print "<LI> An analysis of the motif occurrences in your ",
      "sequences. A maximum of $maxseqs sequences will be returned.</LI>\n";
    }
    else {
      print "<LI> An analysis of the motif occurrences in your ",
      "sequences.</LI>\n";
    }
  }
  if ($search && ($search eq "on")) {
    if (($maxseqs ne "") and ($maxseqs > 0)) {
      print "<LI> A list of homologous sequences from the ",
      "<B>$db_name</B> database. A maximum of $maxseqs sequences will be returned.</LI>\n";
    }
    else {
      print "<LI> A list of homologous sequences from the ",
      "<B>$db_name</B> database.</LI>\n";
    }
  }
  print "</UL>\n",
  $entry->end_html;
}


##############################################################################
# Prompt default option
##############################################################################
sub prompt_default_option{
  my($entry, $db_filename, $job_log_filename, $new_meme_filename,
    $new_seq_filename, $userinput_motif, $userinput_seq, $maxseqs, @stat)
      = @_;

  #my $maxseqs = $entry->param('maxseqs');
  print $entry->start_multipart_form(-name=>'form2',
                                     -action=>'metameme.cgi');

  print "By default, Meta-MEME uses reasonable parameter settings.
         You may click the <B>Submit</B> button now to accept these defaults,
         or scroll down to select advanced options.<BR><BR>";
  print $entry->submit(-label=>'Submit');

  print "<HR>";
  print "The following options are grouped according to their relevance to the
         various stages in the Meta-MEME process. Click ",
      $entry->a({href=>"/doc/mhmm.html"}, "here"),
      " for an explanation of the options.";
  print "<HR>";

  # Pass extra parameters to metameme.cgi.
  print $entry->hidden(-name=>'jobnum',
                       -default=>$$);
  print $entry->hidden(-name=>'db',
                       -default=>$db_filename);
  print $entry->hidden(-name=>'logfile',
		       -default=>$job_log_filename);
  print $entry->hidden(-name=>'memefile',
		       -default=>$new_meme_filename);
  print $entry->hidden(-name=>'seqfile',
		       -default=>$new_seq_filename);
  print $entry->hidden(-name=>'statistic',
                       -default=>[@stat]);
  print $entry->hidden(-name=>'userinputmotif',
		       -default=>$userinput_motif);
  print $entry->hidden(-name=>'userinputseq',
		       -default=>$userinput_seq);
  print $entry->hidden(-name=>'maxseqs',
		       -default=>$maxseqs);
  

  # At this point, do not end multipart form yet.
  # Wait until it finishes the advanced-option form.
}


##############################################################################
# Prompt advanced option
##############################################################################
sub prompt_advanced_option{
  my($form, @motifs) = @_;

  print "<H3>Model building</H3>\n";
  print "<UL>\n";
  print "<BR><LI><B>Model topology</B><BR><BR>\n",
  $form->radio_group(-name=>'topology',
		     -"values"=>['linear','complete'],
		     -default=>'linear',
		     -linebreak=>'true',
		     -labels=>{'linear'=>"Build a model with a linear topology.",
			       'complete'=>"Build a model with a complete connected topology."});

  my $num_motifs = $#motifs + 1;
  my @value = (1..$num_motifs);
  print "<BR><LI><B>Motifs</B><BR><BR>\n",
  "Include these motifs in the model: ",
  $form->checkbox_group(-name=>'motif',
			-"values"=>[@value],
			-defaults=>[@motifs]);
  
  print "\n<BR><BR><LI><B>Spacer models</B><BR><BR>\nUse ",
  $form->popup_menu(-name=>'nspacer',
		    -"values"=>[1,2,3,4],
		    -default=>1),
  " state(s) to model each spacer.";
  print "</ul>";

  if (($form->param('search') && ($form->param('search') eq "on"))) {
    print "<HR>";
    print "<H3>Homology detection</H3>\n";
    print "<UL>\n";

    #  
    # Add E-value option
    #
    print "<BR><LI><B>E-value threshold</B><BR><BR>\n",
    "Report sequences with E-value below ",
    $form->textfield(-name=>'threshold', - default=>10, -size=>4),
    ".<BR>";
    #maxprint option is on longer availabe
    #print "<BR><LI><B>Maximum number of sequences</B><BR><BR>\n",
    #"Report scores for at most ",
    #$form->textfield(-name=>'maxprint', -default=>100, -size=>3),
    #" sequences.<BR>";
    if ($memetype =~ /protein/){
	    print "<BR><LI><B>BLOSUM</B><BR><BR>\n",
	       $form->radio_group(-name=>'BLOSUM',
			   -"values"=>['BLOSUM62','BLOSUM45','BLOSUM80'],
			   -default=>'BLOSUM62',
			   -linebreak=>'true',
			   -labels=>{'BLOSUM62'=>"BLOSUM62",
				     'BLOSUM45'=>"BLOSUM45",
				     'BLOSUM80'=>"BLOSUM80"});
    }
    if ($memetype =~ /DNA/){
	    print "<BR><LI><B>PAM</B><BR><BR>\n",
	      "PAM value ",$form->textfield(-name=>'pamvalue',-default=>1,-size=>4),"<BR>";
    }
	  print "</UL>\n";
  }
  print "<HR>\n";

  print $form->submit(-label=>'Submit');

  # Pass along all the options and entries from the submission page
  # and this advanced-option page to the metameme.cgi script
  foreach my $key ($form->param) {
    print $form->hidden(-name=>$key,
                        -default=>[$form->param($key)]);
  }

  print $form->end_multipart_form;

  print "<HR>\n";
} 


###############################################################################
# Return error
###############################################################################
sub return_error {
  my ($global_log_filename, $job_log_filename, $keyword, $message) = @_;
  my ($field, $date_and_time, $global_log, $job_log);

  # Print an error log for this job.
  open($job_log, ">$job_log_filename");
  print($job_log "keyword=$keyword\n");
  print($job_log "message=$message\n");
	#TODO: mult field not currently supported
  #foreach $field ('seq', 'meme', 'database', 'description', 'mult',
  foreach $field ('seq', 'meme', 'database', 'description', 'occur', 'search') {
    my $value = $form->param($field);
		if (!defined $value) {
			$value = "";
		}
    print($job_log "$field=$value\n");
  }
  close($job_log);
  chmod(0664, $job_log_filename);

  # Add an entry to the global log file.
  #chop($date_and_time = `date`);
  $hostname = `/bin/hostname`;
  chomp($hostname);
  $ipaddress = $ENV{REMOTE_ADDR};
  ($date,$time) = &date_and_time();
  open($global_log, ">>$global_log_filename");
  #print($global_log "$date_and_time \t $$ \t <1> ");
  #foreach $field ('seq', 'meme', 'database', 'mult', 'occur', 
	#	  'search') {
    #my $value = $form->param($field);
    #print($global_log "$field=$value; ");
  #}

  print ($global_log "$hostname $$ submit: $date $time start: $date $time end: $date $time <1> $ipaddress "); 
  print($global_log "\n\n");
  close($global_log);
  chmod(0664, $global_log_filename);

  # 'form' has dynamic scope from the caller
  print $form->start_html(-title=>"Meta-MEME error - $keyword",
                          -target=>'_self');
  print $form->h1($keyword),
  $form->hr($message),
  $form->p("If you have questions, please contact ",
	   "<A HREF=\"mailto:$webmaster\">$webmaster</A>.<BR>");
  print "</BLOCKQUOTE>\n";
  print $form->end_html;
  exit(1);
}

###############################################################################
# Decide the type of motif file -- DNA or Protein?
###############################################################################
sub type_check {
  my ($file) = $_[0];
  my (@temp) = split('\n',$file);    
  my ($i,@type_info);
    
  for ($i=0;$i<=$#temp;$i++){
       
	  $temp[$i] =~ /(.+)/;
		if ($temp[$i] =~ /ALPHABET=/){
   		@type_info = split(/ALPHABET= /,$temp[$i]);
   		if (length($type_info[1]) >= 20){
				return "protein";
   		} else {
				return "DNA";
   		}
   		last;
		}
  }
}

#############################################################################
# Decide the version of MEME file --  old or new?
#############################################################################
sub version_check{
    my ($file) = $_[0];
    my $check = "NAME = bgfreq";
    my (@temp) = split('\n',$file);
    my ($i,$mark);
    for ($i=0;$i<=$#temp;$i++){
	if ($temp[$i] =~ /$check/){
	    #return "new";
	    #break;
	    $mark = 1;
	}
    }
    #return "old";
    if (defined $mark && $mark == 1){
	return "new";
    }
    else{
	return "old";
    }
}
	  
############################################################################
# Return Current Date and Time
# This is used to keep LOG the same format as MEME LOG file
############################################################################
sub date_and_time{
    my ($sec,$min,$hour,$mday,$mon,$year,@rest) = localtime(time);
    # I do not understand why we need to add 1 to month
    $date_now = sprintf "%02d/%02d/%02d",$mon+1,$mday,$year-100;
    $time_now = sprintf "%02d:%02d:%02d", $hour,$min,$sec;
    return($date_now,$time_now);
}
