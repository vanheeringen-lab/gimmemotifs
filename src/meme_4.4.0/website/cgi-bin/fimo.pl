#!@WHICHPERL@

# Author: Timothy Bailey

# defaults
$PROGRAM = "FIMO";
$NERRORS = 0;				 # no errors yet
$MAX_UPLOAD_SIZE = 1000000;

use lib qw(@PERLLIBDIR@);
require "Utils.pm";			# must come before "whines"
eval { require CGI; import CGI qw/:standard/; }; whine("$@") if $@;
eval { require XML::Simple; import XML::Simple;}; whine("$@") if $@;
use Globals;
use Validation;
use OpalServices;
use OpalTypes;

# get the directories using the new installation scheme
$dir = "@MEME_DIR@";	# installed directory
$logs = "$dir/LOGS";	# directory for logs
$bin = "$dir/bin";	# binary directory
$service_url = "@OPAL@/FIMO_@S_VERSION@";
$email_contact = '@contact@';
$names_file = "fasta_db.csv";
$short_dna_only = 0;		# FIMO can scan long DNA sequences
$dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file&short_only=$short_dna_only";
$db_doc_url = $dbnames_url . "&doc=1";

# get the parameters for the query
get_params();

# if there is no action specified, print an input form
if (! $action) {
  print_form() unless $NERRORS;
} else {
  check_params();
  submit_opal() unless $NERRORS;
  print_tailers();
}

##############################
# SUBROUTINES 		     #
##############################

#
# print the input form
#
sub print_form {
  my $action = "fimo.cgi";
  my $logo = "../doc/images/fimo_logo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = "Use this form to submit motifs to <B>$PROGRAM</B> to be used in searching a sequence database.";

  #
  # required section
  #

  # required left side: address, motifs
  my $motif_doc = "../doc/meme-format.html";
  my $motif_example = "../examples/sample-dna-motif.meme-io";
  my $example_type = "DNA minimal motif format file";
  my $inline_descr = "MEME motifs from sequences in file <TT>'$inline_name'</TT>.";
  my $req_left = make_address_field($address, $address_verify);
  $req_left .= "<BR>\n";
  $req_left .= make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet,
    $motif_doc, $motif_example, $inline_descr, $example_type);

  # required right side: database
  my $fasta_doc = "../doc/fasta-format.html";
  my $sample_db = "../examples/sample-kabat.seq";
  my $req_right = make_supported_databases_field("database", $dbnames_url, $db_doc_url);
  $req_right .= "<BR>\n<B>or</B>\n<BR>\n";
  $req_right .= make_upload_fasta_field("upload_db", $MAX_UPLOAD_SIZE, $fasta_doc, $sample_db);

  my $required = make_input_table("Required", $req_left, $req_right);

  # 
  # optional arguments
  #

  # optional left: description
  my $descr = "motifs";
  my $opt_left = make_description_field($descr, $description);


  # optional right side:
  my @pv_list = ("1", "0.1", "0.01", "0.001", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8", "1e-9");
  my $pv_selected = "1e-4";
  my $opt_right = make_select_field('<i>p</i>-value output theshold<br/>', 'output_pthresh', 
      $pv_selected, @pv_list);
  $opt_right .= "<BR><BR>\n";
  #my $text = "Scan given <A HREF='$doc' <B>strand</B></A> only (DNA motifs only)";
  my $text = "Scan given strand only";
  $opt_right .= make_checkbox("norc", "1", $text, 0);
  $opt_right .= "<BR>\n";

  my $optional = make_input_table("Optional", $opt_left, '<td style="width: 60%;">' .
      $opt_right . '</br></td>');

  #
  # print final form
  #
  my $form = make_submission_form(
    make_form_header($PROGRAM, "Submission form"),
    make_submission_form_top($action, $logo, $alt, $form_description),
    $required,
    $optional, 
    make_submit_button("Start search", $email_contact),
    make_submission_form_bottom(),
    make_submission_form_tailer()
  );
  print "Content-Type: text/html \n\n$form";


} # print_form

#
# Make verification message in HTML
#
sub make_verification {

  my $content = "<HR>\n<UL>\n";

  $content .= "<LI> Description: <B>$description</B>\n" if $description;
  $content .= "<LI> Motif file name: <B>$motif_file_name</B>\n" if $motif_file_name;
  $content .= "<LI> Number of motifs: <B>$nmotifs</B>\n" if $nmotifs;
  $content .= "<LI> Total motif columns: <B>$total_cols</B>\n" if $total_cols;
  $content .= "<LI> Motif alphabet: <B>$motif_alphabet</B>\n" if $motif_alphabet;
  $content .= "<LI> Database to search: <B>$db_menu_name ($db)</B>\n" if $db_menu_name;
  $content .= "$db_descr\n" if $db_descr;
  $content .= "<LI> Motif match display threshold: <B><I>p</I>-value $output_pthresh</B>\n"; 
  $content .= "<LI> Scan given strand only</B>\n" if $norc; 
  $content .= "</UL>\n";

  return($content);

} # make_verification

#
# get parameters from the input
#
sub get_params {

  # command options
  $options = "--verbosity 0";

  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $description = param('description');
  $inline_motifs = param('inline_motifs');
  $inline_name = param('inline_name');
  $output_pthresh = param('output_pthresh'); 
  $options .= " --output-pthresh $output_pthresh";
  $norc = param('norc');  
  $options .= ' --norc ' if $norc;
  @database = split(/,/, param('database'));
  $upload_db_name = param('upload_db');
  #$motif_file_name = param('motifs'); 
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs');
  
  # check list of DBs
  #open(DB_NAMES, "<$dbnames") || whine("Cannot open list of databases.");
  #close(DB_NAMES);

} # get_params

#
# check the user input
#
sub check_params {

  # set working directory to LOGS directory
  chdir($logs) || &whine("Can't cd to $logs");

  # check that valid email address was provided
  check_address($address);

  # check description field
  check_description($description);

  #check pvalue threshold
  check_pvalue($output_pthresh);

  # create file containing motifs
  $motif_data = upload_motif_file($motif_file_name, $inline_motifs);

  # check motifs and get the alphabet
  ($motifs_found, $motif_alphabet, $nmotifs, $total_cols) = 
    check_meme_motifs("meme-io", $motif_data) if ($motif_data);

  # check norc option
  check_dna_only($motif_alphabet, "Scan given strand only", "with DNA motifs") if ($norc && defined $motif_alphabet);

  # choose the appropriate database
  if ($motif_alphabet) {
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr) =
      get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, @database);
  }

} # check_params

#
# Submit job to OPAL
#
sub submit_opal
{
  $fimo = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  my $args = "motifs $db $options";

  #
  # start OPAL request
  #
  $req = JobInputType->new();
  $req->setArgs($args);
  #
  # create list of OPAL file objects
  #
  @infilelist = ();
  # 1) Motif File
  $inputfile = InputFileType->new("motifs", $motif_data);
  push(@infilelist, $inputfile);
  # 2) Uploaded sequence file
  if ($db_alphabet ne "local") {
    $inputfile = InputFileType->new($db, $uploaded_data);
    push(@infilelist, $inputfile);
  }
  # Email address file (for logging purposes only)
  $inputfile = InputFileType->new("address_file", $address);
  push(@infilelist, $inputfile);
  # Submit time file (for logging purposes only)
  $inputfile = InputFileType->new("submit_time_file", `date -u '+%d/%m/%y %H:%M:%S'`);
  push(@infilelist, $inputfile);

  # Add file objects to request
  $req->setInputFile(@infilelist);

  # Submit the request to OPAL
  $result = $fimo->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  verify_opal_job($result, $address, $email_contact, $verify);

} # submit
