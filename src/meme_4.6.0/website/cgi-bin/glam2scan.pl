#!@WHICHPERL@
##
## $Id:$
## $Log:$
## $Rev:$
##
# Author: Timothy Bailey

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);

use CGI qw(:standard);
use MIME::Base64;
use SOAP::Lite;

use Globals;
use OpalServices;
use OpalTypes;
use MemeWebUtils;

# get the directories using the new installation scheme
my $logs = "@MEME_DIR@/LOGS";		# log file directory
my $bin = "@MEME_DIR@/bin";		# binary directory
my $service_url = "@OPAL@/GLAM2SCAN_@S_VERSION@";	
#$email_contact = "&#x67;&#108;&#x61;&#x6d;&#50;&#x40;&#105;&#x6d;&#x62;&#x2e;&#117;&#x71;&#46;&#101;&#100;&#x75;&#46;&#97;&#117;";
my $email_contact = '@contact@';
my $index_path = "@MEME_DIR@/etc/fasta_db.index";
my $query_url = "get_db_list.cgi?db_names=fasta_db.csv&mode=xml&catid=~category~";
my $doc_url =   "get_db_list.cgi?db_names=fasta_db.csv&mode=doc";

# constant globals
my $PROGRAM = "GLAM2SCAN";
my $MAX_UPLOAD_SIZE = 1000000;
my $MAX_ALIGNMENTS = 200;

# globals from get_params
my ($options, $dna_options, $action, $address, $address_verify, $description, $naligns, 
  $alphabet, $inline_aln, $inline_name, $upload_db_name, $both_str, $aln_file_name);
my $database_id;
# globals from check_params
my ($aln_data, $motif_alphabet, $db, $uploaded_data, $db_alphabet, $target_alphabet, 
  $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr, $args);
# utilities
my $utils = new MemeWebUtils($PROGRAM, $bin);

# get the parameters for the query
get_params();

# if there is no action specified, print an input form
if (! $action) {
  print_form() unless $utils->has_errors();
} else {
  check_params();
  submit_to_opal() unless $utils->has_errors();
  $utils->print_tailers();
}

##############################
# SUBROUTINES 		     #
##############################

#
# print the input form
#
# Note: if query was included when this was called, that part of form
# is omitted.
#
sub print_form
{
  my $action = "glam2scan.cgi";
  my $logo = "../doc/images/glam2scan_logo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = "Use this form to submit a GLAM2 motif to <B>$PROGRAM</B> to be used in searching a sequence database.";

  #
  # required section
  #
  # required left side: address, motif, alphabet
  my $aln_doc = "../doc/glam2_ref.html";
  my $aln_example = "../examples/sample-glam2-aln.dna";
  my $example_type = "GLAM2 nucleotide motif";
  my $inline_descr = "GLAM2 motif <TT>'" . (defined $inline_name ? $inline_name : "") . "'</TT>.";
  my $req_left = $utils->make_address_field($address, $address_verify);
  $req_left .= "<BR>\nSpecify the name of your motif file--a GLAM2 output file (HTML or text format).\n";
  $req_left .= $utils->make_motif_field($PROGRAM, "aln", $inline_aln, $alphabet,
    $aln_doc, $aln_example, $inline_descr, $example_type);
  unless ($alphabet) {
    my $text = "<BR>\n<BR>\nSpecify the <B>alphabet</B> of your motif:";
    my @values = (",", "n,nucleotide", "p,protein");
    $req_left .= $utils->make_select_field($text, "alphabet", "", @values);
  }

  # required right side: database
  my $fasta_doc = "../doc/fasta-format.html";
  my $sample_db = "../examples/sample-kabat.seq";
  my $req_right .= $utils->make_supported_databases_field("database", "Sequence", $index_path, $query_url, $doc_url);
  $req_right .= "<BR>\n<B>or</B>\n<BR>\n";
  $req_right .= $utils->make_upload_fasta_field("upload_db", $MAX_UPLOAD_SIZE, $fasta_doc, $sample_db);

  my $required = $utils->make_input_table("Required", $req_left, $req_right);

  # 
  # optional arguments
  #

  # optional left: description, -n
  my $descr = "motif";
  my $opt_left = $utils->make_description_field($descr, $description);

  $opt_left .= qq {
<BR>
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="-n" VALUE="25">
<B>Number of alignments</B> to report (<= $MAX_ALIGNMENTS)
  }; # end_quote


  # optional right side: dna-only

  # print DNA-ONLY options if motif is nucleotide or undefined
  my $opt_right = "";
  if (!defined $alphabet || $alphabet eq "n") {
    my $text = "<B>Examine both strands</B> - forward and reverse complement";
    my $options .= $utils->make_checkbox("-2", "1", $text, 1);
    $opt_right .= $utils->make_dna_only($options)
  } # print DNA-ONLY options


  my $optional = $utils->make_input_table("Optional", $opt_left, $opt_right);

  #
  # print final form
  #
  my $form = $utils->make_submission_form(
    $utils->make_form_header($PROGRAM, "Submission form"),
    $utils->make_submission_form_top($action, $logo, $alt, $form_description),
    $required,
    $optional,
    $utils->make_submit_button("Start search", $email_contact),
    $utils->make_submission_form_bottom(),
    $utils->make_submission_form_tailer()
  );
  print "Content-Type: text/html\n\n$form";

} # print_form

#
# Create an HTML verification message
#
sub make_verification
{
  my $content = "<HR><UL>\n";

  $content .= "<LI> Description: <B>$description</B>\n" if $description;
  $content .= "
    <LI> Alignment file name: <B>$aln_file_name</B>
    <LI> Alphabet: <B>$motif_alphabet</B>
    <LI> Database to search: <B>$db_menu_name ($db)</B>
    $db_descr
  ";
  $content .= "<LI> Searching: " . ($both_str ? "<B>forward and reverse complement strands</B>\n" : "<B>forward strand only<B>\n") 
      if ($motif_alphabet eq "DNA");
  $content .= "</UL>\n";

  return($content);

} # make_verification

#
# get parameters from the input
#
sub get_params
{

  # command options
  $options = "";
  $dna_options = "";

  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $description = param('description');
  $naligns = param('-n'); $options .= " -n $naligns" if $naligns; 
  ($alphabet) = defined param('alphabet') ? split(/,/, param('alphabet')) : ();
  $inline_aln = param('inline_aln');
  $inline_name = param('inline_name'); # $options .= " -mf '$inline_name'" if $inline_name;
  $database_id = param('database');
  $upload_db_name = param('upload_db');

  $both_str = param('-2'); $dna_options .= " -2" if $both_str;

  # uploaded file
  if ($inline_aln) {
    $aln_file_name = $inline_name;
  }
  elsif (param('motifs')) {
    #Motif search from names this field 'motifs'.
    $aln_file_name = param('motifs');
  }
  else {
    $aln_file_name = param('aln');
  }

  # get the list of databases; done here so that whine works
  #open(DB_NAMES, "<$dbnames") || $utils->whine("Cannot open list of databases.");

} # get_params

#
# check the user input
#
sub check_params
{
  # set working directory to LOGS directory
  chdir($logs) || $utils->whine("Can't cd to $logs");

  # check that valid email address was provided
  $utils->check_address($address);

  # check description field
  $utils->check_description($description);

  # check number of alignments
  if ($naligns < 1 || $naligns > $MAX_ALIGNMENTS) {
    $utils->whine("You must specify a number of alignments between 1 and $MAX_ALIGNMENTS, inclusive.");
  }

  # create file containing motif 
  $aln_data = $utils->upload_motif_file($aln_file_name, $inline_aln);

  # check that alphabet was provided
  if (!$alphabet) {
    $utils->whine("You must specify the alphabet of your motif.  Please enter it on the form and retry.");
  }

  # Choose the appropriate database
  if ($alphabet) {
    $motif_alphabet = $alphabet eq "p" ? "PROTEIN" : "DNA"; 
    my $translate_dna = 0;		# GLAM2SCAN can't scan DNA with protein
    my $short_dna_only = 0;		# GLAM2SCAN can scan long DNA
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr) =
      $utils->get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, $index_path, $database_id);
  }

  # finish GLAM2SCAN argument list
  $options .= $dna_options if ($alphabet eq "n"); 
  $args = "$options $alphabet aln $db";

} # check_params

#
# Submit the job to OPAL
#
sub submit_to_opal
{
  my $service = OpalServices->new(service_url => $service_url);

  #
  # start OPAL request
  #
  my $req = JobInputType->new();
  $req->setArgs($args);

  #
  # create list of OPAL file objects
  #
  my @infilelist = ();
  # Alignment File
  my $inputfile = InputFileType->new("aln", $aln_data);
  push(@infilelist, $inputfile);
  # Uploaded sequence file
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

  # output HTTP header
  print "Content-type: text/html\n\n";

  # Submit the request to OPAL
  my $result = $service->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  $utils->verify_opal_job($result, $address, $email_contact, $verify);

} # submit_to_opal

