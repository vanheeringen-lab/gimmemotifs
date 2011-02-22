#!@WHICHPERL@

# Author: Charles Grant

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);

use CGI qw(:standard);
use XML::Simple;

use Globals;
use OpalServices;
use OpalTypes;
use MemeWebUtils;

# get the directories using the new installation scheme
my $logs = "@MEME_DIR@/LOGS";	# directory for logs
my $bin = "@MEME_DIR@/bin";	# binary directory
my $service_url = "@OPAL@/MCAST_@S_VERSION@";
my $email_contact = '@contact@';
my $index_path = "@MEME_DIR@/etc/fasta_db.index";
my $query_url = "get_db_list.cgi?db_names=fasta_db.csv&mode=xml&catid=~category~";
my $doc_url =   "get_db_list.cgi?db_names=fasta_db.csv&mode=doc";

# constant globals
my $PROGRAM = "MCAST";
my $MAX_UPLOAD_SIZE = 1000000;

# globals for get_params
my ($action, $address, $address_verify, $description, $inline_motifs, 
  $inline_name, $motif_hit_pthresh, $options, $max_gap, $pseudocount_weight, 
  $evalue_output_thresh, $upload_db_name, $motif_file_name);
my $database_id;

# globals for check_params
my ($motif_data, $motifs_found, $motif_alphabet, $nmotifs, $total_cols, $db, 
  $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, 
  $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr);

# utils
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
sub print_form {
  my $action = "mcast.cgi";
  my $logo = "../doc/images/mcast-logo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = "Use this form to submit motifs to <B>$PROGRAM</B> to be used in searching a sequence database.";

  #
  # required section
  #

  # required left side: address, motifs
  my $motif_doc = "../doc/meme-format.html";
  my $motif_example = "../examples/sample-dna-motif.meme-io";
  my $example_type = "DNA minimal motif format file";
  my $inline_descr = "MEME motifs from sequences in file <TT>'".(defined $inline_name ? $inline_name : "")."'</TT>.";
  my $req_left = $utils->make_address_field($address, $address_verify);
  $req_left .= "<BR>\n";
  # the alphabet could be extracted from inline_motifs, but this is closest to 
  # what was happening when it was not declared at all.
  my $alphabet = undef; 
  $req_left .= $utils->make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet,
    $motif_doc, $motif_example, $inline_descr, $example_type);

  # required right side: database
  my $fasta_doc = "../doc/fasta-format.html";
  my $sample_db = "../examples/sample-kabat.seq";
  my $req_right = $utils->make_supported_databases_field("database", "Sequence", $index_path, $query_url, $doc_url);
  $req_right .= "<BR>\n<B>or</B>\n<BR>\n";
  $req_right .= $utils->make_upload_fasta_field("upload_db", $MAX_UPLOAD_SIZE, $fasta_doc, $sample_db);

  my $required = $utils->make_input_table("Required", $req_left, $req_right);

  # 
  # optional arguments
  #

  # optional left: description
  my $descr = "motifs";
  my $opt_left = $utils->make_description_field($descr, $description);


  # optional right side:
  # Text box for motif hit p-value threshold
  my $opt_right = make_text_field(
                  '<a href="../doc/mcast.html#p-thresh"><i>P</i>-value threshold for motif hits</a>',
                  "motif_hit_pthresh",
                  "0.0005",
                  "5"
                );
  # Text box for max allowed distance between adj. hits.
  $opt_right .= "<p />\n";
  $opt_right .= make_text_field(
                  '<a href="../doc/mcast.html#max-gap">Maximum allowed distance between adjacent hits</a>',
                  "max_gap",
                  "50",
                  "5"
                );
  # E-value output threshold
  my @ev_list = ("0.01", "0.1", "1", "10", "100", "1000", "10000");
  my $ev_selected = "10";
  $opt_right .= "<p />\n";
  $opt_right .= $utils->make_select_field(
                  '<a href="../doc/mcast.html#e-thresh">Print matches with <i>E</i>-values less than this value</a> ', 
                  'evalue_output_thresh', 
                  $ev_selected, 
                  @ev_list
  );
  # Text box for pseudocount weight.
  $opt_right .= "<p />\n";
  $opt_right .= make_text_field(
                  '<a href="../doc/mcast.html#bg-weight">Pseudocount weight</a>',
                  "pseudocount_weight",
                  "4",
                  "5"
                );

  my $optional = $utils->make_input_table("Optional", $opt_left, '<td style="width: 80%;">' .
      $opt_right . '</td>');

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
  print "Content-Type: text/html \n\n$form";


} # print_form

#
# make a text field
#
sub make_text_field {
  my ($title, $name, $value, $width) = @_;

  my $content = qq {
$title
<input class="maininput" type="text" size="$width" name="$name" VALUE="$value">
<BR>
  }; # end quote
  return($content);
} # make_text_field

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
  $content .= "<LI> Motif hit p-value threshold: $motif_hit_pthresh</B>\n"; 
  $content .= "<LI> Maximum gap between motifs: $max_gap</B>\n"; 
  $content .= "<LI> E-value output threshold: $evalue_output_thresh</B>\n"; 
  $content .= "<LI> Psuedocount weight: $pseudocount_weight</B>\n"; 
  $content .= "</UL>\n";

  return($content);

} # make_verification

#
# get parameters from the input
#
sub get_params {


  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $description = param('description');
  $inline_motifs = param('inline_motifs');
  $inline_name = param('inline_name');
  # command options
  $motif_hit_pthresh = param('motif_hit_pthresh'); 
  $options = "";
  $options .= " -p-thresh $motif_hit_pthresh" if defined $motif_hit_pthresh;
  $max_gap = param('max_gap'); 
  $options .= " -max-gap $max_gap" if defined $max_gap;
  $pseudocount_weight = param('pseudocount_weight'); 
  $options .= " -bgweight $pseudocount_weight" if defined $pseudocount_weight;
  $evalue_output_thresh = param('evalue_output_thresh'); 
  $options .= " -e-thresh $evalue_output_thresh" if defined $evalue_output_thresh;
  $database_id = param('database');
  $upload_db_name = param('upload_db');
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs');
  
} # get_params

#
# check the user input
#
sub check_params {

  # set working directory to LOGS directory
  chdir($logs) || $utils->whine("Can't cd to $logs");

  # check that valid email address was provided
  $utils->check_address($address);

  # check description field
  $utils->check_description($description);

  # create file containing motifs
  $motif_data = $utils->upload_motif_file($motif_file_name, $inline_motifs);

  # check motifs and get the alphabet
  ($motifs_found, $motif_alphabet, $nmotifs, $total_cols) = 
    $utils->check_meme_motifs("meme-io", $motif_data) if ($motif_data);

  # I'm guessing mcast doesn't have the option of translating dna into protein
  my $translate_dna = 0;

  # choose the appropriate database
  if ($motif_alphabet) {
    my $short_dna_only = 0;		# MCAST can scan long DNA sequences
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr) =
      $utils->get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, $index_path, $database_id);
  }

} # check_params

#
# Submit job to OPAL
#
sub submit_to_opal
{
  my $mcast = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  my $args = "motifs $db $options";

  #
  # start OPAL request
  #
  my $req = JobInputType->new();
  $req->setArgs($args);
  #
  # create list of OPAL file objects
  #
  my @infilelist = ();
  # 1) Motif File
  my $inputfile = InputFileType->new("motifs", $motif_data);
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

  # output HTTP header
  print "Content-type: text/html\n\n";

  # Submit the request to OPAL
  my $result = $mcast->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  $utils->verify_opal_job($result, $address, $email_contact, $verify);

} # submit_to_opal
