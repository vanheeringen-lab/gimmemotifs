#!@WHICHPERL@
# Author: Timothy Bailey
##
## $Id: mast.pl 5280 2011-01-05 11:10:34Z james_johnson $
##
## $Log:$
## $Rev: 5280 $ $Date: 2011-01-05 03:10:34 -0800 (Wed, 05 Jan 2011) $ $Author: james_johnson $

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
my $logs = "@MEME_DIR@/LOGS";	# directory for logs
my $bin = "@MEME_DIR@/bin";	# binary directory
my $service_url = "@OPAL@/MAST_@S_VERSION@";	
my $email_contact = '@contact@';
my $index_path = "@MEME_DIR@/etc/fasta_db.index";
my $query_url = "get_db_list.cgi?db_names=fasta_db.csv&mode=xml&catid=~category~&short_only=1";
my $doc_url =   "get_db_list.cgi?db_names=fasta_db.csv&mode=doc&short_only=1";

# constant globals
my $PROGRAM = "MAST";
my $MAX_UPLOAD_SIZE = 1000000;

# globals for get_params
my ($options, $dna_options, $action, $address, $address_verify, $description, 
  $alphabet, $inline_motifs, $inline_name, $upload_db_name, $translate_dna, $use_seq_p, $ev, 
  $use_seq_comp, $strands, $strands_text, $mev, $motif_file_name);
my $database_id;

# globals for check_params
my ($motif_data, $motifs_found, $motif_alphabet, $nmotifs, $total_cols, $db, 
  $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, 
  $short_seqs, $db_menu_name, $db_long_name, $db_descr);

#utils
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
  my $action = "mast.cgi";
  my $logo = "../images/mast.png";
  my $alt = "$PROGRAM logo";
  my $form_description = "Use this form to submit motifs to <B>$PROGRAM</B> to be used in searching a sequence database.";

  #
  # required section
  #

  # required left side: address, motifs
  my $motif_doc = "../mast-input.html#motif-file";
  my $motif_example = "../examples/sample-dna-motif.mast-io";
  my $example_type = "DNA motifs";
  my $inline_descr = "MEME motifs from sequences in file <TT>'".(defined $inline_name ? $inline_name : "")."'</TT>.";
  my $req_left = $utils->make_address_field($address, $address_verify);
  $req_left .= "<BR>\n";
  $req_left .= $utils->make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet, 
    $motif_doc, $motif_example, $inline_descr, $example_type);

  # required right side: database

  my $sample_db = "../examples/sample-kabat.seq";
  my $req_right = $utils->make_supported_databases_field("database", "Sequence", $index_path, $query_url, $doc_url);
  $req_right .= "<BR>\n<B>or</B>\n<BR>\n";
  $req_right .= $utils->make_upload_fasta_field("upload_db", $MAX_UPLOAD_SIZE, $sample_db);

  my $required = $utils->make_input_table("Required", $req_left, $req_right);

  # 
  # optional arguments
  #

  # optional left: description, ev, mev, use_seq_comp
  my $descr = "motifs";
  my $opt_left = $utils->make_description_field($descr, $description);

  # make EV field
  my $ev_doc = "../mast-input.html#ev";
  my $ev_text = "<A HREF='$ev_doc'><B>Display sequences</A></B> with <I>E</I>-value below:";
  my $ev_name = "ev";
  my $ev_selected = "10";
  my @ev_list = ("0.01", "0.1", "1", "10", "20", "50", "100", "200", "500", "1000");
  $opt_left .= "<BR>\n<P>\n" .
    $utils->make_select_field($ev_text, $ev_name, $ev_selected, @ev_list) .
    "</P>\n";

  # make MEV field
  my $mev_doc = "../mast-input.html#mev";
  my $mev_text = "<A HREF='$mev_doc'><B>Ignore motifs</A></B> if <I>E</I>-value above:";
  my $mev_name = "mev";
  my $mev_selected = "use all motifs";
  my @mev_list = ("use all motifs", "100", "50", "20", "10", "5", "2", "1", "0.5", "0.2", "0.1", "0.05", "0.02", "0.01", "0.005", "0.002", "0.001", "1e-5", "1e-10", "1e-50", "1e-100");
  $opt_left .= "<P>\n" .
    $utils->make_select_field($mev_text, $mev_name, $mev_selected, @mev_list) .
    "</P>\n";

  # make seq_comp checkbox
  my $seq_comp_doc = "../mast-input.html#use_seq_comp";
  my $seq_comp_descr = "Use individual <A HREF='$seq_comp_doc'><B> sequence composition</B></A><BR>in <I>E</I>- and <I>p</I>-value calculation";
  $opt_left .= $utils->make_checkbox("use_seq_comp", "-comp", $seq_comp_descr, 0);

  # optional right side: scale, translate, dna-only

  # make scale threshold checkbox
  my $use_seq_p_doc = "../mast-input.html#use_seq_p";
  my $use_seq_p_descr = "<A HREF='$use_seq_p_doc'><B>Scale motif display threshold</B></A> by sequence length";
  my $opt_right = $utils->make_checkbox("use_seq_p", "-seqp -mt 0.01", $use_seq_p_descr, 0);

  # make translate DNA checkbox
  my $translate_dna_doc = "../mast-input.html#dna";
  my $translate_dna_descr = "<A HREF='$translate_dna_doc'><B>Search nucleotide</B></A> database with protein motifs";
  $opt_right .= "<P>\n" .
    $utils->make_checkbox("dna", "-dna", $translate_dna_descr, 0) .
    "</P>\n";

  # DNA options
  if (!defined $alphabet || $alphabet eq "ACGT") {
    # make strands menu
    my $strands_doc = "../mast-input.html#strands";
    my $strands_text = "Treatment of <A HREF='$strands_doc'><B>reverse complement</B></A> strands:<BR>";
    my $strands_name = "strands";
    my $strands_selected = ",combine with given strand";
    my @strands_list = ($strands_selected, "-sep,treat as separate sequence", "-norc,none");
    my $strands_options = $utils->make_select_field($strands_text, $strands_name, $strands_selected, @strands_list);
    $opt_right .= $utils->make_dna_only($strands_options);
  } # dna options

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
  $content .= "<LI> Database to search (uploaded): <B>$upload_db_name </B>\n" if $upload_db_name;
  $content .= "$db_descr\n" if $db_descr;
  $content .= "<LI> Treatment of reverse complement strands: <B>$strands_text</B>\n" if $strands_text && $target_alphabet eq "DNA";
  $content .= "<LI> Searching <B>nucleotide</B> database with <B>protein</B> motifs\n" if $translate_dna;
  $content .= "<LI> Displaying sequences with <B><I>E</I>-value < $ev </B>\n" if $ev;
  $content .= "<LI> Adjusting p-values and E-values for <B>sequence composition</B>\n" if $use_seq_comp;
  my $pvalue_text = $use_seq_p ? "sequence" : "motif";
  $content .= "<LI> Motif display threshold: <B>$pvalue_text <I>p</I>-value < 0.0001 </B>\n"; 
  $content .= "<LI> Using only motifs with <B><I>E</I>-value < $mev</B>\n" if $mev;

  $content .= "</UL>\n";

  return($content);

} # make_verification

#
# get parameters from the input
#
sub get_params {

  # command options
  $options = "";
  $dna_options = "";

  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $description = param('description');
  $alphabet = param('alphabet');
  $inline_motifs = param('inline_motifs');
  $inline_name = param('inline_name'); $options .= " -mf '$inline_name'" if $inline_name;
  $database_id = param('database');
  $upload_db_name = param('upload_db');
  $translate_dna = param('dna'); $dna_options .= " $translate_dna" if $translate_dna;
  $use_seq_p = param('use_seq_p'); $options .= " $use_seq_p" if $use_seq_p;
  $ev = param('ev'); $options .= " -ev $ev" if $ev;
  $use_seq_comp = param('use_seq_comp'); $options .= " $use_seq_comp" if $use_seq_comp;
  ($strands, $strands_text) = defined param('strands') ? split(/,/, param('strands')) : (); $dna_options .= " $strands" if $strands;
  $mev = defined param('mev') ? param('mev') + 0 : 0; $options .= " -mev $mev" if $mev;

  # name and file handle of motif file
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs');

  $options .= " -remcorr";
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

  # choose the appropriate database
  if ($motif_alphabet) {
    my $short_dna_only = 1;	# mast can't scan long DNA
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr) =
      $utils->get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, $index_path, $database_id);
  }

  # check validity of options
  if ($translate_dna && $motif_alphabet ne "PROTEIN") {
    $utils->whine("You specified searching a nucleotide database with protein motifs, but your motifs are DNA.  Go back and uncheck that box.\n");
  }

} # check_params

#
# Submit job to OPAL
#
sub submit_to_opal
{
  my $mast = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  $options .= $dna_options if ($target_alphabet eq "DNA");

  #my $args = "motifs $db -df '$db_menu_name' $options";
  my $args = "motifs $db -df '$db_root' $options";

  #
  # start OPAL request
  #
  my $req = JobInputType->new();
  $req->setArgs($args);

  #
  # create list of OPAL file objects
  #
  my @infilelist = ();
  # Motif File
  my $inputfile = InputFileType->new("motifs", $motif_data);
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
  my $result = $mast->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  $utils->verify_opal_job($result, $address, $email_contact, $verify);

} # submit_to_opal
