#!@WHICHPERL@
# Author: Timothy Bailey
##
## $Id: mast.pl 4514 2010-05-04 02:28:44Z james_johnson $
##
## $Log:$
## $Rev: 4514 $ $Date: 2010-05-04 12:28:44 +1000 (Tue, 04 May 2010) $ $Author: james_johnson $

use lib qw(@PERLLIBDIR@);
use Globals;
use Validation;
use CGI qw/:standard/;          # use the CGI package
use SOAP::Lite;
use MIME::Base64;
use OpalServices;
use OpalTypes;
use CGI::Carp qw( fatalsToBrowser );
require "Utils.pm";

# get the directories using the new installation scheme
$dir = "@MEME_DIR@";	# installed directory
$logs = "$dir/LOGS";	# directory for logs
$bin = "$dir/bin";	# binary directory
$service_url = "@OPAL@/MAST_@S_VERSION@";	
$email_contact = '@contact@';
$short_dna_only = 1;            # MAST cannot scan long DNA sequences
$names_file = 'fasta_db.csv';
$dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file&short_only=$short_dna_only";
$db_doc_url = $dbnames_url . "&doc=1";

# defaults
$PROGRAM = "MAST";
$NERRORS = 0;				 # no errors yet
$MAX_UPLOAD_SIZE = 1000000;

# get the parameters for the query
get_params();

# if there is no action specified, print an input form
if (! $action) {
  print_form() unless $NERRORS;
} else {
  check_params();
  submit() unless $NERRORS;
  print_tailers();
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
  my $inline_descr = "MEME motifs from sequences in file <TT>'$inline_name'</TT>.";
  my $req_left = make_address_field($address, $address_verify);
  $req_left .= "<BR>\n";
  $req_left .= make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet, 
    $motif_doc, $motif_example, $inline_descr, $example_type);

  # required right side: database

  my $sample_db = "../examples/sample-kabat.seq";
  my $req_right = make_supported_databases_field("database", $dbnames_url, $db_doc_url);
  $req_right .= "<BR>\n<B>or</B>\n<BR>\n";
  $req_right .= make_upload_fasta_field("upload_db", $MAX_UPLOAD_SIZE, $sample_db);

  my $required = make_input_table("Required", $req_left, $req_right);

  # 
  # optional arguments
  #

  # optional left: description, ev, mev, use_seq_comp
  my $descr = "motifs";
  my $opt_left = make_description_field($descr, $description);

  # make EV field
  my $doc = "../mast-input.html#ev";
  my $ev_text = "<A HREF='$doc'><B>Display sequences</A></B> with <I>E</I>-value below:";
  my $ev_name = "ev";
  my $ev_selected = "10";
  my @ev_list = ("0.01", "0.1", "1", "10", "20", "50", "100", "200", "500", "1000");
  $opt_left .= "<BR>\n<P>\n" .
    make_select_field($ev_text, $ev_name, $ev_selected, @ev_list) .
    "</P>\n";

  # make MEV field
  my $doc = "../mast-input.html#mev";
  my $mev_text = "<A HREF='$doc'><B>Ignore motifs</A></B> if <I>E</I>-value above:";
  my $mev_name = "mev";
  my $mev_selected = "use all motifs";
  my @mev_list = ("use all motifs", "100", "50", "20", "10", "5", "2", "1", "0.5", "0.2", "0.1", "0.05", "0.02", "0.01", "0.005", "0.002", "0.001", "1e-5", "1e-10", "1e-50", "1e-100");
  $opt_left .= "<P>\n" .
    make_select_field($mev_text, $mev_name, $mev_selected, @mev_list) .
    "</P>\n";

  # make seq_comp checkbox
  my $doc = "../mast-input.html#use_seq_comp";
  my $descr = "Use individual <A HREF='$doc'><B> sequence composition</B></A><BR>in <I>E</I>- and <I>p</I>-value calculation";
  $opt_left .= make_checkbox("use_seq_comp", "-comp", $descr, 0);

  # optional right side: scale, translate, dna-only

  # make scale threshold checkbox
  my $doc = "../mast-input.html#use_seq_p";
  my $descr = "<A HREF='$doc'><B>Scale motif display threshold</B></A> by sequence length";
  my $opt_right = make_checkbox("use_seq_p", "-seqp -mt 0.01", $descr, 0);

  # make translate DNA checkbox
  my $doc = "../mast-input.html#dna";
  my $descr = "<A HREF='$doc'><B>Search nucleotide</B></A> database with protein motifs";
  $opt_right .= "<P>\n" .
    make_checkbox("dna", "-dna", $descr, 0) .
    "</P>\n";

  # DNA options
  if (!defined $alphabet || $alphabet eq "ACGT") {
    # make strands menu
    my $doc = "../mast-input.html#strands";
    my $text = "Treatment of <A HREF='$doc'><B>reverse complement</B></A> strands:<BR>";
    my $name = "strands";
    my $selected = ",combine with given strand";
    my @list = ($selected, "-sep,treat as separate sequence", "-norc,none");
    my $options = make_select_field($text, $name, $selected, @list);
    $opt_right .= make_dna_only($options);
  } # dna options

  my $optional = make_input_table("Optional", $opt_left, $opt_right);

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
  $pvalue_text = $use_seq_p ? "sequence" : "motif";
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
  @database = split(/,/, param('database'));
  $upload_db_name = param('upload_db');
  $translate_dna = param('dna'); $dna_options .= " $translate_dna" if $translate_dna;
  $use_seq_p = param('use_seq_p'); $options .= " $use_seq_p" if $use_seq_p;
  $ev = param('ev'); $options .= " -ev $ev" if $ev;
  $use_seq_comp = param('use_seq_comp'); $options .= " $use_seq_comp" if $use_seq_comp;
  ($strands, $strands_text) = split(/,/, param('strands')); $dna_options .= " $strands" if $strands;
  $mev = param('mev') + 0; $options .= " -mev $mev" if $mev;

  # name and file handle of motif file
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs');

  $options .= " -remcorr";
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

  # create file containing motifs
  $motif_data = upload_motif_file($motif_file_name, $inline_motifs);

  # check motifs and get the alphabet
  ($motifs_found, $motif_alphabet, $nmotifs, $total_cols) = 
    check_meme_motifs("meme-io", $motif_data) if ($motif_data);

  # choose the appropriate database
  if ($motif_alphabet) {
    my $short_dna_only = 1;	# mast can't scan long DNA
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $db_menu_name, $db_long_name, $db_descr) =
      get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, @database);
  }

  # check validity of options
  if ($translate_dna && $motif_alphabet ne "PROTEIN") {
    whine("You specified searching a nucleotide database with protein motifs, but your motifs are DNA.  Go back and uncheck that box.\n");
  }

} # check_params

#
# Submit job to OPAL
#
sub submit
{
  $mast = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  $options .= $dna_options if ($target_alphabet eq "DNA");

  #my $args = "motifs $db -df '$db_menu_name' $options";
  my $args = "motifs $db -df '$db_root' $options";

  #
  # start OPAL request
  #
  $req = JobInputType->new();
  $req->setArgs($args);

  #
  # create list of OPAL file objects
  #
  @infilelist = ();
  # Motif File
  $inputfile = InputFileType->new("motifs", $motif_data);
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

  # Submit the request to OPAL
  $result = $mast->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  verify_opal_job($result, $address, $email_contact, $verify);

} # submit
