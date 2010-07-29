#!@WHICHPERL@

# Author: Timothy Bailey

use lib qw(@PERLLIBDIR@);
use Globals;
use Validation;
use CGI qw/:standard/;          # use the CGI package
use SOAP::Lite;
use MIME::Base64;
use OpalServices;
use OpalTypes;
require "Utils.pm";

# get the directories using the new installation scheme
$dir = "@MEME_DIR@";	# installed directory
$logs = "$dir/LOGS";	# directory for logs
$bin = "$dir/bin";	# binary directory
$service_url = "@OPAL@/GOMO_@S_VERSION@";	
$email_contact = '@contact@';
$names_file = "gomo_db.csv";
$dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file&short_only=$short_dna_only";
$db_doc_url = $dbnames_url . "&doc=1";

# defaults
$PROGRAM = "GOMO";
$NERRORS = 0;				 # no errors yet
$THRESHOLD_DEFAULT = 0.05;
$THRESHOLD_MAX = 0.5;
$SHUFFLE_SCORES_DEFAULT = 100;
#$MAX_UPLOAD_SIZE = 1000000;

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
#############################

#
# print the input form
#
sub print_form {
  my $tab = " ";
  my $action = "gomo.cgi";
  my $logo = "../images/gomo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = make_description_message(11, $tab); 

  #
  # required section
  #

  # required left side: address, motifs
  my $motif_doc = "../doc/meme-format.html";
  my $motif_example = "../examples/sample-dna-motif.meme-io";
  my $example_type = "DNA minimal motif format file";
  my $inline_descr = "MEME motifs from sequences in file <tt>'$inline_name'</tt>.";
  my $req_left = make_address_field($address, $address_verify, 11, $tab);
  $req_left .= make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet, 
    $motif_doc, $motif_example, $inline_descr, $example_type, 11, $tab);

  # required right side: database
  my $req_right = make_supported_databases_field("database", $dbnames_url, $db_doc_url, 0, 12, $tab);
  
  my $multi_species = "<a href=\"../help_multiple_genomes.html\"><b>Use multiple genomes:</b></A>";
  my $checked = "yes";
  my @values = ("Whenever avaliable,yes", "Single genome only<br />,no");
  $req_right .= make_radio_field($multi_species, "multispecies", $checked, \@values, 12, $tab);
 
  $req_right .= make_significance_threshold(12, $tab);

  my $required = make_input_table("", $req_left, $req_right, 8, $tab);

  # 
  # optional arguments
  #

  # optional left: description
  my $descr = "motifs";
  my $opt_left = make_description_field($descr, $description);

  #
  # finally print form
  #
  my $form = make_submission_form(
    make_form_header($PROGRAM, "Submission form", "", 0, $tab),
    make_submission_form_top($action, $logo, $alt, $form_description, 2, $tab),
    $required,
    "", 
    make_submit_button("Start search", $email_contact, 8, $tab),
    make_submission_form_bottom(6, $tab),
    make_submission_form_tailer(0, $tab),
    1,    #indent_lvl
    $tab  #tab character
  );
  print "Content-Type: text/html \n\n$form";


} # print_form


#
# Make Significance threshold
#
sub make_significance_threshold {
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content =
      $indent."<br />\n".
      $indent."<br />\n".
      $indent."<b>Significance threshold</b>:<br />\n".
      $indent."<i>q</i>-value &le; <input class=\"maininput\" type=\"text\" ".
          "size=\"4\" name=\"threshold\" value=\"$THRESHOLD_DEFAULT\" /> ".
          "(maximum $THRESHOLD_MAX)<br />\n".
      $indent."<br />\n";
  return($content);
}

#
# Make description message
#
sub make_description_message {
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $content =
      $indent."Use this form to submit motifs to <b>$PROGRAM</b>.\n".
      $indent."$PROGRAM will use the DNA binding motifs to find putative target genes\n".
      $indent."and analyze their associated GO terms. A list of significant GO terms\n".
      $indent."that can be linked to the given motifs will be produced.\n".
      $indent."<a href=\"../gomo-intro.html\">Click here for more documentation.</a>\n";
  return($content);
}


#
# Make verification message in HTML
#
sub make_verification {

  my $content = "<hr />\n<ul>\n";

  $content .= "<li> Description: <b>$description</b></li>\n" if $description;
  $content .= "<li> Motif file name: <b>$motif_file_name</b></li>\n" if $motif_file_name;
  $content .= "<li> Number of motifs: <b>$nmotifs</b></li>\n" if $nmotifs;
  $content .= "<li> Total motif columns: <b>$total_cols</b></li>\n" if $total_cols;
  $content .= "<li> Motif alphabet: <b>$motif_alphabet</b></li>\n" if $motif_alphabet;
  $content .= "<li> Significance threshold:<b> $threshold </b></li>\n" if $threshold;
  #$content .= "<li> TF score function:<b> $score_function </b></li>\n" if $score_function;
  $content .= "<li> Database to search: " if ($db_menu_name || $db_descr);
  $content .= "<b>$db_menu_name ($db)</b>" if $db_menu_name;
  $content .= "\n$db_descr" if $db_descr;
  $content .= "</ul>\n";

  return($content);

} # make_verification

#
# get parameters from the input
#
sub get_params {
  $options = '';		# global
  
  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $description = param('description');
  $inline_motifs = param('inline_motifs');
  $inline_name = param('inline_name'); 
  #  $output_pthresh = param('output_pthresh');
  @database = split(/,/, param('database'));
  #  $upload_db_name = param('upload_db');
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs'); 
  $threshold = param('threshold'); $options .= " --t $threshold" if $threshold;
  $shuffle_scores = param('shuffle_scores'); 
  $shuffle_scores = $SHUFFLE_SCORES_DEFAULT unless $shuffle_scores;
  $options .= " --shuffle_scores $shuffle_scores"; 
  $is_multispecies = (param('multispecies') eq 'yes');
  #$score_function = param('score_function');

} # get_params

#
# check the user input
#
sub check_params {
  my $local_extra_dbs;
  # set working directory to LOGS directory
  chdir($logs) || &whine("Can't cd to $logs");

  # check that valid email address was provided
  check_address($address);

  # check description field
  check_description($description);
  
  # check that threshold is OK
  &check_threshold;
  
  # create file containing motifs
  $motif_data = upload_motif_file($motif_file_name, $inline_motifs);

  # check motifs and get the alphabet
  ($motifs_found, $motif_alphabet, $nmotifs, $total_cols) = 
    check_meme_motifs("meme-io", $motif_data) if ($motif_data);

  # choose the appropriate database
  if ($motif_alphabet) {
    my $short_dna_only = 0;	# gomo scan long DNA
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, 
        $short_seqs, $db_menu_name, $db_long_name, $db_descr, $local_extra_dbs) =
        get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, @database);
    if ($is_multispecies) {
      @extra_dbs = @$local_extra_dbs;
    } else {
      @extra_dbs = ();
    }
  }


} # check_params

#
# check that threshold is OK
#
sub check_threshold {
  whine("The threshold you specified ($threshold) is too small.<br />") if $threshold <= 0;
  whine("The threshold you specified ($threshold) is too large.<br />") if $threshold > $THRESHOLD_MAX;
} # check_threshold

#
# Submit job to OPAL
#
sub submit
{
  $gomo = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  #@extra_dbs = ("styph.na", "sflex.na");
  my $args = "motifs $db @extra_dbs END_OF_DBS $options";

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
  $result = $gomo->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  verify_opal_job($result, $address, $email_contact, $verify);

} # submit

