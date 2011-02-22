#!@WHICHPERL@

# Author: Timothy Bailey

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);

use CGI qw/:standard/;          # use the CGI package
use MIME::Base64;
use SOAP::Lite;

use Globals;
use OpalServices;
use OpalTypes;
use MemeWebUtils;

# get the directories using the new installation scheme
my $logs = "@MEME_DIR@/LOGS";	# directory for logs
my $bin = "@MEME_DIR@/bin";	# binary directory
my $service_url = "@OPAL@/GOMO_@S_VERSION@";	
my $email_contact = '@contact@';
my $index_path = "@MEME_DIR@/etc/gomo_db.index";
my $query_url = "get_db_list.cgi?db_names=gomo_db.csv&mode=xml&catid=~category~";
my $doc_url =   "get_db_list.cgi?db_names=gomo_db.csv&mode=doc";

# constant globals
my $PROGRAM = "GOMO";
my $THRESHOLD_DEFAULT = 0.05;
my $THRESHOLD_MAX = 0.5;
my $SHUFFLE_SCORES_DEFAULT = 100;

# globals from get_params
my ($options, $action, $address, $address_verify, $description, $inline_motifs, 
  $inline_name, $motif_file_name, $threshold, $shuffle_scores, $is_multispecies);
my $database_id;

# globals from check_params
my ($motif_data, $motifs_found, $motif_alphabet, $nmotifs, $total_cols, $db, 
  $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, 
  $short_seqs, $db_menu_name, $db_long_name, $db_descr, $local_extra_dbs);
my @extra_dbs;

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
  my $inline_descr = "MEME motifs from sequences in file <tt>'".(defined $inline_name ? $inline_name : "")."'</tt>.";
  my $req_left = $utils->make_address_field($address, $address_verify, 11, $tab);
  #we could determine the alphabet from the inline_motifs (check_params does this) 
  #but this worked fine before when this variable wasn't even decared so for now I won't change it...
  my $alphabet = undef; 
  $req_left .= $utils->make_motif_field($PROGRAM, "motifs", $inline_motifs, $alphabet, 
    $motif_doc, $motif_example, $inline_descr, $example_type, 11, $tab);

  # required right side: database
  my $req_right = $utils->make_supported_databases_field("database", "Sequence", $index_path, $query_url, $doc_url);
  
  my $multi_species = "<a href=\"../help_multiple_genomes.html\"><b>Use multiple genomes:</b></A>";
  my $checked = "yes";
  my @values = ("Whenever avaliable,yes", "Single genome only<br />,no");
  $req_right .= $utils->make_radio_field($multi_species, "multispecies", $checked, \@values, 12, $tab);
 
  $req_right .= make_significance_threshold(12, $tab);

  my $required = $utils->make_input_table("", $req_left, $req_right, 8, $tab);

  # 
  # optional arguments
  #

  # optional left: description
  my $descr = "motifs";
  my $opt_left = $utils->make_description_field($descr, $description);

  #
  # finally print form
  #
  my $form = $utils->make_submission_form(
    $utils->make_form_header($PROGRAM, "Submission form", "", 0, $tab),
    $utils->make_submission_form_top($action, $logo, $alt, $form_description, 2, $tab),
    $required,
    "", 
    $utils->make_submit_button("Start search", $email_contact, 8, $tab),
    $utils->make_submission_form_bottom(6, $tab),
    $utils->make_submission_form_tailer(0, $tab),
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
  $database_id = param('database');
  #  $upload_db_name = param('upload_db');
  $motif_file_name = $inline_motifs ? $inline_name : param('motifs'); 
  $threshold = param('threshold'); $options .= " --t $threshold" if $threshold;
  $shuffle_scores = param('shuffle_scores'); 
  $shuffle_scores = $SHUFFLE_SCORES_DEFAULT unless $shuffle_scores;
  $options .= " --shuffle_scores $shuffle_scores"; 
  $is_multispecies = (defined param('multispecies') ? param('multispecies') eq 'yes' : 0);
  #$score_function = param('score_function');

} # get_params

#
# check the user input
#
sub check_params {
  my $local_extra_dbs;
  # set working directory to LOGS directory
  chdir($logs) || $utils->whine("Can't cd to $logs");

  # check that valid email address was provided
  $utils->check_address($address);

  # check description field
  $utils->check_description($description);
  
  # check that threshold is OK
  &check_threshold;
  
  # create file containing motifs
  $motif_data = $utils->upload_motif_file($motif_file_name, $inline_motifs);

  # check motifs and get the alphabet
  ($motifs_found, $motif_alphabet, $nmotifs, $total_cols) = 
    $utils->check_meme_motifs("meme-io", $motif_data) if ($motif_data);

  # uploads are not a possibility on gomo due to the complexity of setting up a promoter to gene map
  # this used to be a get_params global but it was disabled and seemed to be working without a declaration
  # and so I'll just set it to undef for now.
  my $upload_db_name = undef;
  #another variable that was never defined... sigh
  my $translate_dna = 0;

  # choose the appropriate database
  if ($motif_alphabet) {
    # gomo doesn't care if the sequence is marked as long 
    # (though the typical use case would be promoter sequences which are short)
    my $short_dna_only = 0; 
    ($db, $uploaded_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, 
        $short_seqs, $db_menu_name, $db_long_name, $db_descr, $local_extra_dbs) =
        $utils->get_db_name($motif_alphabet, $upload_db_name, $MAX_UPLOAD_SIZE, $translate_dna, $short_dna_only, $index_path, $database_id);
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
  $utils->whine("The threshold you specified ($threshold) is too small.<br />") if $threshold <= 0;
  $utils->whine("The threshold you specified ($threshold) is too large.<br />") if $threshold > $THRESHOLD_MAX;
} # check_threshold

#
# Submit job to OPAL
#
sub submit_to_opal
{
  my $gomo = OpalServices->new(service_url => $service_url);

  #
  # create the command
  #
  #@extra_dbs = ("styph.na", "sflex.na");
  my $args = "motifs $db @extra_dbs END_OF_DBS $options";

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
  my $result = $gomo->launchJob($req);

  # Give user the verification form and email message
  my $verify = make_verification();
  $utils->verify_opal_job($result, $address, $email_contact, $verify);

} # submit_to_opal

