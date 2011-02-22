#!@WHICHPERL@
# Author: Timothy Bailey
##
## $Id: meme.pl 5280 2011-01-05 11:10:34Z james_johnson $
##
## $Log:$
## $Rev: 5280 $ $Date: 2011-01-05 03:10:34 -0800 (Wed, 05 Jan 2011) $ $Author: james_johnson $
## new meme form-submission file, which integrates with opal
##
##

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

our $upload_negfile_name; #why is this decared as 'our'?

# get the directories using the new installation scheme
my $dir = "@MEME_DIR@";		# installed directory
my $logs = "$dir/LOGS";		# directory for temp files
my $bin = "$dir/bin";		# directory for executables 
my $service_url = "@OPAL@/MEME_@S_VERSION@";
my $email_contact = '@contact@';
my $short_dna_only = 1;            # only short sequence bfiles
my $names_file = 'fasta_db.csv';
my $dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file&short_only=$short_dna_only";
my $db_doc_url = $dbnames_url . "&doc=1";

# defaults
my $PROGRAM = "MEME";
my $utils = new MemeWebUtils($PROGRAM, $bin);

#params
my ($options, $dna_options, $action, $address, $description, $datafile_name, 
  $data, $dist, $nmotifs, $minsites, $maxsites, $minw, $maxw, 
  $update_bfile_name, $evt, $shuffle, $pal, $posonly, $purgeoptions, 
  $dustoptions, $psp_type);
my @bfile;
my $upload_bfile_name;
my ($fasta_data, $alphabet, $num, $min, $max, $ave, $total);
my $bfile_data;
my $args;
my ($neg_data, $neg_alphabet, $neg_num, $neg_min, $neg_max, $neg_ave, $neg_total);
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

###############################################################################
####		SUBROUTINES:
###############################################################################

#
# print the form
#
sub print_form {

  my $action = "meme.cgi";
  my $logo = "../images/meme.png";
  my $alt = "$PROGRAM logo";
  my $form_description = qq {
Use this form to submit DNA or protein sequences to $PROGRAM. 
$PROGRAM will analyze your sequences for similarities among them and 
produce a description (<A HREF="../meme-intro.html"><B>motif</B></A>) 
for each pattern it discovers.
  }; # end quote

  #
  # required left side: address and sequences fields
  #
  my $seq_doc = "../help_sequences.html#sequences";
  my $alpha_doc = "../help_alphabet.html";
  my $format_doc = "../help_format.html";
  my $filename_doc = "../help_sequences.html#filename";
  my $paste_doc = "../help_sequences.html#actual-sequences";
  my $sample_file = "../examples/At.fa";
  my $sample_alphabet = "Protein";
  my $address_verify = undef; #TODO why isn't this set anywhere?
  my $target = undef; #TODO investigate what this is for and why it isn't declared!!!!
  my $req_left = $utils->make_address_field($address, $address_verify);
  $req_left .= $utils->make_upload_sequences_field("datafile", "data", $MAXDATASET,
    $seq_doc, $alpha_doc, $format_doc, $filename_doc, $paste_doc, $sample_file,
    $sample_alphabet, $target);

  #
  # required right side
  #

  # make distribution buttons
  my $distr_doc = "../meme-input.html#distribution";
  my $distr = "How do you think the occurrences of a single motif are <A HREF='$distr_doc'><B>distributed</B></A> among the sequences?\n<BR>\n";
  my $checked = "zoops";
  my @values = (
    "<B>One per sequence</B>\n<BR>,oops", 
    "<B>Zero or one</B> per sequence\n<BR>,zoops",
    "<B>Any number</B> of repetitions\n<BR>,tcm"
  );
  my $req_right = $utils->make_radio_field($distr, "dist", $checked, \@values);
  #$req_right .= "<B>Note:</B>The maximum number of occurrences of a motif is limited to $MAXSITES.";

  # make width fields
  $req_right .= qq {
<BR>
<!-- width fields -->
<BR>
$PROGRAM will find the optimum <A HREF="../help_width.html#motif_width"><B>width</B></A> 
of each motif within the limits you specify here: 
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="minw" VALUE=6>
<B>Minimum</B> width (>= $MINW)
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="maxw" VALUE=50>
<B><A HREF="../help_width.html#max_width">Maximum</B></A> width (<= $MAXW)
<BR>
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=2 NAME="nmotifs" VALUE=3>
Maximum <A HREF="../meme-input.html#nmotifs"><B>number of motifs</B></A>
to find
  }; # end quote

  # finish required fields
  my $required = $utils->make_input_table("Required", $req_left, $req_right);

  #
  # optional fields
  #

  # optional left: description, nsites, shuffle fields 
  my $descr = "sequences";
  my $opt_left = $utils->make_description_field($descr, $description);

  # make nsites fields
  $opt_left .= qq {
<BR>
<!-- nsites fields -->
<BR>
$PROGRAM will find the optimum 
<A HREF="../meme-input.html#nsites"><B>number of sites</B></A>
for each motif within the limits you specify here: 
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="minsites" VALUE="">
<B>Minimum</B> sites (>= $MINSITES)
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="maxsites" VALUE="">
<B>Maximum</B> sites (<= $MAXSITES)
<BR>
<BR>
  }; # end quote
  # shuffle checkbox
  my $doc = "../help_filters.html#shuffle";
  my $text = "<A HREF='$doc'><B>Shuffle</B></A> sequence letters";
  $opt_left .= $utils->make_checkbox("shuffle", 1, $text, 0);

  # optional right fields: bfile, shuffle, dna-only: strand, pal, PSP

# FIXME: started to implement background model list
  #my $opt_right = "Background model for sequences--select <B>one</B> of the following:<BR>\n";
  #%my $opt_right .= make_supported_databases_field("bfile", $dbnames_url, $db_doc_url);
# add the options for creating a PSP including uploading a negative set file

  my $opt_right .= qq {
<font style="background-color: yellow">NEW</font>
Perform <A HREF="../help_psp.html#discriminative"><B>discriminative</B></A>
motif discovery &ndash; Enter the name of a file containing
&lsquo;<A HREF="../help_psp.html#negativeset"><B>negative sequences</B>&rsquo;</A>:
<BR>
<div id="upload_neg_div">
  <input class="maininput" name="upload_negfile" type="file" size=18>
  <a onclick="clearFileInputField('upload_neg_div')" href="javascript:noAction();"><b>Clear</b></a>
</div>
<br><br>
Enter the name of a file containing a 
<A HREF="../meme-input.html#bfile"><B>background Markov model</B></A>:
<BR>
<div id="upload_bfile_div">
  <input class="maininput" name="upload_bfile" type="file" size=18>
  <a onclick="clearFileInputField('upload_bfile_div')" href="javascript:noAction();"><b>Clear</b></a>
</div>
  }; # end quote

  my $psp_doc = "../help_psp.html#psp";

  my $triples_doc =  "../help_psp.html#triples";

  # DNA-ONLY options (there is no definition of alphabet anywhere that I can find hence this will always run)
  #if (!defined $alphabet || $alphabet eq "ACGT") 
  {
    my ($doc, $text, $options);
    $doc = "../meme-input.html#posonly";
    $text = "Search given <A HREF='$doc' <B>strand</B></A> only";
    $options = $utils->make_checkbox("posonly", "1", $text, 0); 
    $options .= "<BR>\n";
    $doc = "../meme-input.html#pal";
    $text = "Look for <A HREF='$doc'><B>palindromes</B></A> only";
    $options .= $utils->make_checkbox("pal", "-pal", $text, 0); 
    $opt_right .= $utils->make_dna_only($options);
  } # dna options

  my $optional = $utils->make_input_table("Options", $opt_left, $opt_right, undef, undef); # add 1 at end to make visibility toggle
  
  # filtering options
  my $filter_purge = "";
  my $filter_dust = "";

  $filter_purge = qq {
<!-- purge fields -->
To use purge to remove similar sequences specify the maximum 
<A HREF="http://folk.uio.no/einarro/Presentations/blosum62.html"><B>blosum62
relatedness score</B></A> (100-200 recommended): 
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="purgescore" VALUE="">
<BR>
<BR>     
  }; # end quote

  $filter_dust = qq {
<!-- purge fields -->
To use <A HREF="http://www.liebertonline.com/doi/abs/10.1089/cmb.2006.13.1028"><B>dust</B></A>
to filter low complexity regions from your sequences set the cut-off
 (10-20 recommended): 
<BR>
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="dustcutoff" VALUE="">
<BR>
<BR>     
  }; # end quote

  # add the dust and purge input area
  my $filters = $utils->make_input_table("Filters", $filter_purge, $filter_dust);
  

  #
  # print final form
  #
  my $form = $utils->make_submission_form(
    $utils->make_form_header($PROGRAM, "Submission form"),
    $utils->make_submission_form_top($action, $logo, $alt, $form_description),
    $required, 
    #$optional.$filters, 
    $optional, 
    $utils->make_submit_button("Start search", $email_contact),
    $utils->make_submission_form_bottom(),
    $utils->make_submission_form_tailer()
  );
  print "Content-Type: text/html\n\n$form";
  
} # print_form

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
  $description = param('description');
  $datafile_name = param('datafile');
  $data = param('data');
  $dist = param('dist'); $options .= " -mod $dist" if $dist;
  $nmotifs = param('nmotifs'); $options .= " -nmotifs $nmotifs" if $nmotifs;
  $minsites = param('minsites'); $options .= " -minsites $minsites" if $minsites;
  $maxsites = param('maxsites'); $options .= " -maxsites $maxsites" if $maxsites;
  $minw = param('minw'); $options .=" -minw $minw" if $minw;
  $maxw = param('maxw'); $options .=" -maxw $maxw" if $maxw;
  @bfile = (param('bfile') ? split(/,/, param('bfile')) : ());
  $upload_bfile_name = param('upload_bfile'); 
  $options .= " -bfile uploaded_bfile" if $upload_bfile_name; 
  $evt = param('evt'); $options .= " -evt $evt" if $evt;
  $shuffle = param('shuffle'); $shuffle = 0 unless ($shuffle);
  $pal = param('pal'); $dna_options .= " -pal" if $pal;
  $posonly = param('posonly'); $dna_options .= " -revcomp" unless $posonly;

  # extra options
  $options .= " -time $MAXTIME -maxsize $MAXDATASET";

  # filter options
  $purgeoptions = param('purgescore');
  $dustoptions = param('dustcutoff');

  # PSP options
  $upload_negfile_name = param('upload_negfile');
  if ($upload_negfile_name) {
      $psp_type = param('psp');
  }

} # get_params

#
# Check the parameters on the form.
#
sub check_params {

  # change working directory to LOGS
  chdir($logs) || $utils->whine("Can't cd to $logs");

  # check that valid email address was provided
  $utils->check_address($address, $email_contact);

  # check description field
  $utils->check_description($description);

  # get FASTA sequences and information about them
  ($fasta_data, $alphabet, $num, $min, $max, $ave, $total)
    = $utils->get_sequence_data($data, $datafile_name, $MAXDATASET, $shuffle, $purgeoptions,
       $dustoptions);
  $alphabet = lc($alphabet);

  # get the background model file
  if ($upload_bfile_name) {
    $bfile_data = "";
    while (<$upload_bfile_name>) {
      s/\r\n/\n/g;                                # Windows -> UNIX eol
      s/\r/\n/g;                                  # MacOS -> UNIX eol
      $bfile_data .= $_;
    }
    if (! $bfile_data) { 
      $utils->whine("The background file '$upload_bfile_name' is empty or non-existent."); 
    }
  }

  # get the negative sequences file if specified; since there is no text box
  # for negative files, the first argument to get_sequence_data is undefined
  # options that don't apply left as undefined
  my $psp_args = ""; # want this empty if no negative sequences supplied
  my $psp_meme_options = "";
  if ($upload_negfile_name) {
    $utils->whine("Any number of repetitions not allowed ".
	      "for discriminative motif discovery") if($dist eq "tcm");
    ($neg_data, $neg_alphabet, $neg_num, $neg_min, $neg_max, $neg_ave, $neg_total)
	      =  $utils->get_sequence_data(undef, $upload_negfile_name, $MAXDATASET*4);
    $utils->whine ("negative sequence alphabet ($neg_alphabet) should be same".
	      " as postive sequence alphabet ($alphabet)")
	      unless lc($neg_alphabet) eq lc($alphabet);
    $neg_alphabet = "prot" if lc($neg_alphabet) eq "protein";
    # use the minw and maxw settings for MEME for finding the PSP but
    # trim to the allowed range for PSPs
    # the actual width set by the PSP finder is that with the highest
    # score before normalizing; allow X or N or other nonspecific residue/base
    # codes (but score any sites containing them as zero)
    my $psp_minw = $minw < $MINPSPW ? $MINPSPW : $minw;
    my $psp_maxw = $maxw > $MAXPSPW ? $MAXPSPW : $maxw;
    $utils->whine(
       "min width ($psp_minw) must be <= max width ($psp_maxw) for PSPs"
    ) unless ($psp_minw <= $psp_maxw);
    $psp_args = "-neg negfile -alpha $neg_alphabet -minw $psp_minw -maxw $psp_maxw -arbitraryOK";
    $psp_args .= ' -maxrange' if ($neg_alphabet eq "prot");
    $psp_args .= ' -triples' if ($neg_alphabet eq "prot");
    $psp_meme_options = " -psp priors.psp";
  }

  #FIXME tlb; started to implement bfile menu
  # get the chosen background model name (from menu)
  # my $short_dna_only = 1;     #  short sequence dbs only
  # ($db, $bfile_data, $db_alphabet, $target_alphabet, $db_root, $prot_db, $nuc_db, $short_seqs, $bfile_menu_name, $bfile_long_name, $bfile_descr) =
  #   get_db_name($alphabet, $upload_bfile_name, $MAX_UPLOAD_SIZE, 0, $short_dna_only, @bfile);

  # check that TCM specified if only one sequence
  if (($num == 1) && !($dist eq "tcm")) {
    $utils->whine("
      You must specify <I>Any number of repetitions</I> under the 
      <I>distribution</I> option
      since your dataset contains only one sequence.  Alternatively, you might 
      wish to break your sequence into several sequences.
    ");
  }

  # check that number of motifs is OK
  if (not defined $nmotifs) {
    $utils->whine(
      "You must define the <i>number of motifs</i>."
    );
  } elsif ($nmotifs < 1 || $nmotifs > $MAXMOTIFS) {
    $utils->whine(
      "You must specify <I> number of motifs</I> between 1 and $MAXMOTIFS, inclusive."
    );
  }

  # check that number of sites is OK
  &check_nsites($minsites, $maxsites, $dist, $num);

  # check that width is OK
  &check_width;

  # remove spaces, non-ASCII and single quotes from $datafile_name
  $datafile_name = "pasted sequences" unless($datafile_name); 
  $datafile_name =~ s/[ \'\x80-\xFF]/\_/g;

  # add dna options if alphabet permits
  $options .= $dna_options if ($alphabet eq "dna");

  # add PSP options if negative file specified
  $options .= $psp_meme_options if ($upload_negfile_name);

  # finish MEME argument list
  $psp_args .= " -pos sequences" if ($psp_args);
  $args = $psp_args ." meme sequences -sf $datafile_name -$alphabet $options";
 
} # check_params

#
# check that number of sites is OK
#
sub check_nsites {
  my (
    $minsites, 			# minimum nsites
    $maxsites, 			# maximum nsites
    $dist, 			# type of distribution
    $num			# number of sequences
  ) = @_;

  if ($minsites && $minsites < $MINSITES) {
    $utils->whine("You must specify a minimum number of sites >= $MINSITES");
  }
  if ($maxsites && $maxsites > $MAXSITES) {
    $utils->whine("You must specify a maximum number of sites <= $MAXSITES");
  }
  if (($minsites && $maxsites) && $minsites > $maxsites) {
    $utils->whine(
     "The minimum number of sites is larger than the maximum number of sites");
  }
  # Make sure there aren't too many sequences in OOPS mode 
  # Note: MAXSITES and other globals are in lib/perl/Globals.pm
  if ($dist eq "oops" && $num > $MAXSITES) {
    $utils->whine(
      "Your dataset must contain no more than $MAXSITES sequences when you specify",
      "that the motif is <I>distributed</I> one per sequence.",
      "Please input a dataset with no more than $MAXSITES sequences or chose a",
      "different motif distribution."
    );
  }
} # check_nsites

#
# check that width is OK
#
sub check_width {

  if ($minw < $MINW) {
    $utils->whine("The minimum width you specified ($minw) is too small.<BR>")
  }
  if ($maxw > $MAXW) {
    $utils->whine("The maximum width you specified ($maxw) is too large.<BR>")
  }
  if ($minw > $maxw) {
    $utils->whine("The minimum width you specified ($minw) is larger than the
	maximum width you specified ($maxw).<BR>")
  }
} # check_width


#
# make the verification message in HTML
#
sub make_verification {

  my $content = "<HR><UL>\n";

  $content .= "<LI> Description:<B> $description </B>\n" if $description;
  $content .= "<LI> Sequence file:<B> $datafile_name </B>\n";

  $content .= "<LI> Distribution of motif occurrences: <B>";
  if ($dist eq "oops") {
    $content .= "One per sequence</B>\n";
  } elsif ($dist eq "zoops") {
    $content .= "Zero or one per sequence</B>\n";
  } else {
    $content .= "Any number of repetitions</B>\n";
  }

  $content .= "<LI> Number of different motifs:<B> $nmotifs </B>\n" if $nmotifs;
  $content .= "<LI> Minimum number of sites:<B> $minsites </B>\n" if $minsites;
  $content .= "<LI> Maximum number of sites:<B> $maxsites </B>\n" if $maxsites;
  $content .= "<LI> Minimum motif width:<B> $minw </B>\n" if $minw;
  $content .= "<LI> Maximum motif width:<B> $maxw </B>\n" if $maxw;
  $content .= "<LI> Searching: <B>given strand only</B>\n" if $posonly;
  $content .= "<LI> Looking for: <B>palindromes only</B>\n" if $pal;
  $content .= "<LI> <B>Shuffling</B> letters in input sequences\n" if $shuffle;
  #$content .= "<LI> Background model file: <B>$bfile_menu_name</B>\n" if $bfile_menu_name; #TODO not implemented
  $content .= "<LI> Background model file: <B>$upload_bfile_name</B>\n" if $upload_bfile_name;

  if ($upload_negfile_name) {
    $content .= "<LI> Negative sequence file: <B>$upload_negfile_name</B><BR>";
    $content .= " PSP options: <B>";
    $content .= "spaced triples, " if ($neg_alphabet eq "prot");
    $content .= "discriminative";
    $content .= " prior</B>";
  }

  $content .= "
    <LI> Statistics on your dataset:
      <TABLE BORDER>
	<TR> <TD> type of sequence <TH ALIGN=RIGHT> $alphabet
	<TR> <TD> number of sequences <TH ALIGN=RIGHT> $num
	<TR> <TD> shortest sequence (residues) <TH ALIGN=RIGHT> $min
	<TR> <TD> longest sequence (residues) <TH ALIGN=RIGHT> $max
	<TR> <TD> average sequence length (residues) <TH ALIGN=RIGHT> $ave
	<TR> <TD> total dataset size (residues) <TH ALIGN=RIGHT> $total";
  if ($upload_negfile_name) {
      $content .= "
        <TR> <TD> number of <I>negative</I> sequences <TH ALIGN=RIGHT> $neg_num
	<TR> <TD> shortest <I>negative</I> sequence (residues) <TH ALIGN=RIGHT> $neg_min
	<TR> <TD> longest <I>negative</I> sequence (residues) <TH ALIGN=RIGHT> $neg_max
	<TR> <TD> average <I>negative</I> sequence length (residues) <TH ALIGN=RIGHT> $neg_ave
	<TR> <TD> total <I>negative</I> dataset size (residues) <TH ALIGN=RIGHT> $neg_total";
  }
  $content .= "
      </TABLE>
    </UL>
  ";

  return($content);

} # make_verification

#
# Submit job to webservice via OPAL
#
sub submit_to_opal
{
  my $service = OpalServices->new(service_url => $service_url);

  #
  # start OPAL requst
  #
  my $req = JobInputType->new();
  $req->setArgs($args);

  #
  # create list of OPAL file objects
  #
  my @infilelist = ();
  # 1) Sequence file
  my $inputfile = InputFileType->new("sequences", $fasta_data);
  push(@infilelist, $inputfile);
  # 2) Background file
  if ($upload_bfile_name) {
    $inputfile = InputFileType->new("uploaded_bfile", $bfile_data);
    push(@infilelist, $inputfile);
  }
  # 2.1) PSP negative set file
  if ($upload_negfile_name) {
    $inputfile = InputFileType->new("negfile", $neg_data);
    push(@infilelist, $inputfile);
  }
  # 3) Email address file (for logging purposes only)
  $inputfile = InputFileType->new("address_file", $address);
  push(@infilelist, $inputfile);
  # 4) Submit time file (for logging purposes only)
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
