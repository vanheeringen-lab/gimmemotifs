#!@WHICHPERL@
# Author: Timothy Bailey
# Modified: Philip Machanick (for ChIP-seq)
##
## $Id: meme.pl 4921 2010-09-06 01:56:55Z phillip_machanick $
##
## $Log:$
## $Rev: 4921 $ $Date: 2010-09-06 11:56:55 +1000 (Mon, 06 Sep 2010) $ $Author: phillip_machanick $
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

use Scalar::Util qw(looks_like_number);

my $DNA_ALPHABET = "ACGTURYKMSWBDHVN"; #includes IUPAC ambiguous characters

# get the directories using the new installation scheme
my $dir = "@MEME_DIR@";		# installed directory
my $logs = "$dir/LOGS";		# directory for temp files
my $bin = "$dir/bin";		# directory for executables 
my $service_url = "@OPAL@/MEMECHIP_@S_VERSION@";
my $email_contact = '@contact@';
my $short_dna_only = 1;            # only short sequence bfiles
my $names_file = 'fasta_db.csv';
my $dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file&short_only=$short_dna_only";
my $db_doc_url = $dbnames_url . "&doc=1";

# defaults
my $PROGRAM = "MEME";
my $utils = new MemeWebUtils($PROGRAM."CHIP", $bin);

#params
my ($options, $dna_options, $action, $address, $description, $datafile_name, 
  $data, $dist, $nmotifs, $minsites, $maxsites, $minw, $maxw,
  $evt, $pal, $posonly);
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

  my $action = "meme-chip.cgi";
  my $logo = "../doc/images/memechip_logo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = qq {
This version of $PROGRAM is designed especially for discovering <A
HREF="../meme-intro.html"><B>motifs</B></A>
in <b>LARGE</b> sets of DNA sequences such as those produced by ChIP-seq
experiments. Unlike the MEME page, <b>there is no upper limit on the number
or size of your sequences.</b>
A motif is a description of a common pattern in sequences,
such as a transcription factor binding site.  Use this form to submit DNA
sequences to $PROGRAM and a number of other tools that aid in finding motifs
and evaluating the quality of the found motifs.
For more information about MEME-ChIP,
<span style="color:red; font-weight: bold">see this <a
href="../doc/meme-chip-tutorial.html">Tutorial</a></span>.
  }; # end quote

  #
  # required left side: address and sequences fields
  #
  my $seq_doc = "../help_sequences.html#sequences";
  my $alpha_doc = "../help_alphabet.html";
  my $format_doc = "../help_format.html";
  my $filename_doc = "../help_sequences.html#filename";
  my $paste_doc = "../help_sequences.html#actual-sequences";
  my $sample_file = "../examples/sample-dna-Klf1.fa";
  my $sample_alphabet = "DNA";
  my $address_verify = undef; #TODO why isn't this set anywhere?
  my $target = undef; #TODO investigate what this is for and why it isn't declared!!!!
  my $req_left = $utils->make_address_field($address, $address_verify);
  $req_left .= $utils->make_upload_sequences_field("datafile", "data", undef,
    $seq_doc, $alpha_doc, $format_doc, $filename_doc, $paste_doc, $sample_file,
    $sample_alphabet, $target);

  #
  # required right side
  #

  # make distribution buttons
  my $distr_doc = "../meme-input.html#distribution";
  my $distr = "How do you think the occurrences of a single motif are <A HREF='$distr_doc'><B>distributed</B></A> among the DNA sequences?\n<BR>\n";
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
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=3 NAME="maxw" VALUE=$DEFAULTMAXW>
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
  my $descr = "DNA sequences";
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

# FIXME: started to implement background model list
  #my $opt_right = "Background model for sequences--select <B>one</B> of the following:<BR>\n";
  #%my $opt_right .= make_supported_databases_field("bfile", $dbnames_url, $db_doc_url);
# removed options for creating a PSP including uploading a negative set file

  my $opt_right .= qq {
Enter the name of a file containing a 
<A HREF="../meme-input.html#bfile"><B>background Markov model</B></A>:
<BR>
<div id="upload_bfile_div">
  <input class="maininput" name="upload_bfile" type="file" size=18>
  <a onclick="clearFileInputField('upload_bfile_div')" href="javascript:noAction();"><b>Clear</b></a>
</div>
  }; # end quote

  # DNA-ONLY options (only do DNA here but leave comment in as placeholder in case this changes)
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
    $opt_right .= $options;
  } # dna options

  # last 1 turns on showhide
  my $showhide = 1;
  my $optional = $utils->make_input_table("Options", $opt_left, $opt_right, undef, undef,
					  ); # add $showhide if you want the options to toggle hidden-visible
  #
  # print final form
  #
  my $form = $utils->make_submission_form(
    $utils->make_form_header($PROGRAM."CHIP", "Submission form"),
    $utils->make_submission_form_top($action, $logo, $alt, $form_description),
    $required, 
    $optional , #. $advanced, 
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
  $pal = param('pal'); $dna_options .= " -pal" if $pal;
  $posonly = param('posonly'); $dna_options .= " -revcomp" unless $posonly;

  # extra options
  $options .= " -time $MAXTIME";

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
    = $utils->get_fasta_data($data, $datafile_name,
			     $DNA_ALPHABET, "DNA");

  $utils->whine("MEME cannot process sequences of length shorter than $MINMEMESEQW: fix your data and resubmit.")
		if $min <  $MINMEMESEQW;
  $alphabet = "dna"; # always DNA in this case

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

  # add dna options if alphabet permits (placeholder: only DNA here)
  $options .= $dna_options; 

  # finish backend argument list
  my $psp_args = ""; # placeholder: not doing PSP in this version
  # FIXME: replace array contents below by DB(s) from user interface
  my @databases = ("JASPAR_CORE_2009.meme");

  #################################################
  # sequences is the name of the file Opal will set up for the back
  # end script in submit_to_opal

  # Add command lines in the order they should run in the back end here.
  # General rule: |<name>| rewritten in back end, <name>: a variable whose
  # value is the rest of the arguments up to another "<..>" or <name>:
  # must rewrite vars in "||" such as |MEMESEQUENCES|, |AMESEQUENCES|, in the backend
  # TOMTOM uses the same DB as AME: if this changes, add the DB spec
  # onto the end of the TOMTOM args and FIX IN meme-chip back end script
  # FIXME: SpaMo taken out until a good version is available.
  $args = " <MEME> |MEMESEQUENCES| -sf $datafile_name -$alphabet $options".
      " <|UP|> ".
      " <TOMTOMMEME> -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc ".
      "|MEMEDIR|/meme.html |MOTIFDB|".
      " <MEMEMAST> |MEMEDIR|/meme.html sequences -ev |NSEQS|".
      " <AMA> --sdbg 0 |MEMEDIR|/meme.html sequences".
      " <|DOWN|>".
      " <DREME> -p |DREMESEQUENCES|".
      " <|UP|> ".
      " <TOMTOMDREME> -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc ".
      " |DREMEDIR|/dreme.txt |MOTIFDB|".
      " <DREMEMAST> |DREMEDIR|/dreme.txt sequences -ev |NSEQS| -mt DMASTTHRESH:0.0005".
#      " <SPAMO> -png -r 1 sequences |MEMEDIR|/meme.html ||". # same DB as AME FIXME: reinstate
      " <|DOWN|>".
      " <AME> --fix-partition |AMESIZE| --bgformat 0 |AMESEQUENCES|".
      " MOTIFDB:@databases".
      " <END>";
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
  #$content .= "<LI> Background model file: <B>$bfile_menu_name</B>\n" if $bfile_menu_name; #TODO not implemented
  $content .= "<LI> Background model file: <B>$upload_bfile_name</B>\n" if $upload_bfile_name;

  $content .= "
    <LI> Statistics on your dataset:
      <TABLE BORDER>
	<TR> <TD> type of sequence <TH ALIGN=RIGHT> $alphabet
	<TR> <TD> number of sequences <TH ALIGN=RIGHT> $num
	<TR> <TD> shortest sequence (bases) <TH ALIGN=RIGHT> $min
	<TR> <TD> longest sequence (bases) <TH ALIGN=RIGHT> $max
	<TR> <TD> average sequence length (bases) <TH ALIGN=RIGHT> $ave
	<TR> <TD> total dataset size (bases) <TH ALIGN=RIGHT> $total";
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
  # 3) Email address file (for logging purposes only)
  $inputfile = InputFileType->new("address_file", $address);
  push(@infilelist, $inputfile);
  # 4) Submit time file (for logging purposes only)
  $inputfile = InputFileType->new("submit_time_file", `date -u '+%d/%m/%y %H:%M:%S'`);
  push(@infilelist, $inputfile);

  # Add file objects to request
  $req->setInputFile(@infilelist);

  `echo about to launch > /tmp/verify`;

  # Submit the request to OPAL
  my $result = $service->launchJob($req);

  `echo about to verify >> /tmp/verify`;

  # Give user the verification form and email message
  my $verify = make_verification();

  print "Content-Type: text/html\n\n";
  $utils->verify_opal_job($result, $address, $email_contact, $verify);

} # submit_to_opal
