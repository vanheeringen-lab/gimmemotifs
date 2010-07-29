#!@WHICHPERL@
##
## tomtom.cgi made from tomtom.pl by make
# Author: Timothy Bailey

# Combination CGI to:
# 1) Create the tomtom input form
# 2) Parse input form and execute query on local host

# defaults
$PROGRAM = "TOMTOM";
$NERRORS = 0;                            # no errors yet

use lib qw(@PERLLIBDIR@);
use Validation;
require "Utils.pm";
eval { require CGI; import CGI qw/:standard/; }; whine("$@") if $@;

$dir = "@MEME_DIR@";		# installed directory
$db_dir = "$dir/db/motif_databases";       # directory containing dbs
$bin_dir = "$dir/bin";		# directory for executables
$out_dir = "../output";         # directory for results
$names_file = "motif_db.csv";
$dbnames_url = "@SITE_URL@/cgi-bin/get_db_list.cgi?db_names=$names_file";
$db_doc_url = $dbnames_url . "&doc=1";
$REFRESH = 20;
$email_contact = '@contact@';

# get the parameters for the query
get_params();

# if there is no action specified, print an input form
if (! $action) {
  print_form() unless $NERRORS;
} else {
  run_query() unless $NERRORS;
  system("chmod 0777 *"); # make all files deletable
}

################################################################################
# Subroutines
################################################################################

#
# get parameters from the input
#
sub get_params {
  # retrieve the fields from the form
  $action = param('target_action');
  $address = param('address');
  $address_verify = param('address_verify');
  $query = param('query');
  $query_name = param('query_name');
  $example_query = param('example_query');
  $description = param('description');
  @database = split(/,/, param('target_db'));
  $dist = param('dist');
  $qthresh = param('-q-thresh');
} # get_params

#
# print the tomtom input form
# 
# Note: if query was included when this was called, that part of form
# is omitted.
#
sub print_form 
{
  my $action = "tomtom.cgi";
  my $logo = "../doc/images/tomtom_logo.png";
  my $alt = "$PROGRAM logo";
  my $form_description = qq {  
<B>
Use this form to use <A HREF="../tomtom-intro.html"><CODE>TOMTOM</CODE></A>
to compare your DNA motif against a database of known motifs
(e.g., JASPAR or Transfac).
</B>
<BR>
<CODE>TOMTOM</CODE> will rank the motifs in the target
database by the <I>q</I>-value of the similarity score.
<CODE>TOMTOM</CODE> outputs the <I>q</I>-value of
each match, and an alignment of LOGOs representing your query
and the target motif it matches.  
You can click on the target motif ID in your results to get more information.
(A Transfac license is required for viewing Transfac details.)
  }; # end quote

  #
  # required section
  #

  # See if the query has been input already
  # Put it in a hidden field if it has.
  my $req_left, $req_right;
  if ($query) {
    $req_left = qq {
<H3>
TOMTOM will compare your uploaded motif 
with each motif in a DNA <B>motif database</B>.
</H3>
<INPUT TYPE="hidden" NAME="query_name" VALUE="$query_name">
<INPUT TYPE="hidden" NAME="query" VALUE="
$query
">
    }; # end quote

  } else {

    $req_left = qq {
You can search a <B>single DNA motif</B> against
a DNA motif database.  The format for your motif should be as shown
on the right, with each column showing the counts of bases A, C, G and T
in that position of the motif.
<BR>
<BR>
<B>Paste your DNA motif here:</B>
<BR>
<TEXTAREA NAME="query" ROWS="4" COLS="50"></TEXTAREA>
<BR>
    }; # end quote

    $req_right = qq {
<table border="1">
<B>
<caption>Example Motif</caption>
<tr>
<td>
<pre>
 3  3 19  0  1  0  2 26  5
 8  0  0  1  0  1 23  1 15
14  0  9 27 26  4  3  0  4
 3 25  0  0  1 23  0  1  4
</pre>
<input type="hidden" name="example_query" value="
 3  3 19  0  1  0  2 26  5
 8  0  0  1  0  1 23  1 15
14  0  9 27 26  4  3  0  4
 3 25  0  0  1 23  0  1  4
">
</B>
</td>
</tr>
</table>
<BR>
<input type="submit" name="target_action" value="Search with Example Motif">
    }; # end quote

  }

  $req_left .= make_supported_databases_field("target_db", $dbnames_url, $db_doc_url, 1);

  $req_left .= qq {
<BR>
<BR>
<b>Select the Motif Column Comparison Function:</b>
<BR>
<!-- <Select> -->
<select name="dist">
<option selected value="pearson">Pearson correlation coefficient (pearson)
<option value="ed">Euclidean distance (ed)
<option value="sandelin">Sandelin-Wasserman similarity function (sandelin)
<!-- 
<option value="allr">Average log-likelihood ratio (allr)
<option value="chi">Pearson chi square test (chi)
<option value="fish">Fisher-Irwin exact test (fisher)
<option value="kullback">Kullback-Leibler divergence (kullback)
-->
</select>
<BR>
<BR>
<!-- q-value -->
<b>Significance threshold:</b> <i>q</i>-value &le
<INPUT CLASS="maininput" TYPE="TEXT" SIZE=2 NAME="-q-thresh" VALUE=0.5>
  }; # end quote

  $req_right .= qq {
<BR>
<BR>
<TABLE>
<B>
<CAPTION>Example Alignment</CAPTION>
</B>
<TR>
<TD>
<IMG SRC="../images/logo_alignment.png" WIDTH="200" BORDER="1">
</TD>
</TR>
</TABLE>
  }; # end quote

  my $required = make_input_table("Required", $req_left, $req_right);

  # add the required fields that are not used here but checked for
  # by the javascript
  $required .= qq {
<INPUT TYPE="hidden" name="address" value="none">
<INPUT TYPE="hidden" name="address_verify" value="none">
  }; # end quote

  #
  # optional arguments
  #

  my $descr = "motif";
  my $opt_left = make_description_field($descr, $description);

  my $opt_right = "";

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
# run a tomtom query nice-d
#
sub run_query {

  #umask 0777;

  # create a temporary directory name
  srand( time() ^ ($$ + ($$ << 15)) );
  my $random1 = int(rand(1) * 100000000);
  my $random2 = int(rand(1) * 100000000);
  $res_dir = "$out_dir/tomtom_${random1}_${random2}";
  mkdir "$res_dir";
  chmod 0777, "$res_dir";

  # change to the output directory
  chdir($res_dir);

  # create a link to tomtom
  symlink("$bin_dir/tomtom", "tomtom");
  symlink("$bin_dir/motif2meme.pm", "motif2meme.pm");
  # create a link to the database
  symlink("$db_dir/".$database[1], $database[1]);

  $final_out = "tomtom_web.html";

  open OUT, ">$final_out";
  print OUT <<EOF
<HTML>
<HEAD><TITLE>$PROGRAM Job Status Report</TITLE>
<META HTTP-EQUIV="expires" CONTENT="Sun 1 Jan 1997 00:00:00 GMT">
<META HTTP-EQUIV="Refresh" CONTENT="$REFRESH">
</HEAD>
<h2>This page will be updated every $REFRESH seconds
and will be replaced with your results automatically.
<p>
Your search may take a few moments.
</p>
You may bookmark this page and return to it later.
</h2>
</BODY></HTML>
EOF
    ;
  close OUT;

  # set up the sample query
  if ($action eq "Search with Example Motif") {
    $query = $example_query;
    $query_name = 'example_query';
  }

  # set up the query from the input box
  unless (defined $query_name) {
    $query_name = 'query';	# jaspar2meme won't work with just .pfm as filename
  }
  
  # create jaspar formatted query
  $jaspar = "$query_name.pfm";
  if (! open(OUT, ">$jaspar")) {
    whine("Unable to create query.");
    return(1);
  }

  # check query
  unless ($query =~ /\S/) {
    whine("You must enter a motif.");
    return(1);
  }

  print OUT $query;
  close OUT;
  chmod 0777, "$jaspar";

  # Make meme formatted query
  $meme = "query.meme";
  $err = "jaspar2meme.err";
  $description =~ s/\s+/_/g;			# needed due to MEME BL format restrictions
  $status = system("$bin_dir/jaspar2meme -pseudo 1e-6 -pfm -descr \'$description\' . 1> $meme 2> $err");
  chmod 0777, "$meme";
  $error = `cat $err`; unlink "$err";
  if ($status) {
    &whine("An error occurred converting your motif to the proper format.<BR>
      JASPAR2MEME returned: <pre>$error</pre>"
    );
    return(1);
  }

  # set up the target db parameters
  $target_db_type = $database[0];
  unless ($target_db_type) {
    whine("You must choose a motif database to search."); 
    return(1);
  }
  $target_db = "-target $database[1]";
  if ($database[6]) {
    $target_db .= " -target-url $database[6]";
  }

  # check that q-value threshold is legal
  if ($qthresh <= 0) {
    whine("You must choose a <i>q</i>-value threshold greater than 0.");
    return(1);
  }
  if ($qthresh > 1) {
    whine("You must choose a <i>q</i>-value threshold no greater than 1.");
    return(1);
  }

  # print the redirection
  print <<EOF
Content-Type: text/html

<html>
<head>
<title>
TOMTOM - Job Status Report
</title>
<meta name="author" content="Timothy Bailey">
</head>
<body>
EOF
    ;

  print "<center><h2>TOMTOM Job Status Report</h2></center>\n";
  print "<h2> Your results will be <a href='$res_dir/$final_out'> here</a>.<h2>\n";
  print "<p><h2>Running...</h2></p>\n";

  # run tomtom
  $err = "tomtom.err";
  $results = "tomtom.html";
  $status = system("nice -n +19 $bin_dir/tomtom -oc . -query $meme $target_db -dist $dist -q-thresh $qthresh -no-ssc 2> $err");
  print "<h2>Done.</h2>\n";
  $error = `cat $err`; unlink "$err";
  if ($status) {
    print "<BR>An error occurred running TOMTOM.<BR>
      TOMTOM returned: <pre>$error</pre>"
  } else {
    if (! rename ($results, $final_out)) {
      print "<BR>An error occurred in renaming output file.<BR>"
    }
  }
  print " </body> </html>\n";

} # run_query
