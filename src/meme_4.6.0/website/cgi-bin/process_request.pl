#!@WHICHPERL@
##
## $Id: process_request.pl 3941 2009-08-10 06:06:11Z tbailey $
##
## $Log: process_request.pl,v $
# Author: Timothy Bailey

use lib qw(@PERLLIBDIR@);
#use Globals;
$dir = "@MEME_DIR@";		# installed directory
$bin_dir = "$dir/bin";		# directory for executables
$logs_dir = "$dir/logs";	# directory for executables

#
# Process a request submitted from a MEME results form.
# Requests include:
#	Submit to MAST
#	View hidden field
#	Submit to BLOCKS
#	Submit to Meta-MEME
#	Submit to Tomtom
#	Submit to FIMO
#	Submit to GOMO
#	Submit to GLAM2 and GLAM2SCAN
#	Print a man page
#
#

# requires
# use the CGI package
use CGI qw/:standard/;			# use the CGI package
use HTTP::Request::Common qw(POST);
use LWP::UserAgent;

# URLs used
$mast_url = '@SITE_URL@/cgi-bin/mast.cgi';
$tomtom_url = '@SITE_URL@/cgi-bin/tomtom.cgi';
$glam2_url = '@SITE_URL@/cgi-bin/glam2.cgi';
$glam2scan_url = '@SITE_URL@/cgi-bin/glam2scan.cgi';
$metameme_url = '@SITE_URL@/cgi-bin/mhmm_process_request.cgi';
$gomo_url = '@SITE_URL@/cgi-bin/gomo.cgi';
$fimo_url = '@SITE_URL@/cgi-bin/fimo.cgi';
$blocks_url = '@BLOCKS_URL@';

#
# get the action (type of submit button pressed) and branch on it
#
$action = param('action');			# which submit button pressed
$nmotifs = param('nmotifs');			# number of motifs
#print header;

if ($action eq "MAST") {
  search('MAST', 1, $nmotifs);
} elsif (
  $action =~ /^MAST PSSM (\d+)/ ||
  $action =~ /^Scan PSSM (\d+)/ 	# backwards compatability
  ) {
  search('MAST', $1, $1);
} elsif ($action eq "FIMO") {
  search('FIMO', 1, $nmotifs);
} elsif ($action =~ /^FIMO PSPM (\d+)/) {
  search('FIMO', $1, $1);
} elsif ($action eq "GOMO") {
  search('GOMO', 1, $nmotifs);
} elsif ($action =~ /^GOMO PSPM (\d+)/) {
  search('GOMO', $1, $1);
} elsif ($action =~ /^BLOCKS/) {
  submit_block(1..$nmotifs);
} elsif ($action =~ /^Submit BLOCK (\d+)/) {
  submit_block($1);
} elsif ($action =~ /^View (\w+) (\d+)/) {
  hidden_field($1, $2);
} elsif ($action =~ /^Submit BLOCK (\d+)/) {
  submit_block($1);
} elsif (
  $action =~ /^TOMTOM PSPM (\d+)/ ||
  $action =~ /^Compare PSPM (\d+)/ ||	# backwards compatability
  $action =~ /^COMPARE PSPM (\d+)/	# backwards compatability
  ) {
  submit_pspm_to_tomtom($1);
} elsif ($action =~ /^View motif summary/) {
  hidden_field('motif-summary');
} elsif ($action =~ /MetaMEME/){
  forward_query("$metameme_url");
} elsif ($action =~ /re-run GLAM2/){
  $params{'-n'} = 2 * param('-n'); 	# double number of iterations
  $params{'-n_attributes'} = "READONLY";	# don't allow -n to be changed by user
  forward_query("$glam2_url")
} elsif ($action =~ /^Scan alignment (\d+)/) {
  $params{'inline_name'} = "alignment_$1";
  $params{'inline_aln'} = param('aln'.$1);
  forward_query("$glam2scan_url")
} elsif ($action =~ /MEME Man Page/){
  print_man("meme");
} elsif ($action =~ /MAST Man Page/){
  print_man("mast");
} else {
  print_header($action);			# print a response header
  print "action: $action\n";
  print "<H1>Feature not implemented yet.</H1>\n";
  exit(1);
}

# all done!
exit(0);

#
# forward a query to another service
#
sub forward_query {
  my($cgi) = @_;
  print "Content-type: text/html\n\n";
  # Get the contents of the query
  foreach $pname (param()) { $params{$pname} = param($pname) unless $params{$pname};}
  undef $params{'action'};
  # post the query
  $ua = LWP::UserAgent->new();
  my $req = POST "$cgi", [%params];
  my $request = $ua->request($req);
  my $content = $request->content;
  print $content;
} # forward_query

#
# Submit a search query
#
sub search {
  my($program, $start, $stop) = @_;
  my(%params);

  print "Content-type: text/html\n\n";		# start form

  # set up url
  if ($program eq 'FIMO') {
    $url = $fimo_url;
  } elsif ($program eq 'GOMO') {
    $url = $gomo_url;
  } elsif ($program eq 'MAST') {
    $url = $mast_url;
  } elsif ($program eq 'MAST') {
  } else {
    print_header("$program: Search not implemented yet.");
    exit(1);
  } 

  # check input
  if ($stop - $start < 0) {
    print_header("$program Search");
      print "<html><body><PRE>There are no valid motifs in your MEME file.</PRE></body></html>\n";
    exit(1);
  }

  # set up the parameters to pass
  $params{'version'} = param('version');
  $params{'alphabet'} = param('alphabet');
  $params{'bgfreq'} = param('bgfreq');
  $params{'inline_name'} = param('name');
  $params{'inline_motifs'} = param('version') . "\n\nALPHABET= " . param('alphabet') . "\n\nstrands: + - \n\nBackground letter frequencies (from dataset with add-one prior applied):" . param('bgfreq') . "\n";
  for ($i=$start; $i<=$stop; $i++) { 
    $params{'inline_motifs'} .= "\nMOTIF ".$i." m" . $i . "\n";
    $params{'inline_motifs'} .= param('motifblock'.$i);
    $params{'inline_motifs'} .= param('pssm'.$i); 
    $params{'inline_motifs'} .= "\n"; 
    $params{'inline_motifs'} .= param('pspm'.$i);
  }
  $params{'address'} = $params{'address_verify'} = param('address');
  
  # post the query
  my $ua = LWP::UserAgent->new();
  my $req = POST "$url", [%params];
  my $request = $ua->request($req);
  my $content = $request->content;
  print $content;
} # search 

#
# view a hidden field
#
sub hidden_field {
  my($name, $number) = @_;

  print_header($action);
  if ($name eq "BLOCK") {
    $field = param('BLOCKS'.$number);		# get BLOCK
  } elsif ($name eq "FASTA") {
    $field = param('BLOCKS'.$number);
    $field = block2fasta($field);		# convert to FASTA
  } elsif ($name eq "RAW") {
    $field = param('BLOCKS'.$number);
    $field = block2raw($field);			# convert to raw format 
  } elsif ($name eq "PSSM") {
    $field = param('pssm'.$number);		# get PSSM
  } elsif ($name eq "PSPM") {
    $field = param('pspm'.$number);		# get PSPM
  } elsif ($name eq "alignment") {
    $field = param('aln'.$number);		# get PSPM
  } elsif ($name eq "motif-summary") {
    $field = param('motif-summary');		# get motif-summary
  } else { 
    print "Unknown hidden field type: $name\n";
    exit(1)
  }
  print "<PRE>$field</PRE>";
} # hidden_field

#
# Convert a BLOCK to RAW sequence format
#
sub block2raw {
  my($block) = @_;
  my($i, @lines, @words, $raw);

  $fasta = "";					# return value
  @lines = split(/\n/, $block); 		# split block into lines
  for ($i = 2; $i<$#lines; $i++) {
    last if $lines[$i] =~ /^\/\//;
    @words = split(/\s+/, $lines[$i]);		# split line into words
    # get sequence line
    $raw .= "$words[3]\n";
  }
  $raw;
} # block2raw

#
# Convert a BLOCK to FASTA
#
sub block2fasta {
  my($block) = @_;
  my($i, @lines, @words, $fasta, $start);

  $fasta = "";					# return value
  @lines = split(/\n/, $block); 		# split block into lines
  for ($i = 2; $i<$#lines; $i++) {
    last if $lines[$i] =~ /^\/\//;
    @words = split(/\s+/, $lines[$i]);		# split line into words
    # get id line and sequence line
    $start = substr($words[2], 0, length($words[2])-1);
    $fasta .= ">$words[0] ( start= $start )\n$words[3]\n";
  }
  $fasta;
} # block2fasta

#
# Submit a block to the blocks processor
#
sub submit_block {
  my(@numbers) = @_;
  my($blocks);

  print_header("Submit BLOCKS");
 
  # get the BLOCK(S)
  $blocks = "";
  foreach $number (@numbers) {
    $blocks .= param('BLOCKS'.$number);
  }

  $ua = LWP::UserAgent->new();
  my $req = POST $blocks_url, [ sequences => $blocks ];
  my $request = $ua->request($req);
  #$content = $ua->request($req)->as_string;
  $content = $request->content;
  # put in the absolute url's : this is FRAGILE!
  $content =~ s#HREF=\"#HREF=\"http://blocks.fhcrc.org#g;
  $content =~ s#ACTION=\"/blocks-bin#ACTION=\"http://blocks.fhcrc.org/blocks-bin#g;
  print $content;
} # submit_block

#
# Compare a PSPM to the TOMTOM motif database search tool
#
sub submit_pspm_to_tomtom {
  my($number) = @_;
  my($pspm, @fields, $i, $n, $w, $row, $col);

  print "Content-type: text/html\n\n";

  # get the motif PSPM 
  $_ = param('pspm'.$number);
  @fields = split;

  # ignore 1st 10 entries
  $nsites = $fields[7];               # save number of sites
  for ($i=0; $i<=$#fields-10; $i++) { $fields[$i] = $fields[$i+10] * $nsites; }
  $#fields = $#fields - 10;

  # rotate PSPM 90 degrees (natural format) as a string with newlines
  $n = $#fields + 1;		# number of entries in motif
  $w = $n/4;			# motif width
  $pspm = "";
  for ($row=0; $row<4; $row++) {
    for ($col=0; $col<$w; $col++) {
      $pspm .= " " . $fields[($col*4) + $row];
    }
    $pspm .= "\n";		# terminate row with newline
  }

  # create the request
  $ua = LWP::UserAgent->new();
  my $req = POST "$tomtom_url",
    Content_Type => 'multipart/form-data',
    Content => [ 'query' => $pspm, 'query_name' => "Motif_$number" ];
  my $request = $ua->request($req);
  $content = $request->content;

  # display the page
  print $content;
} # submit_pspm_to_tomtom


#
# Print a man page
#
sub print_man {
  my($command) = @_;
  print <<END; 
  Content-type: text/plain

END
  chdir("$logs_dir");
  $bin = "$bin_dir";
  @tmp = `$bin/$command`;
  print @tmp[3..$#tmp];
}

#
# start a response form
#
sub print_header {
  my($action) = @_;
  print <<END; 
Content-type: text/html

<HTML>
<TITLE> MEME SUITE - $action</TITLE>
<BODY BACKGROUND=\"../images/bkg.jpg\">
END
} # print_header
