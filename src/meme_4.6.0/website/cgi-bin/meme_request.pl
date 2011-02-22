#!@WHICHPERL@

use strict;
use lib qw(@PERLLIBDIR@);
use CGI qw/:standard/;
use HTTP::Request::Common qw(POST);
use LWP::UserAgent;
use File::Temp qw/ tempfile /;

my ($blocks_url) = '@BLOCKS_URL@';
my ($fimo_url) = '@SITE_URL@/cgi-bin/fimo.cgi';
my ($gomo_url) = '@SITE_URL@/cgi-bin/gomo.cgi';
my ($mast_url) = '@SITE_URL@/cgi-bin/mast.cgi';
my ($tomtom_url) = '@SITE_URL@/cgi-bin/tomtom.cgi';


my $query = CGI->new();

# discover what the user wants to do
my ($param_name, $program, $motif);
foreach $param_name ($query->param()) {
  if ($param_name =~ m/^do_(MAST|FIMO|GOMO|BLOCKS|TOMTOM|LOGO)_(all|\d+)$/) {
    $program = $1;
    $motif = $2;
    last;
  }
}

# get the range of motifs
my ($nmotifs, $start, $end);
$nmotifs = int($query->param("nmotifs"));
if ($motif eq "all") {
  $start = 1;
  $end = $nmotifs;
} elsif ($motif) {
  $start = int($motif);
  $end = $start;
}

# generate content
my ($content, $content_len);
if ($nmotifs > 0) {
  if ($motif && $program) {
    if ($program eq "MAST") {
      $content = search($mast_url, $start, $end);
    } elsif ($program eq "FIMO") {
      $content = search($fimo_url, $start, $end);
    } elsif ($program eq "GOMO") {
      $content = search($gomo_url, $start, $end);
    } elsif ($program eq "BLOCKS") {
      $content = submit_block($blocks_url, $start, $end);
    } elsif ($program eq "TOMTOM") {
      $content = search($tomtom_url, $start, $end);
    } elsif ($program eq "LOGO") {
      #this produces an image which requires different headers so it calls exit and does not return
      generate_logo($start);
    } elsif ($program) {
      $content = error_page("Unknown program", "The program \"$program\" can't be handled by this script.");
    }
  } else {
    $content = error_page("No action?", "The script couldn't find a parameter of the form \"do_[Program]_[Motif#|all]\" and so can't perform an action."); 
  }
} else {
  $content = error_page("No motifs?", "The number of motifs parameter resolved to a non-positive number. Maybe it wasn't specified?");
}
$content_len = length($content);

# output content
print "Content-type: text/html", "\n";
print "Content-length: $content_len", "\n\n";
print $content;
exit(0);

sub get_motif {
  my ($start, $stop) = @_;
  my ($motif, $i);

  $motif = $query->param('version') . "\n\nALPHABET= " . $query->param('alphabet') . 
      "\n\nstrands: + - \n\nBackground letter frequencies (from dataset with add-one prior applied):\n" . 
      $query->param('bgfreq') . "\n";
  for ($i = $start; $i <= $stop; $i++) { 
    $motif .= "\nMOTIF ".$i." m" . $i . "\n";
    $motif .= $query->param('motifblock'.$i);
    $motif .= $query->param('pssm'.$i); 
    $motif .= "\n"; 
    $motif .= $query->param('pspm'.$i);
  }
  $motif;
}

#
# Test if all parameters required for a search are avaliable
#
sub can_search {
  my ($start, $stop) = @_;
  my ($all_params) = $query->param('version') && $query->param('alphabet') 
      && $query->param('bgfreq') && $query->param('name') && $query->param('bgfreq');
  if ($all_params) {
    my ($i);
    for ($i = $start; $i <= $stop; $i++) {
      $all_params = $query->param('motifblock'.$i) && $query->param('pssm'.$i) && $query->param('pspm'.$i);
      last unless $all_params;
    }
  }
  $all_params;
}

#
# Submit a search query
#
sub search {
  my($url, $start, $stop) = @_;
  my(%params, $content, $i);

  if (can_search($start, $stop)) {
    # set up the parameters to pass
    $params{'version'} = $query->param('version');
    $params{'alphabet'} = $query->param('alphabet');
    $params{'bgfreq'} = $query->param('bgfreq');
    $params{'inline_name'} = $query->param('name');
    $params{'inline_motifs'} = get_motif($start, $stop);
    
    # post the query
    my($ua) = LWP::UserAgent->new();
    my($req) = POST "$url", [%params];
    my($request) = $ua->request($req);
    $content = $request->content;
  } else {
    $content = error_page("Missing required variable","");
  }
  $content;
} # search 

sub error_page {
  my($title, $message) = @_;
  my($content) = 
      "<html>\n".
      "\t<head>\n".
      "\t\t<title>$title</title>\n".
      "\t</head>\n".
      "\t<body>\n".
      "\t\t<h1>$title</h1>\n".
      "\t\t<p>$message</p>\n".
      "\t\t<table border=\"1\">\n".
      "\t\t\t<tr>\n".
      "\t\t\t\t<th>Parameter Name</th>\n".
      "\t\t\t\t<th>Parameter Value</th>\n".
      "\t\t\t</tr>\n";
  my ($param_name, $param_value);
  foreach $param_name ($query->param()) {
    $param_value = $query->param($param_name);
    $content .=
        "\t\t\t<tr>\n".
        "\t\t\t\t<td>$param_name</td>\n".
        "\t\t\t\t<td>$param_value</td>\n".
        "\t\t\t</tr>\n";
  }
  $content .= 
      "\t\t</table>\n".
      "\t</body>\n".
      "</html>";
  $content;
}

#
# Submit a block to the blocks processor
#
sub submit_block {
  my($blocks_url, $start, $end) = @_;
  my($content, $blocks, $i);

  # get the BLOCK(S)
  $blocks = "";
  for ($i = $start; $i <= $end; $i++) {
    $blocks .= param('BLOCKS'.$i);
  }

  my $ua = LWP::UserAgent->new();
  my $req = POST $blocks_url, [ sequences => $blocks ];
  my $request = $ua->request($req);
  my $response = $request->content;
  # put in the absolute url's : this is FRAGILE!
  $response =~ s#HREF=\"/#HREF=\"http://blocks.fhcrc.org/#g;
  $response =~ s#ACTION=\"/blocks-bin#ACTION=\"http://blocks.fhcrc.org/blocks-bin#g;
  $content = 
      "<html>\n".
      "\t<head>\n".
      "\t\t<title>MEME Suite - Submit block</title>\n".
      "\t</head>\n".
      "\t<body style=\"background-image:url('../images/bkg.jpg');\">\n".
      "$response\n".
      "\t</body>\n".
      "</html>";

  $content;
} # submit_block

#
# Compare a PSPM to the TOMTOM motif database search tool
#
sub submit_pspm_to_tomtom {
  my($tomtom_url, $number) = @_;
  my(@fields, $nsites, $i, $n, $w, $pspm, $row, $col);

  # get the motif PSPM 
  $_ = param('pspm'.$number);
  @fields = split;

  # ignore 1st 10 entries
  $nsites = $fields[7];# save number of sites
  for ($i = 0; $i <= $#fields - 10; $i++) { 
    $fields[$i] = $fields[$i+10] * $nsites; 
  }
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
  my $ua = LWP::UserAgent->new();
  my $req = POST "$tomtom_url",
    Content_Type => 'multipart/form-data',
    Content => [ 'query' => $pspm, 'query_name' => "Motif_$number" ];
  my $request = $ua->request($req);

  my $content = $request->content;

  # return the page
  $content;
} # submit_pspm_to_tomtom

sub generate_logo {
  my ($number) = @_;
  my ($bin) = "@MEME_DIR@/bin";

  binmode STDOUT;

  #print "Content-type: text/plain", "\n\n";
  #print "Hello World";
  #exit(0);
  my ($content) = "";


  my ($motifs) = get_motif($number, $number);

  #print "Content-type: text/plain", "\n\n";
  #print $motifs;
  #exit(0);

  # write the motifs to a temporary file that will be deleted when perl exits
  my ($fh, $tmpname);
  ($fh, $tmpname) = tempfile(UNLINK => 1);
  print($fh $motifs);
  close($fh);

  #print "Content-type: text/plain", "\n\n";
  #print "cat $tmpname", "\n";
  #print `cat $tmpname`;
  #exit(0);

  # get parameters
  my ($ssc, $png, $rc, $width, $height, $outname);
  $ssc = ($query->param('logossc_' . $number) eq "true");
  $png = ($query->param('logoformat_' . $number) eq "png");
  $rc  = ($query->param('logorc_' . $number) eq "true");
  $width = ($query->param('logowidth_'.$number));
  $height = ($query->param('logoheight_'.$number));

  if ($width =~ m/^\s*(\d+(.\d+)?)\s*$/) {
    $width = "-w $1";
  } else {
    $width = "";
  }

  if ($height =~ m/^\s*(\d+(.\d+)?)\s*$/) {
    $height = "-h $1";
  } else {
    $height = "";
  }

  # generate the name that should be specified for downloading
  $outname = "logo";
  if ($rc) {
    $outname .= "_rc";
  }
  if ($ssc) {
    $outname .= "_ssc";
  }
  $outname .= $number;
  $outname .= ".";
  if ($png) {
    $outname .= "png";
  } else {
    $outname .= "eps";
  }

  #print "Content-type: text/plain", "\n\n";
  #print "ceqlogo -k AA -i $tmpname", "\n";
  #print `@MEME_DIR@/bin/ceqlogo -k AA -i $tmpname`;
  #print "cat $tmpname", "\n";
  #print `cat $tmpname`;
  #exit(0);

  # find out the format
  my ($alphabet_size) = length($query->param('alphabet'));
  my ($kind);
  if ($alphabet_size < 20) {
    $kind = "NA";
  } else {
    $kind = "AA";
  }

  # setup ssc toggle
  my ($ssc_toggle) = "";
  if ($ssc) {
    $ssc_toggle = "-E -S";
  }


  # setup reverse complement flag as required
  # note that only nucleic acid can be reverse complemented
  my ($rc_flag) = "";
  if ($rc and $alphabet_size < 20) {
    $rc_flag = "-r";
  }

  my $fineprint = "MEME ";
  if ($ssc) {
    $fineprint = $fineprint."(with SSC)";
  } else {
    $fineprint = $fineprint."(no SSC)";
  }
  my (undef,$min,$hour,$day,$month,$year,undef,undef,undef) = localtime(time);
  $year += 1900; # year contains years since 1900
  $month += 1; # month counted from zero.
  $fineprint .= sprintf("%d.%d.%d %02d:%02d", $day, $month, $year, $hour, $min);

  # generate the image content
  if ($png) {
    $content = `$bin/ceqlogo -Y -N $ssc_toggle -k $kind $rc_flag -d \"$fineprint\" -t \"\" -x \"\" $width $height -i $tmpname | $bin/convert - png:-`;
  } else {
    $content = `$bin/ceqlogo -Y -N $ssc_toggle -k $kind $rc_flag -d \"$fineprint\" -t \"\" -x \"\" $width $height -i $tmpname`;
  }
  #print "Content-type: text/plain\n\n";
  #print $content;
  #exit(0);

  # this header tells the browser to download the file
  print "Content-Disposition: attachment; filename=$outname", "\n"; 
  if ($png) {
    print "Content-type: image/png", "\n";
  } else {
    print "Content-type: application/postscript", "\n";
  }
  print "\n";
  print $content;
  exit(0);
}
