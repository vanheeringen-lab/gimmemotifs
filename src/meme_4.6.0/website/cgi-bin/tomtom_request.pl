#!@WHICHPERL@ -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use List::Util qw(sum);
use File::Temp qw(tempfile);

my $bin = "@MEME_DIR@/bin";
my $q = CGI->new;
dispatch();

#
# dispatch the handling to the approprate
# subroutine dependent on the mode parameter
#
sub dispatch {
  my $mode = $q->param('mode');
  if (not defined $mode) {
    unknown_request();
  } else {
    if ($mode =~ m/^logo$/i) {
      tomtom_logo();
    } else {
      unknown_request();
    }
  }
}

sub unknown_request {
  print $q->header,
        $q->start_html('Unrecognised request'),
        $q->h1('Unrecognised request'),
        $q->p('The value specified for the parameter mode tells this script what to do. '.
            'The value that was specified was unrecognised and so this script doesn\'t know what to do.');
  print $q->h2('Parameters');
  my @rows = $q->Tr($q->th('name'), $q->th('value'));
  my $param;
  foreach $param ($q->param) {
    push @rows, $q->Tr($q->td($param), $q->td($q->pre($q->param($param))));
  }
  print $q->table(@rows);
  print $q->end_html;
}

sub tomtom_logo {
  my $version = $q->param('version');
  die("Input version must be up to three dot separated numbers.") unless ($version =~ m/^\d+(\.\d+(\.\d+)?)?$/);
  my $bgline = parse_background($q->param('background')); 
  my $target_id = $q->param('target_id');
  die("Input target_id must not contain spaces.") unless ($target_id =~ m/^\S+$/);
  my $target_length = $q->param('target_length');
  die("Input target_length must be a positive integer.") unless ($target_length =~ m/^\d+$/);
  my $target_pspm = $q->param('target_pspm');
  die("Input target_pspm must be defined.") unless defined $target_pspm;
  my $target_rc_str = $q->param('target_rc');
  die("Input target_rc must be either 1 or 0.") unless ($target_rc_str =~ m/^[01]$/);
  my $target_rc = int($target_rc_str);
  
  my $query_id = $q->param('query_id');
  die("Input query_id must not contain spaces.") unless ($query_id =~ m/^\S+$/);
  my $query_length = $q->param('query_length');
  die("Input query_length must be a positive integer.") unless ($query_length =~ m/^\d+$/);
  my $query_pspm = $q->param('query_pspm');
  die("Input query_pspm must be defined.") unless defined $query_pspm;
  my $query_offset_str = $q->param('query_offset');
  die("Input query_offset must be a number.") unless ($query_offset_str =~ m/^[-]?\d+$/);
  my $query_offset = int($query_offset_str);

  my $error_bars = $q->param('error_bars');
  die("Input error_bars but be 0 or 1.") unless ($error_bars =~ m/^[01]$/);
  my $small_sample_correction = $q->param('small_sample_correction');
  die("Input small_sample_correction must be 0 or 1.") unless ($small_sample_correction =~ m/^[01]$/);
  my $flip = $q->param('flip');
  die("Input flip must be 0 or 1.") unless ($flip =~ m/^[01]$/);


  my $img_w = "";
  my $image_width = $q->param('image_width');
  if (defined $image_width && $image_width =~ m/^\d+$/) {
    $img_w = "-w $image_width ";
  }
  my $img_h = "";
  my $image_height = $q->param('image_height');
  if (defined $image_height && $image_height =~ m/^\d+$/) {
    $img_h = "-h $image_height ";
  }
  my $image_type = $q->param('image_type');
  die("Input image_type must be either png or eps.") unless ($image_type =~ m/^(png|eps)$/);

  my $q_file = meme_nucleotide_file($version, $bgline, $query_id, $query_pspm);
  my $t_file = meme_nucleotide_file($version, $bgline, $target_id, $target_pspm);
  
  my $eb = ($error_bars ? "-E " : "");
  my $ssc = ($small_sample_correction ? "-S " : "");
  my $t_rc = (($flip ? !($target_rc) : $target_rc) ? "-r " : "");
  my $q_rc = ($flip ? "-r " : "");
  my $q_off = ($flip ? $target_length - ($query_length + $query_offset) : $query_offset);
  my $t_shift = ($q_off < 0 ? "-s ".(-$q_off)." " : "");
  my $q_shift = ($q_off > 0 ? "-s ".($q_off)." " : "");

  my $fine_text = "Tomtom " . date_time_string();
  
  my $ceqlogo_cmd = "$bin/ceqlogo -Y -N $eb $ssc -k NA $img_w $img_h -d \"$fine_text\" -t \"$target_id\" -x \"$query_id\" $t_rc -i $t_file $t_shift $q_rc -i $q_file $q_shift";

  my $image;
  my $mime_type;
  if ($image_type =~ m/^eps$/) {
    $image = `$ceqlogo_cmd`;
    $mime_type = "application/postscript";
  } elsif ($image_type =~ m/^png$/) {
    $image = `$ceqlogo_cmd | $bin/convert - png:-`;
    $mime_type = "image/png";
  }

  my $outname = "logo_" . $query_id . "_" . $target_id . "." . $image_type;
  binmode STDOUT;
  
  print $q->header(-type=>$mime_type, -Content_Disposition=>"attachment; filename=$outname");
  print $image;
}

sub parse_background {
  my $bgtext = shift;
  die("parse_background: expected background probabilities string.\n") unless defined $bgtext;
  die("parse_background: expected 4 probabilities in a space delimited string.") unless ($bgtext =~ m/^0\.\d+\s+0\.\d+\s+0\.\d+\s+0\.\d+$/);
  my @bgprobs = split(/\s+/, $bgtext);
  die("parse_background: error expected the split to be of length 4 but it wasn't.") unless (@bgprobs == 4);
  my $probsum = sum(@bgprobs);
  die("parse_background: expected probabilities to sum to one (or close).") unless abs($probsum - 1) < 0.1;
  return "A $bgprobs[0] C $bgprobs[1] G $bgprobs[2] T $bgprobs[3]";
}

sub meme_nucleotide_file {
  my $motif = meme_nucleotide_motif(@_);
  # write the motifs to a temporary file that will be deleted when perl exits
  my ($fh, $tmpname);
  ($fh, $tmpname) = tempfile(UNLINK => 1);
  print $fh $motif;
  close $fh;
  return $tmpname;
}

sub meme_nucleotide_motif {
  my ($version, $bgline, $id, $pspm) = @_;

  return "MEME version $version\n" .
         "ALPHABET= ACGT\n" .
         "strands: + -\n" .
         "Background letter frequencies (from an unknown source):\n" .
         $bgline . "\n" .
         "MOTIF $id\n" .
         $pspm;
}

sub date_time_string {
  my (undef,$min,$hour,$day,$month,$year,undef,undef,undef) = localtime(time);
  $year += 1900; # year contains years since 1900
  $month += 1; # month counted from zero.
  return sprintf("%d.%d.%d %02d:%02d", $day, $month, $year, $hour, $min);
}
