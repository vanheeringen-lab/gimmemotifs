#!@WHICHPERL@
# FILE: mhmm_process_request.cgi
# AUTHOR: Haoyuan Zhu
# CREATE DATE: 11-3-2001
# PROJECT: MHMM
# DESCRIPTION: Sumbit to MHMM server from MEME output page

use CGI qw/:standard/;
use HTTP::Request::Common qw(POST);
use LWP::UserAgent;
use lib qw(@PERLLIBDIR@);
use mhmm_globals;
use Globals;

# get the directory for executables relative to working directory
$scratch = "$MEME_LOGS";
$bin = "$MEME_BIN";
$cgidir = "$MEME_WEB";
# email address of maintainer; this is set by the install script
$maint = "$SITE_MAINTAINER";

$ori_seqfile = "ori_seqfile.$$";
$job_log = "oldmeme.log.$$";

$cgi = new CGI;

print header;
print "<BLOCKQUOTE>";

chdir ("$scratch") || die "Cannot cd to $scratch";

$process_dna = 0;
$process_protein = 0;
$meme_head = "MEME version 3.0";
$bg_head = "Background letter frequencies (from dataset with add-one prior applied):";
$motif_summary_begin = "<INPUT TYPE = HIDDEN NAME = motif-summary VALUE = \"";
$motif_summary_end = "\">";
$nmotifs = param('nmotifs');
$name = param('name');
$root = "..";

# Get the info. MHMM needs
$alphabet = param('alphabet');
if (length($alphabet) >= 20){
	$process_protein = 1;
}
else{
	$process_dna = 1;
}

$bgfreq = param('bgfreq');
$motif_summary = param('motif-summary');

for ($i=1;$i<=$nmotifs;$i++){
    $pspm[$i] = param('pspm'.$i);
}

#print $motif_summary;

#open (OUT,"> submit_meme.$process_record") 
#open (OUT, "</misc/www/projects/MEME/METAMEME/beta/website/temp/submit_meme.$record")
#      || die "Cannot open $scratch/submit_meme.$process_record";
# Work with old MEME output version

if (defined $bgfreq && $bgfreq eq ""){
    #do 'oldmeme_process_request.cgi';
    #exit(0);
    #print "<HTML>";
    #if (!$cgi->param){
    #$old_motif_summary = $motif_summary;
    #print $motif_summary;
    $original = param("originalseq");
    if (!$original){
	print_intro();
    }
    else{
	verify($cgi);
	#print $meme_head;
    }
    #else{
	#process_mhmm();
    #}
}

if (defined $bgfreq && $bgfreq ne ""){
    process_mhmm();
}

sub process_mhmm{
    # Open the METAMEME submission form
    open(METAFORM,"< $root/submit-verify.html")
	|| die "Cannot open METAMEME submission form submit-verify.html\n";
    
    #read the submission form and remove motif file,  
    #put in inline name and motif

    while(defined($line=<METAFORM>)){
	
	if ($line =~ /name=\"seqtype\"/){
	    if ($process_protein == 1){ 
		print "<b>protein&nbsp</b>";
		print "<input type = hidden name = seqtype value = \"protein\">";
	    }
	    elsif ($process_dna == 1){
		print "<b>DNA&nbsp</b>";
		print "<input type = hidden name = seqtype value = \"DNA\">";
	    }
	    else{
		print $line;
	    }
	    
	}
	
	elsif($line =~ /MEME motifs/){
	    print "<b>Your motif file</b>: MEME results on <b>$name</b>.\n";
	}
	
	elsif ($line=~ /name=\"meme\"/){
	    print "<INPUT TYPE = HIDDEN NAME = inline_motifs VALUE = \n\"";
	    print_meme();
	    print "\">\n";
	    print_motif_summary();
	    print_meme_name();
	}
	
	else{
	    print $line;
	}
    }
    close (METAFORM);
}

sub print_meme{
    print "$meme_head\n\n";
    print "ALPHABET= ";
    print "$alphabet\n\n";
    print "$bg_head$bgfreq\n\n";
    for ($i=1;$i<=$#pspm;$i++){
	print "MOTIF $i\n";
	print $pspm[$i];
	print "\n\n";
    }
    print "Combined\n\n";
}

sub print_motif_summary{
    print $motif_summary_begin;
    #chomp $motif_summary;
    print $motif_summary;
    print "$motif_summary_end\n";
}

sub print_meme_name{
    print "<INPUT TYPE = Hidden name = memename value = \"";
    print "Your MEME result on $name";
    print "\">\n";
}

sub print_intro{
    print $cgi->start_html(-"title"=>"Input your sequence file",
			   -"style"=>{'src'=>"../metameme.css"},
			   -target=>'_self',
			   -bgcolor=>"White",
			   -"script"=>{-language=>'JAVASCRIPT',
				     -src=>'../validate.js'}),
    "<IMG ALIGN='LEFT' SRC='../images/meta-meme.gif' Height='81' Width='444'>",
    "<BR><BR><BR><BR><BR>",
    "<H3>Data submission Form</H3>",
    "It appears that your motif file is from an old version of ",
    $cgi->a({href=>"http://www.sdsc.edu/MEME"},"MEME"),
    ". It lacks some information that Meta-MEME requires. ",
    "Please use this form to upload the original sequences that you input to ",
    $cgi->a({href=>"http://www.sdsc.edu/MEME"},"MEME"),
    ".",
    "<P>",
    "<HR>",
    $cgi->start_multipart_form(-name => 'form1',
			       -onSubmit=>"return validateForm()"),
    "Please input the original sequence (<B>$name</B>) you gave to ",
    $cgi->a({href=>"http://www.sdsc.edu/MEME"},"MEME"),
    "<BR>",
   
    "<input type=\"file\" name=\"originalseq\" size=\"40\" maxlength=\"80\">",
    "<BR><BR>",
    "<input type=\"hidden\" name=\"nmotifs\" value=\"$nmotifs\">",
    "<input type=\"hidden\" name=\"name\" value=\"$name\">",
    "<input type=\"hidden\" name=\"alphabet\" value=\"$alphabet\">",
    "<input type=\"hidden\" name=\"bgfreq\" value=\"$bgfreq\">",
    "<input type=\"hidden\" name=\"motif-summary\" value=\"$motif_summary\">",
    $cgi->submit(-label=>'Submit'),
    "&nbsp",
    $cgi->reset(-name=>'Clear Input'),
    ;
    
    for ($i=1;$i<=$nmotifs;$i++){
	print "<input type=\"hidden\" name=\"pspm$i\" value=\"$pspm[$i]\">";
	
    }

    print $cgi->end_form;
    print "<HR>";

    print $cgi->a({href=>"../metameme-intro.html"}, "Return"),
    " to the Meta-MEME home page.<BR><BR>",
    "<I>Please send comments and questions to ",
    $cgi->a({href=>"mailto:$webmaster"},
	    "$webmaster"),
    "</I></BLOCKQUOTE>",

    $cgi->end_html;
}    

sub verify{
    local ($form) = $_[0];
    
    if ($form->param('originalseq') eq "") {
	&return_error($global_log, $job_log, "No sequence file given",
		      "You must provide a input sequence");
    }

    chmod(0444,$form->param('originalseq'));
    foreach my $infile ($form->param('originalseq')) {
	if (-z $infile) {
	    &return_error($global_log, $job_log, "Empty file",
                    "Meta-MEME could not read the following file:<BR>
                     &nbsp; &nbsp; &nbsp; $infile<BR>
                     Make sure the name is correct and that you have read
                     access to the file.");
	}
    }

    my $oriseq_file = $form->tmpFileName($form->param('originalseq'));
    my $oriseq = `cat $oriseq_file`;
    if (open($ori_seqfile, ">$ori_seqfile")) {
	$sub = "";
	$oriseq =~ s/\r/$sub/g;
	print($ori_seqfile $oriseq);
	close($ori_seqfile);
    }
    chmod(0444, $ori_seqfile);
    # print $ori_seqfile;
    # Convert the data to FASTA format.
    my $ori_readseq_file = "ori_readseq.$$";
    my $ori_readseq_error = "ori_readseq.error.$$";
    $status = system("$bindir/readseq -a -f=8 $ori_seqfile \\
                    1> $ori_readseq_file 2> $ori_readseq_error");
    $error = `cat $ori_readseq_error`;
    unlink($ori_readseq_file, $ori_readseq_error);
    if ($status) {
	&return_error($global_log, $job_log, "READSEQ error", 
		      "An error occurred when the READSEQ program attempted
                       to convert your dataset to FASTA format.<BR>
                       READSEQ returned: $error");
    }
    #print $motif_summary;    
    #print $record{"motif_summary"};
    #print $new_name;
    #print $form->param('name');
    #@temp = split('\n',$motif_summary);
    #print $temp[0],$temp[5];

    # Run 'name_length_bg.fasta' program to get information on converted data
    my $getinfo_file = "getinfo.$$";
    my $getinfo_error = "getinfo.error.$$";
    my $status = system("$bindir/name_length_bg.fasta $ori_seqfile \\
                         1> $getinfo_file 2> $getinfo_error");
    my $error = `cat $getinfo_error`;
    my $getinfo = `cat $getinfo_file`;
    if ($error || $status || ($getinfo eq "")){
	&return_error($global_log, $job_log, "Sequence error",
		      "After converting to FASTA format using the READSEQ
                       program, the following errors in your dataset were
                       detected: <BR>$error.\n");
    }
    @seq_info = split('\n',$getinfo);
    @ori_name = split(' ',$seq_info[0]);
    @ori_length = split(' ',$seq_info[1]);
    $bgfreq = "\n".$seq_info[2];

    foreach (@ori_length){
	if ($_ == 0){
	    &return_error($global_log,$job_log,"Sequence error",
			  "Your sequence file appears to contain one or more zero-length
                           sequences.<BR> Please check to be sure that your data is 
                           <A HREF=../help/help_format.html> formatted</A> properly.");
	}
    }

    @motif_summary = split('\n',$motif_summary);
    
    if ($#motif_summary != $#ori_name+1){
	&return_error($global_log, $job_log, "Sequence error",
		      "It appears the sequence file you gave is not the same
                       as as you gave to MEME.\n");
    }
    
    $motif_summary = "\n";
    for ($i=0;$i<=$#ori_name;$i++){
	@temp = split(' ',$motif_summary[$i+1]);
       
	#if ($ori_name[$i] ne $temp[0]){
	if (!($temp[0] =~ $ori_name[$i])){
	     &return_error($global_log, $job_log, "Sequence error",
			  "It seems the sequence file you gave is not the same
                          as as you gave to MEME.\n");
	}
	else{
	    for ($j=0;$j<3;$j++){
		$motif_summary = $motif_summary.$temp[$j]." ";
	    }
	    if ($#temp % 3 != 0){
		$motif_summary = $motif_summary.$ori_length[$i]." ";
	    }
	    for ($j=3;$j<=$#temp;$j++){
		$motif_summary = $motif_summary.$temp[$j];
		if ($j < $#temp){
		    $motif_summary = $motif_summary." ";
		}
		else{
		    $motif_summary = $motif_summary."\n";
		}
	    }
	}
    }
    
}


###############################################################################
# Return error
###############################################################################
sub return_error {
  my ($global_log, $job_log, $keyword, $message) = @_;
  my ($field, $date_and_time);

  # Print an error log for this job.
  open($job_log, ">$job_log");
  print($job_log "keyword=$keyword\n");
  print($job_log "message=$message\n");
  foreach $field ('originalseq') {
    my $value = $form->param($field);
    print($job_log "$field=$value\n");
  }
  close($job_log);
  chmod(0666, $job_log);

  # Add an entry to the global log file.
  chop($date_and_time = `date`);
  open($global_log, ">>$global_log");
  print($global_log "$date_and_time \t $$ \t <1> ");
  foreach $field ('originalseq') {
    my $value = $form->param($field);
    print($global_log "$field=$value; ");
  }

  print($global_log "\n\n");
  close($global_log);
  chmod(0666, $global_log);

  # 'form' has dynamic scope from the caller
  print $form->start_html(-"title"=>"Meta-MEME error - $keyword",
                          -target=>'_self');
  print $form->h1($keyword),
  $form->hr($message),
  $form->p("If you have questions, please contact ",
	   "<A HREF=\"mailto:$webmaster\">$webmaster</A>.<BR>");
  print "</BLOCKQUOTE>\n";
  print $form->end_html;
  exit(1);
}













