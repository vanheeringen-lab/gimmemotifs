#!/bin/perl -w

=head1 NAME

treg_comparator.pl

=head1 SYNOPSIS

treg_comparator.pl -q <query matrix file> -f <query matix format (default raw)> -d <matrix set file> -g <vertebrate> -r <Jaspar> -ori <orientation> -c <cutoff> -out <output file>

=head1 DESCRIPTION

This script is a stand alone program version of the T-Reg Comparator web application,
which is available under http://treg.molgen.mpg.de. 

=head1 EXAMPLE

treg_comparator.pl -q example_input.meme -f MEME -d matrix_set_Jaspar -g vertebrate -r JASPAR -ori both -c 0.9 -out result

treg_comparator.pl -q pwms.FASTAlike  -f raw -d matrix_set_Jaspar -g vertebrate -r JASPAR  -ori both -c 0.9 -out result

=head1 OPTIONS

=over 8

=item -h	

help

=item -q

Query matrix file in raw or MEME format	
If you want to compare matrices of the matrix set (option -d) among each other, 
you can just ommit the -q option.

=item -f	

Type or format of the matrices in the query matrix file
'MEME' or 'raw'

=item -d	

Matrix file against which the query matrices are planned to be compared.
This file has to have the following format:
	M00001#E47#vertebrate#TRANSFAC#1#0.2 0.3548828125 ...
Each weight matrix is coded in one line, first comes the name terminated by "#" 
then comes the name of the factor, the species group, the set name like JASPAR 
and then all matrix entries (float) separated by a whitespace. The matrix entries 
start with the frequency of "A", followed by that of "C", "G", and "T" in the 
first position, then the same fo the second position and so on. 

=item -g	

Restriction of the matrix set to a species group, default is "all", meaning no restriction.
Other than "all", each string that matches values in the third position of the matrix file 
can be entered. (like "vertebrate" in the example above)

=item -r	

Restirction of the matrix set to specialized subsets, default is "all" meaning 
no restriction. "input" means that the query matrices are compared among each 
other and no matrix file (option -d) is needed. Other than "all" and "input", 
you can enter each string that matches values in the third position of the 
matrix file. (like "TRANSFAC" in the example above)

=item -ori	

Orientation takes the values "given", "complement", "both". 
default is "both"

=item -c	

Cutoff for comparison
Default is 0.8

=item -out	

Output file

=back

=head1 AUTHOR

Stefan Roepcke, roepcke@molgen.mpg.de

=cut

use strict;
use POSIX;
use Storable qw(dclone);
use Getopt::Long;
use CGI qw/:standard :html3/;
use lib qw(./treg_comparator);
use treg_utilities;

#####
#use lib qw(./treg_comparator);
#use treg_utilities;

use FindBin ();
use lib "$FindBin::RealBin";
use treg_utilities;
#####

my $DEBUG = 0;

my $help=0;
my ($query_file_name, $query_format, $reg_pmat_file_name, $species_group, $source_set, $ori, $cutoff, $output_file_name) = ("","","","all","all","both", 0.8,"result"); 
&GetOptions("h!" => \$help,
	"q=s" => \$query_file_name,
	"f=s" => \$query_format,
	"d=s" => \$reg_pmat_file_name,
	"g=s" => \$species_group,
	"r=s" => \$source_set,
	"ori=s" => \$ori,
	"c=f"  => \$cutoff,
	"out=s" => \$output_file_name );
exec('perldoc',$0) if $help;
my $format = "";

print STDOUT "Call parameters: $query_file_name, $query_format, $reg_pmat_file_name, $species_group, $source_set, 
	$cutoff, $ori, $output_file_name\n"  if $DEBUG;

my (%pwm_name, %class, %species, %source, %pmat_orientation, %consensus, %RC_consensus, %reg_pmat_TREG, %RC_reg_pmat_TREG, %pmat_length);
my (%query_set_ids, %db_set_ids);

### 1. Parse query file ###
if( $query_file_name ne "" )  {
	open( QU, "$query_file_name" )  or  die "Can't open $query_file_name : $!\n";
	my @rows = <QU>;
	if( uc($query_format) eq "MEME" )  {
		my $first_line = $rows[0] or die "Empty PWM input file!\n";
		if( $first_line =~ m/\*\*\*/ or uc($query_format) eq "MEME" )  {   ### Format of the file looks like MEME
			process_meme_PWMs( \@rows, \%query_set_ids, \%reg_pmat_TREG, \%RC_reg_pmat_TREG, \%source, \%consensus, \%RC_consensus, \%pmat_length ); 
			# fills %query_set_ids, generates reg_pmat, fills %reg_pmat_TREG, %RC_reg_pmat_TREG 
		}
	}  elsif( $query_format eq "raw" )  {
		process_raw_PWMs( \@rows, \%query_set_ids, \%reg_pmat_TREG, \%RC_reg_pmat_TREG, \%source, \%consensus, \%RC_consensus, \%pmat_length );
	}  else  { print STDERR "Can't read your format: $query_format! Please enter 'raw' or 'MEME'.\n"; }

	close(QU);
	foreach  my $id (keys %query_set_ids)  { $pwm_name{$id} = $id; }
}  ### else: matrices of the matrix file are compared against each other

###########################################################################################################################################
### 2. Load matrix file ###
#M00001#E47#vertebrate#TRANSFAC#0.2 0.3548828125 ...
# PWMs herein are required to be regularized
if( $source_set eq "input" )  { 
	undef(%db_set_ids); %db_set_ids = %query_set_ids; 
}  else  { 
	open( M, $reg_pmat_file_name )  or  print STDOUT "Can't open $reg_pmat_file_name !\n";
	while( my $line = <M> )  {
		chomp $line;
		my ($id, $name, $sp, $src, $reg_pmat_str) = split( "#", $line );  ##, $o
		#my $o;
		if( defined $id )  { 
			if( defined $name )  { $pwm_name{$id} = $name; }  else  { $pwm_name{$id} = $id; }
			if( defined $sp and $sp =~ m/[a-zA-Z]/ )   { 
				$species{$id} = $sp;  ##['vertebrate','insect','yeast','plant', 'all']
			}  else  { $species{$id} = "not specified"; }
			if( defined $src and $src =~ m/[a-zA-Z]/ )  {  ### 'JASPAR', 'all', ...
				$source{$id} = $src; 
				if( uc($source_set) =~ uc($src) or $source_set eq 'all' or $source_set eq 'input' )  { 
					$db_set_ids{$id} = 'TRUE';
				} 
			}  
			else  { $source{$id} = "not specified"; $db_set_ids{$id} = 'TRUE'; }
			#if( defined $o )  { $pmat_orientation{$id} = $o; }  else  { $pmat_orientation{$id} = 0; }
			if( defined $reg_pmat_str )  {
				$reg_pmat_str =~ s/^\s+//; 
				my @reg_pmat = split( /\s+/, $reg_pmat_str );
				#print STDOUT  $id, "\t", join( " ", @reg_pmat ), "\n"  if $DEBUG;
				if( scalar(@reg_pmat) > 0 )  { 
					$reg_pmat_TREG{$id} = \@reg_pmat; 
					$consensus{$id} = get_iupac_from_pmat(\@reg_pmat);
					$pmat_length{$id} = scalar(@reg_pmat) / 4;
					my @tmp = revcom_matrix( \@reg_pmat );
					$RC_reg_pmat_TREG{$id} = \@tmp;
					$RC_consensus{$id} = get_iupac_from_pmat(\@reg_pmat);
				} 
			}  else  { print STDOUT "No values in reg_pmat for $id !\n"; }
		}
	}
	close(M);
	if( $query_file_name eq "" )  { %query_set_ids = %db_set_ids; }
}

print STDOUT "Query Ids: ", keys(%query_set_ids), "\n"  if $DEBUG;
print STDOUT "DB Idssss: ", keys(%db_set_ids), "\n"  if $DEBUG;



###########################################################################################################################################
### 3. Loop over all query PWMs ... loop over all DB subjects
my $result_file = "results.html";
open( RES, ">$output_file_name.html") or  print STDERR "Can't open $output_file_name.html !\n";
open( TXT, ">$output_file_name.txt") or  print STDERR "Can't open $output_file_name.txt !\n";
print TXT  "Query_ID", "\tQuery_Consensus", "\tSubject_Name", "\tSource_DB", "\tSubject_ID", "\tLength", 
	"\tOrientation", "\tOffset", "\tDivergence", "\tOverlap", "\tSubject_Consensus", "\n";
open( CS, ">$output_file_name.csv") or  print STDERR "Can't open $output_file_name.csv !\n";
print CS  "Query_ID", ";Query_Consensus", ";Subject_Name", ";Source_DB", ";Subject_ID", ";Length", ";Orientation", 
	";Offset", ";Divergence", ";Overlap", ";Subject_Consensus", "\n";

foreach  my $query_id (keys %query_set_ids)  {
	#print STDOUT ".";
	print STDOUT "\nQuery Matrix: ", $query_id, ": ", $consensus{$query_id}, ", $source_set: ", join( "#", time), "\n"  if $DEBUG;
	print RES "<h3>Query Matrix ", $query_id, "</h3>";
	print RES "<p>Consensus: ", $consensus{$query_id}, ",  Length: ", $pmat_length{$query_id}, "\n";
	my (%tablecontent_hash, %txt_tablecontent_hash, %cstxt_tablecontent_hash, %divergence_hash);
	my $tablecontent = "";
	my $txt_tablecontent = "";
	my $cstxt_tablecontent = "";
	foreach  my $subject (keys %db_set_ids)  {
		if( uc($species_group) =~ uc($species{$subject})  or  $species_group eq "all" )  {
			print STDOUT "Subject: ", $subject, "\n"  if $DEBUG;
			if( ($ori =~ "given" or $ori =~ "both") and $query_id ne $subject )  {
				my $orientation = "as given";
				my ($divergence, $offset, $overlap, $info_cont1, $info_cont2) = compare_partial( $reg_pmat_TREG{$subject}, $reg_pmat_TREG{$query_id} );
				print STDOUT "Results: ", $divergence, "\t", $offset, "\t", $overlap, "\t", $info_cont1, "\t", $info_cont2, "\n"  if $DEBUG;
				#my ($divergence, $offset, $overlap) = compare_partial( $reg_pmat_TREG{$subject}, $reg_pmat_TREG{$query_id} );
				#print STDOUT "Results: ", $divergence, "\t", $offset, "\t", $overlap, "\n"  if $DEBUG;
				if( $divergence < $cutoff )  {
					my $rowcontent = "";
					my $txt_rowcontent = $query_id . "\t" . $consensus{$query_id} . "\t";
					my $cstxt_rowcontent = $query_id . ";" . $consensus{$query_id} . ";";
					#my $source = "TREG";
					## Jaspar link: javascript:Start('http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?ID=MA0101&rm=present')
					my $ref_str;
					my $counter = 303;
					my $x = $subject;
					$x =~ s/^M(0){1,4}//;
					if( $source{$subject} =~ "TRANSFAC" )  { $ref_str = "<a href=\"http://www.biobase.de/cgi-bin/biobase/transfac/8.4/bin/getTFProf.cgi\?$subject\">$subject</a>"; }
					elsif( $source{$subject} =~ "PUBLIC" and $counter <= 302 )  { 
						$ref_str = "<a href=\"http://www.cbil.upenn.edu/cgi-bin/tess/tess?request=MTX-DBRTRV-Accno&key=$subject\">$subject</a>"; 
					}  elsif( $source{$subject} =~ "PUBLIC" and $counter > 302 )  { 
						$ref_str = "<a href=http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=\"$subject\">$subject</a>"; 
						#http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=M00001
					}  elsif( $source{$subject} =~ "JASPAR" )  { $ref_str = "<a href=\"http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?ID=$subject&rm=present\">$subject</a>"; 
					}  else  { 
						if( defined $source{$subject} )   { print STDOUT  "source: ", $source{$subject}, "\n"  if $DEBUG; }
						else  { print STDOUT  "NO source: ", $subject, "\n"  if $DEBUG; }
						if( $source{$subject} ne "Query" )  { $source{$subject} = "--"; }
						$ref_str = "-";
					}
					my $d = sprintf("%.3f", $divergence);
					$rowcontent .= td($pwm_name{$subject}); $txt_rowcontent .= $pwm_name{$subject} . "\t"; $cstxt_rowcontent .= $pwm_name{$subject} . ";";
					#$rowcontent .= td($class{$subject}); $txt_rowcontent .= $class{$subject} . "\t"; $cstxt_rowcontent .= $class{$subject} . ";";
					$rowcontent .= td($source{$subject}); $txt_rowcontent .= $source{$subject} . "\t"; $cstxt_rowcontent .= $source{$subject} . ";";
					$rowcontent .= td($ref_str); 		$txt_rowcontent .= $subject . "\t"; $cstxt_rowcontent .= $subject . ";";
					$rowcontent .= td($pmat_length{$subject}); $txt_rowcontent .= $pmat_length{$subject} . "\t"; $cstxt_rowcontent .= $pmat_length{$subject} . ";";
					$rowcontent .= td($orientation); 	$txt_rowcontent .= $orientation . "\t"; $cstxt_rowcontent .= $orientation . ";";
					$rowcontent .= td($offset); 		$txt_rowcontent .= $offset . "\t"; $cstxt_rowcontent .= $offset . ";";
					$rowcontent .= td($d); 			$txt_rowcontent .= $d . "\t"; $cstxt_rowcontent .= $d . ";";
					$rowcontent .= td($overlap); 		$txt_rowcontent .= $overlap . "\t"; $cstxt_rowcontent .= $overlap . ";";
					$rowcontent .= td($consensus{$subject}); $txt_rowcontent .= $consensus{$subject} . "\t"; $cstxt_rowcontent .= $consensus{$subject} . ";";
					$tablecontent .= TR($rowcontent); $txt_tablecontent .= $txt_rowcontent . "\n";  $cstxt_tablecontent .= $cstxt_rowcontent . "\n";
					
					my $k = $pwm_name{$subject} . "_" . $orientation;
					$tablecontent_hash{$k} = TR($rowcontent);
					$txt_tablecontent_hash{$k} = $txt_rowcontent . "\n";
					$cstxt_tablecontent_hash{$k} = $cstxt_rowcontent . "\n";
					$divergence_hash{$k} = $d;
				}
			}
			if( $ori =~ "complement" or $ori =~ "both" )  {
				my $orientation = "reverse-complement";
				my ($divergence, $offset, $overlap, $info_cont1, $info_cont2) = compare_partial( $RC_reg_pmat_TREG{$subject}, $reg_pmat_TREG{$query_id} );
				print STDOUT "Results RC: ", $divergence, "\t", $offset, "\t", $overlap, "\t", $info_cont1, "\t", $info_cont2, "\n"  if $DEBUG;
				if( $divergence < $cutoff )  {
					my $rowcontent = "";
					my $txt_rowcontent = $query_id . "\t" . $consensus{$query_id} . "\t";
					my $cstxt_rowcontent = $query_id . ";" . $consensus{$query_id} . ";";
					#my $source = "TREG";
					## Jaspar link: javascript:Start('http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?ID=MA0101&rm=present')
					my $ref_str;
					my $counter = 303;
					my $x = $subject;
					$x =~ s/^M(0){1,4}//;
					if( $source{$subject} =~ "TRANSFAC" )  { $ref_str = "<a href=\"http://www.biobase.de/cgi-bin/biobase/transfac/8.4/bin/getTFProf.cgi\?$subject\">$subject</a>"; }
					elsif( $source{$subject} =~ "PUBLIC" and $counter <= 302 )  { 
						$ref_str = "<a href=\"http://www.cbil.upenn.edu/cgi-bin/tess/tess?request=MTX-DBRTRV-Accno&key=$subject\">$subject</a>"; 
					}  elsif( $source{$subject} =~ "PUBLIC" and $counter > 302 )  { 
						$ref_str = "<a href=http://www.gene-regulation.com/cgi-bin/pub/databases/transfac/getTF.cgi?AC=\"$subject\">$subject</a>"; 
					} elsif( $source{$subject} =~ "JASPAR" )  { $ref_str = "<a href=\"http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?ID=$subject&rm=present\">$subject</a>"; 
					} else  { 
						if( defined $source{$subject} )   { print STDOUT  "source: ", $source{$subject}, "\n"  if $DEBUG; }
						if( $source{$subject} ne "Query" )  { $source{$subject} = "--"; }
						$ref_str = "-";
					}
					my $d = sprintf("%.3f", $divergence);
					
					$rowcontent .= td($pwm_name{$subject}); $txt_rowcontent .= $pwm_name{$subject} . "\t"; $cstxt_rowcontent .= $pwm_name{$subject} . ";";
					#$rowcontent .= td($class{$subject}); $txt_rowcontent .= $class{$subject} . "\t"; $cstxt_rowcontent .= $class{$subject} . ";";
					$rowcontent .= td($source{$subject}); $txt_rowcontent .= $source{$subject} . "\t"; $cstxt_rowcontent .= $source{$subject} . ";";
					$rowcontent .= td($ref_str); 		$txt_rowcontent .= $subject . "\t"; $cstxt_rowcontent .= $subject . ";";
					$rowcontent .= td($pmat_length{$subject}); $txt_rowcontent .= $pmat_length{$subject} . "\t"; $cstxt_rowcontent .= $pmat_length{$subject} . ";";
					$rowcontent .= td($orientation); 	$txt_rowcontent .= $orientation . "\t"; $cstxt_rowcontent .= $orientation . ";";
					$rowcontent .= td($offset); 		$txt_rowcontent .= $offset . "\t"; $cstxt_rowcontent .= $offset . ";";
					$rowcontent .= td($d); 			$txt_rowcontent .= $d . "\t"; $cstxt_rowcontent .= $d . ";";
					$rowcontent .= td($overlap); 		$txt_rowcontent .= $overlap . "\t"; $cstxt_rowcontent .= $overlap . ";";
					$rowcontent .= td($consensus{$subject}); $txt_rowcontent .= $consensus{$subject} . "\t"; $cstxt_rowcontent .= $consensus{$subject} . ";";
					$tablecontent .= TR($rowcontent); $txt_tablecontent .= $txt_rowcontent . "\n";  $cstxt_tablecontent .= $cstxt_rowcontent . "\n";
					
					my $k = $pwm_name{$subject} . "_" . $orientation;
					$tablecontent_hash{$k} = TR($rowcontent);
					$txt_tablecontent_hash{$k} = $txt_rowcontent . "\n";
					$cstxt_tablecontent_hash{$k} = $cstxt_rowcontent . "\n";
					$divergence_hash{$k} = $d;

				}
			}
		}
	}
	if( $tablecontent ne "" )  {
		#This sorts the keys by their associated values:
		my $tablecontent2 = "";
		my $txt_tablecontent2 = "";
		my $cstxt_tablecontent2 = "";
		foreach  my $tr_key (sort { $divergence_hash{$a} <=> $divergence_hash{$b} } keys %divergence_hash)  {
			$tablecontent2 .= $tablecontent_hash{$tr_key};
			$txt_tablecontent2 .= $txt_tablecontent_hash{$tr_key};
			$cstxt_tablecontent2 .= $cstxt_tablecontent_hash{$tr_key};
		}

		print RES  table( {'border' => 1}, "\n",
			TR( th("Name"), th("Source DB"), th("ID, Link to Source"), th("Length"), th("Orientation"), th("Offset"), 
			    th("Divergence"), th("Overlap"), th("Consensus") ), "\n",
			$tablecontent2 ), "\n";   ### , th("IC 1"), th("IC 2")
		print TXT  $txt_tablecontent2;
		print CS   $cstxt_tablecontent2;
	}  else  { print RES  "<p>No similar matrices found. (cutoff: $cutoff) \n"; }     
}
close(RES);
close(TXT);
close(CS);

print STDOUT  "\nThis is the end.\n";

###########################################################################################################################################
###########################################################################################################################################
=pod

=item I<process_raw_PWMs>

        Based on the function with the same name in pfm_meme_query.pl
	Here we don't write to the MX file but return a reference to the result array.
		process_raw_PWMs( \@rows, \%query_set_ids, \%reg_pmat_TREG, \%RC_reg_pmat_TREG, \%source, \%r_consensus, \%r_RC_consensus, \%r_pmat_length ); 
		# fill %query_set_ids
		# generate reg_pmat
		# fill %reg_pmat_TREG, %RC_reg_pmat_TREG 

=cut

sub process_raw_PWMs  {
	my $ref = shift;
	my @rows = @$ref;
	my $r_query_set_ids = shift;  ## reference to %query_set_ids
	my $r_reg_pmat_TREG = shift;  ## reference to %reg_pmat_TREG
	my $r_RC_reg_pmat_TREG = shift;  ## reference to %RC_reg_pmat_TREG
	my $r_source = shift;   ## Source hash has to filled for self comparison
	my $r_consensus = shift;
	my $r_RC_consensus = shift; #, \%pmat_length
	my $r_pmat_length = shift;
	my @result_ary = ();  

	my $str = join( " ", @rows); 
	my @elems = split /\s+/, $str;
	#while( <@elems> )  { print LOG  ":", $_; }
	#print LOG "\n";
	my ($name, $prob_mat_string) = ("","");
	my $first_elem = "";
	$first_elem = shift @elems;
	if( $first_elem ne "" )  { 
		if( $first_elem =~ m/^>/ )  { $first_elem =~ s/\>//; }
		$name = $first_elem; 
		while( defined $elems[0] and $elems[0] !~ m/\d/ )  { shift @elems; }
		foreach  my $element (@elems)  {
			if( $element =~ m/[a-zA-Z]/ )  {  ### element contains any alphabetic letter
				if( $prob_mat_string ne "" )  { 
					print STDOUT  $name, "##", $prob_mat_string, "\n"  if $DEBUG;
					$r_query_set_ids->{$name} = 'TRUE';
					$r_source->{$name} = 'Query';
					my @cmat = split /\s+/, $prob_mat_string;
					my @reg_pmat =  get_reg_pmat_from_cmat( \@cmat );
					$r_reg_pmat_TREG->{$name} = \@reg_pmat;
					$r_consensus->{$name} = get_iupac_from_pmat(\@reg_pmat);
					$r_pmat_length->{$name} = scalar( @reg_pmat ) / 4;
					my @rc = revcom_matrix( \@reg_pmat );
					$r_RC_reg_pmat_TREG->{$name} = \@rc;
					$r_RC_consensus->{$name} = get_iupac_from_pmat(\@rc);
					#$query_pwm_count++;
					#if( $pwm_set ne "" )  { $pwm_set .= "#"; }
					#$pwm_set .= $name;
					$prob_mat_string = ""; 
					if( $element =~ m/^>/ )  { $element =~ s/\>//; }
					$name = $element;
				}  
			}  else {
				if( $prob_mat_string eq "" )  { $prob_mat_string = $element; }
				else  { $prob_mat_string .= " " . $element; }
			}
		}
		if( $prob_mat_string ne "" )  { 
			print STDOUT  $name, "##", $prob_mat_string, "\n"  if $DEBUG;
			$r_query_set_ids->{$name} = 'TRUE';
			$r_source->{$name} = 'Query';
			my @cmat = split /\s+/, $prob_mat_string;
			my @reg_pmat =  get_reg_pmat_from_cmat( \@cmat );
			$r_reg_pmat_TREG->{$name} = \@reg_pmat;
			$r_consensus->{$name} = get_iupac_from_pmat(\@reg_pmat);
			$r_pmat_length->{$name} = scalar( @reg_pmat ) / 4;
			my @rc = revcom_matrix( \@reg_pmat );
			$r_RC_reg_pmat_TREG->{$name} = \@rc;
			$r_RC_consensus->{$name} = get_iupac_from_pmat(\@rc);
		}  
	}
	return;
}


=pod

=item I<process_meme_PWMs>

        Based on the function with the same name in pfm_meme_query.pl
	Here we don't write to the MX file but return a reference to the result array.
		process_meme_PWMs( \@rows, \%query_set_ids, \%reg_pmat_TREG, \%RC_reg_pmat_TREG, \%source, \%r_consensus, \%r_RC_consensus, \%r_pmat_length ); 
		# fill %query_set_ids
		# generate reg_pmat
		# fill %reg_pmat_TREG, %RC_reg_pmat_TREG 

=cut

sub process_meme_PWMs  {
	my $ref = shift;
	my @rows = @$ref;
	my $r_query_set_ids = shift;  ## reference to %query_set_ids
	my $r_reg_pmat_TREG = shift;  ## reference to %reg_pmat_TREG
	my $r_RC_reg_pmat_TREG = shift;  ## reference to %RC_reg_pmat_TREG
	my $r_source = shift;   ## Source hash has to filled for self comparison
	my $r_consensus = shift;
	my $r_RC_consensus = shift; #, \%pmat_length
	my $r_pmat_length = shift;
	my @result_ary = ();  
	while( my $line = shift(@rows) )  { 
		my $prob_mat_string = "";
		if( $line =~ /position-specific probability matrix/ and $line =~ /Motif/ )  {
    			$line =~ s/\D//g;
    			my $pwm_id = sprintf "MEME%03u", $line;    
    			if( not defined( $line = shift(@rows) ) )  { error("No more lines!\n"); }
    			if( not defined( $line = shift(@rows) ) )  { error("No more lines!\n"); }
			my $no_positions = 0;
			if( $line !~ "letter" )  { 
				if( not defined( $line = shift(@rows) ) )  { error("No more lines!\n"); } 
				if( $line !~ "letter" )  {}
				my @b = split /\s+/, $line;
				$no_positions = $b[5];
			}  ### This should do for html-output files of MEME
    			else  { 
				my @b = split /\s+/, $line;
				$no_positions = $b[5];
			}
			#example of $line: letter-probability matrix: alength= 4 w= 10 n= 54 E= 3.2e-002
			if( $no_positions > 0 )  {
    				for (my $j = 0; $j < $no_positions; $j++) {
					if( not defined( $line = shift(@rows) ) )  { error("No more lines!\n"); }
					chomp $line;
					$line =~ s/^\s+//g;
					if( $prob_mat_string ne "" )  { $prob_mat_string .= " "; }
					$prob_mat_string .= $line;
   				}
				print STDOUT  $pwm_id, "##", $prob_mat_string, "\n"  if $DEBUG;
				$r_query_set_ids->{$pwm_id} = 'TRUE';
				$r_source->{$pwm_id} = 'Query';
				my @cmat = split /\s+/, $prob_mat_string;
				my @reg_pmat =  get_reg_pmat_from_cmat( \@cmat );
				$r_reg_pmat_TREG->{$pwm_id} = \@reg_pmat;
				$r_consensus->{$pwm_id} = get_iupac_from_pmat(\@reg_pmat);
				$r_pmat_length->{$pwm_id} = scalar( @reg_pmat ) / 4;
				my @rc = revcom_matrix( \@reg_pmat );
				$r_RC_reg_pmat_TREG->{$pwm_id} = \@rc;
				$r_RC_consensus->{$pwm_id} = get_iupac_from_pmat(\@rc);
			}
		}
	}
	return;
}


=pod

=item I<compare_partial>

        ($divergence, $offset, $overlap) = $pssm1->compare_partial($pssm2);

The third method to compare two PSSMs. $pssm1 is compared with all
shifts of $pssm2 for which at least 3 positions overlap. The best hit is reported.
The distance is computed utilizing the Formula of the relative entropy (made symetric).
The distance value is is divided by the number of overlapping positions.
The offset holds the shift of pssm1's start position relative to pssm2's start position.

=cut

sub compare_partial {
    ### Input parameters: reg_pmat1, reg_pmat2
    #print STDOUT  @_, "\n";
    my @reg_p_mat1;
    my @reg_p_mat2;
    my $cutoff = 0;
    my $exchanged = "FALSE";
    if( scalar(@{$_[0]}) > scalar(@{$_[1]}) )  {
    	@reg_p_mat2 = @{$_[0]};
	@reg_p_mat1 = @{$_[1]};
	$exchanged = "TRUE";
    }  else  {
    	@reg_p_mat1 = @{$_[0]};
	@reg_p_mat2 = @{$_[1]};    
    }
### ASSERTION: @reg_p_mat1 <= @reg_p_mat2

    #my @bg = (0.25, 0.25, 0.25, 0.25);

    my $len1 = scalar(@reg_p_mat1) / 4;
    my $len2 = scalar(@reg_p_mat2) / 4;
    #print STDOUT  "\ncompare M1 ($len1): ",  @reg_p_mat1, "\n"  if $DEBUG; 
    #print STDOUT  "\ncompare M2 ($len2): ",  @reg_p_mat2, "\n"  if $DEBUG; 
    
    my $min_div = 1e10;
    my $offset = -1000;  ## start position of matrix 1 relative to start position of 2
    my $overlap = 0;
    my $protrusion = 0;  
    my ($ic1, $ic2) = (0,0);
    if( $len1 > 6 )  { $protrusion = int( $len1 / 2 ); }  # maximally half of the shorter matrix protrudes
    elsif( $len1 < 4 )  { $protrusion = 0; }
    else  { $protrusion = $len1 - 4; }  ## len1: 4,5,6
    for( my $idx = -$protrusion; $idx <= $len2 - $len1 + $protrusion; $idx++) {
	my $divergence = 0.0;
	my $overlap_tmp = 0;
    	my $information_content1 = 0;  ## information content of the matching part of matrix 1
    	my $information_content2 = 0;  ## - sum_i( pi * log( pi ) )
	if( $idx < 0 )  {
		my $i2 = 0;
		for( my $i1 = 0-$idx; defined( $reg_p_mat1[4*$i1+3] ); $i1++ )  {
			$overlap_tmp++;
			my @p1 = @reg_p_mat1[(4*$i1)..(4*$i1+3)];
			my @p2 = @reg_p_mat2[(4*$i2)..(4*$i2+3)];
			$divergence += _calc_divergence(\@p1, \@p2);
			$information_content1 -= _calc_information_content(\@p1); 
			$information_content2 -= _calc_information_content(\@p2); 
			$i2++;
		}
	}  else  {
		my $i2 = $idx;
		for( my $i1 = 0; defined( $reg_p_mat1[4*$i1+3] ) and defined( $reg_p_mat2[4*$i2+3] ); $i1++ )  {
			$overlap_tmp++;
			my @p1 = @reg_p_mat1[(4*$i1)..(4*$i1+3)];
			my @p2 = @reg_p_mat2[(4*$i2)..(4*$i2+3)];
			$divergence += _calc_divergence(\@p1, \@p2); 
			$information_content1 -= _calc_information_content(\@p1); 
			$information_content2 -= _calc_information_content(\@p2); 
			$i2++;
		}
	}
		
	if( $overlap_tmp > 0 and ($divergence/$overlap_tmp) < $min_div and 
		($information_content1 < $overlap_tmp or $information_content2 < $overlap_tmp ) )  {  
	    $min_div =  $divergence / $overlap_tmp;
	    $offset = $idx;
	    $overlap = $overlap_tmp;
	    $ic1 = $information_content1;
	    $ic2 = $information_content2;
	}
    }
    if( $exchanged eq "TRUE" )  {
    	if( $offset != 0 )  { $offset = -$offset; }
    }
    return( $min_div, $offset, $overlap, $ic1, $ic2 );
}


=pod

=item I<_calc_divergence>

        $div = _calc_divergence(\@a1, \@a2);

Takes references to two arrays which should have equal lengths and
represent probability distributions (positive entries which sum up to
1, length can be arbitrary). Returns the divergence, i.e. the sum of
the respective relative entropies.

=cut


sub _calc_divergence {
    my ($ref1, $ref2)= @_;

    my $epsilon = 0.02;

    my $len1 = @$ref1;
    my $len2 = @$ref2;

    $len1 == $len2 or error("arrays must have equal lengths in _calc_divergence\n");

    my $div = 0;
    my $probsum1 = 0;
    my $probsum2 = 0;
    for(my $j=0; $j < $len1; $j++) {
	($ref1->[$j] > 0 && $ref2->[$j] > 0) or
	  error("your probability distribution has non-positive entries in _calc_divergence\n");

	$probsum1 += $ref1->[$j];
	$probsum2 += $ref2->[$j];
	$div += $ref1->[$j] * log($ref1->[$j]/$ref2->[$j])
	  + $ref2->[$j] * log($ref2->[$j]/$ref1->[$j]);
    }

#    print "sum 1: $probsum1\tsum 2: $probsum2\n";
    (abs($probsum1 - 1) < $epsilon && abs($probsum2 - 1) < $epsilon) or
      error("Calculated divergence from arrays which don't look like probability distributions in _calc_divergence! ($probsum1, $probsum2)\n");

    return($div);
}

=pod

=item I<_calc_information_content>

        $div = _calc_information_content(\@a1);

Takes references to an array...

=cut


sub _calc_information_content {
    my ($ref)= @_;

    my $epsilon = 0.02;
    my $len = @$ref;

    my $info = 0;
    my $probsum = 0;
    for( my $j=0; $j < $len; $j++ )  {
	($ref->[$j] > 0 && $ref->[$j] > 0) or
	  error("your probability distribution has non-positive entries in _calc_information_content\n");
	$probsum += $ref->[$j];
	$info += $ref->[$j] * log($ref->[$j]);
    }

    (abs($probsum - 1) < $epsilon) or
      error( "Calculated divergence from arrays which don't look like probability distributions in _calc_information_content! ($probsum)\n");

    return($info);
}


sub error  {
	my ($r) = @_;
	print STDOUT  "Error: $r\n";
	exit;
}
	

