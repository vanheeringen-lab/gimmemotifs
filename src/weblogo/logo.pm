#!/usr/local/bin/perl -w

=head1 NAME

    logo.pm - organizes data in FASTA and CLUSTAL formats into height data.

=head1 SYNOPSIS

    Perl module

=head1 DESCRIPTION

    logo.pm: Takes in strings of aligned sequences and sorts them vertically
             based on height as assigned by the following equations found in
             Schneider and Stephens paper "Sequence Logos: A New Way to Display
             Consensus Sequences":
    
                 height = f(b,l) * R(l)                            (1)

             where f(b,l) is the frequency of base or amino acid "b" at position
             "l". R(l) is amount of information present at position "l" and can
             be quantified as follows:

                 R(l) for amino acids   = log(20) - (H(l) + e(n))    (2a)
                 R(l) for nucleic acids =    2    - (H(l) + e(n))    (2b)

             where log is taken base 2, H(l) is the uncertainty at position "l", 
             and e(n) is the error correction factor for small "n". H(l) is
             computed as follows:

                 H(l) = - (Sum f(b,l) * log[ f(b,l) ])             (3)

             where again, log is taken base 2. f(b,l) is the frequency of base
             "b" at position "l". The sum is taken over all amino acids or
             bases, depending on which the data is.

             Currently, logo.pm uses an approximation for e(n), given by:

                 e(n) = (s-1) / (2 * ln2 * n)                      (4)

             Where s is 4 for nucleotides, 20 for amino acids ; n is the number
             of sequences in the alignment. e(n) also  gives the height of error
             bars.

=cut

package logo;

use strict;
use Carp;

################################################################################
######                             SOME VARIABLES                         ######
################################################################################

my $DEBUG = 0;

my $AA = 0;
my $NA = 1;

my %BASES = ("a" => "adenine",
	     "t" => "thymine",
	     "g" => "guanine",
	     "c" => "cytosine",
	     "u" => "uracil");

# does not include B or Z
my %AMINOACIDS = ("a" => "", "c" => "", "d" => "", "e" => "", "f" => "",
		  "g" => "", "h" => "", "i" => "", "k" => "", "l" => "",
		  "m" => "", "n" => "", "p" => "", "q" => "", "r" => "",
		  "s" => "", "t" => "", "v" => "", "w" => "", "y" => "");

my @data;
my $kind;
my ($seqs_r, $desc_r);

my $CONFIDENCE_LIMIT = 0.90;

################################################################################
######                             SOME FUNCTIONS                         ######
################################################################################

=head1 APPENDIX

=cut

=head2 getHeightData()

 Usage   : my ($height_data_r, $description_data_r, $kind) =
              logo::getHeightData($input_data_r, $params);
 Returns : * REFERENCE TO array of height data
           * REFERENCE TO array of input description strings
           * $AA if the data is amino acid, $NA otherise
 Args    : $input_data_r : input data in CLUSTAL or FASTA formats
         : $params       : hash of parameters

 getHeightData is the entry point into the logo.pm module. $input_data_r is a
 reference to  an array of strings containing FASTA or CLUSTAL data, where all
 lines whose first character is "#" is considered a comment line.

 $params is a hash of parameters with the following keys:
   * smallsampletoggle : 0 to turn off small sample correction, otherwise
                         small sample correction is on
   * input_kind : 0 for amino acids, 1 for nucleic acids; if undefined,
                  logo.pm will attempt to automatically detect whether the
                  input consists of amino or nucleic acid data. If
                  $input_kind is defined, only those residues  defined by
                  $input_kind will be in the output -- all other residues are
                  considered as spaces. For example, if $input_kind is $NA,
                  the residue "I" or "i" are considered spaces, since "I" and
                  "i" are not nucleic acid residues.
   * stretch : stretch all characters so they are flush at the maximum number
               of bits allowed

 Sample use:

  # get FASTA data
  open (FASTA, "$fastafile");
  my @inputdata = <FASTA>;
  close (FASTA);

   my %heightparams = (
		       smallsamplecorrection => 0,
		       input_kind => 0,
		       stretch => 0
		       );

  # get height data
  my ($heightdata_r, $desc_r, $kind) = logo::getHeightData(\@inputdata, \%heightparams);

=cut

# entry point into module
sub getHeightData {

    # $smallsampletoggle is toggle to turn off small sample correction (1 to turn off)
    # $input_kind can be $AA or $NA or undef
    my ($input_data_r, $params) = @_;

    # gary 040119: adjust for formats (Unix is \n, Mac is \r, Windows is \r\n)
    $input_data_r = normalizeData($input_data_r);

    # gets sequences, sets $kind temporarily
    my ($goodlength, $maxreslength, $badline, $validformat);
    ($seqs_r, $desc_r, $maxreslength, $goodlength, $badline, $validformat) =
	getSeqs($input_data_r, $params->{input_kind});

#    for(my $i = 0; $i < scalar @$seqs_r ; $i++) {
#	print STDERR ($desc_r->[$i] . "\n" . $seqs_r->[$i] . "\n");
#    }
#    print STDERR "maxreslength = $maxreslength\n";
#
#    exit(0);

    if ($DEBUG) { print STDERR ("point 1\n");}

    # check for valid format
    if ((defined $validformat) && ($validformat == 1)) {
#	print("returning\n");
	return (undef, undef, undef, undef, undef, 1);
    }

    if ($DEBUG) { print STDERR ("point 2\n");}

    # check for bad length
    if (!$goodlength) {
	return (undef, undef, undef, $goodlength, $badline);
    }

    # reset $kind if in $input_kind
    if (defined $params->{input_kind} && isLegalKind($params->{input_kind}) ) {
	$kind = $params->{input_kind};
    }

    # build data
    buildData(@$seqs_r, $params->{smallsampletoggle}, $params->{stretch}, $maxreslength);

    if ($DEBUG) { print STDERR ("point 3\n");}

#    print STDERR ("data size = ", scalar @data, "\n");
#    foreach (@data) {
#	print STDERR ("$_\n");
#    }
#
#    exit(0);
#
#    print STDERR ("return at 2\n");
    return (\@data, $desc_r, $kind, $goodlength, $badline);
}

sub isLegalKind {
    return ($_[0] =~ /^[01]$/);
}

################################################################################
#
# sub normalizeData($data_r) returns $data_r, with Mac/Unix/Windows newline
# style normalized to standard Unix-style newline style
#
################################################################################
sub normalizeData {
    my ($data_r) = @_;

    # check args
    if (not defined $data_r) {
	die "data_r must be defined\n";
    }

    my @normalized = ();
    foreach my $pseudo_line (@$data_r) {
	my @split_line = split(/[\r\n]+/, $pseudo_line);
	push(@normalized, @split_line);
    }

    return \@normalized;
}


################################################################################
#
# sub getSeqs($data_r, $kind) returns 5 values:
#
# * array reference to sequence strings
# * array reference to sequence names
# * length of sequence
# * 1 if all sequences have the same length, 0 else
# * line number L where sequenceLength(L) != sequenceLength(other lines), else
#   undef
#
################################################################################
sub getSeqs {
    my ($input_data_r, $kind) = @_;

    unless( $input_data_r->[0] ){
        return (undef, undef, undef, undef, undef, 1);
    }

    # skip all comment chars and lines of all spaces
    while ( ($input_data_r->[0] =~ /^\s*\#/) || ($input_data_r->[0] =~ /^\s*$/) ) {
	shift @$input_data_r;
	if( !defined $input_data_r->[0])  {return (undef, undef, undef, undef, undef, 1);}
    }

    if (isFormat_FASTA($input_data_r)) {
	return getSeqs_FASTA($input_data_r, $kind);

    } elsif (isFormat_CLUSTAL($input_data_r)) {
	return getSeqs_CLUSTAL($input_data_r, $kind);

    } elsif (isFormat_FLAT($input_data_r)) {
	return getSeqs_FLAT($input_data_r, $kind);

    } else {
	if ($DEBUG) {print STDERR ("format nothing\n");}
	return (undef, undef, undef, undef, undef, 1);
    }

#    if ($_[0] =~ />/) {
#	return getSeqs_FASTA(@_);
#    } else {
#	return getSeqs_CLUSTAL(@_);
#    }
}

################################################################################
#
# sub isFormat_FASTA($data_r) returns 1 if $data_r is in FASTA format
#
################################################################################s
sub isFormat_FASTA {
    my ($input_data_r) = @_;

    # check args
    if (not defined $input_data_r) {
	Carp::confess("logo::isFormat_FASTA : input_data_r must be defined\n");
    }

    if ($input_data_r->[0] =~ />/) {
	return 1;
    } else {
	return 0;
    }
}

################################################################################
#
# sub isFormat_CLUSTAL($data_r) returns 1 if $data_r is in CLUSTAL format
#
################################################################################
sub isFormat_CLUSTAL {
    my ($input_data_r) = @_;

    # check args
    if (not defined $input_data_r) {
	Carp::confess("logo::isFormat_CLUSTAL : input_data_r must be " .
		      "defined\n");
    }

    my $i=0;

#    # skip spaces or just "*" and "." and ":"
#    while ($input_data_r->[$i] =~ /^[\*\:\s]*$/) {
#	$i++;
#    }

    # if it looks like CLUSTAL W (version) ... , then it must be clustal
    if ($input_data_r->[$i] =~ /^\s*CLUSTAL/) {
	return 1;
    }

    # CLUSTAL looks like: "name        seq"
    if ($input_data_r->[$i] =~ /^\s*(\S+)\s+(\S+)\s*$/) {
	return 1;
    } else {
	return 0;
    }
}

################################################################################
#
# sub isFormat_FLAT($data_r) returns 1 if $data_r is in FLAT format
#
################################################################################
sub isFormat_FLAT {
    my ($input_data_r) = @_;

    # check args
    if (not defined $input_data_r) {
	Carp::confess("logo::isFormat_FLAT : input_data_r must be defined\n");
    }

#    print("checking flat\n");
#    print("first element = -->", $input_data_r->[0], "<--\n");

    if ($input_data_r->[0] =~ /^[a-zA-Z\-]+\s*$/) {
	return 1;
    } else {
	return 0;
    }
}

################################################################################
######                          FORMATTING FUNCTIONS                      ######
################################################################################

# the flat sequence format is as follows:
# sequence1
# sequence2
# sequence3
# ...
# sequenceN
sub getSeqs_FLAT {

    if ($DEBUG) {print STDERR "DOING FLAT\n";}

    my ($input_data_r, $input_kind) = @_;

    my $linelength = 0;
    my $seqCount = 0;
    my $total_residues = 0;
    my (@returnVal, @desc) = ();
    my $prevlinelength = undef;
    my $NA_count = 0;

    foreach my $seq (@$input_data_r) {
#	chomp $seq;
	$seq =~ s/\s+$//;

	my @chars = split(//,$seq);

	my $char;
	foreach (@chars) {
	    $total_residues++;
	    $linelength++;

	    # set $char
	    if (defined $input_kind) {
		if ($input_kind == $AA) {
		    $char = (isAA($_)) ? $_ : "-";
		} else { # == $NA
		    $char = (isNA($_)) ? $_ : "-";
		}
	    } else {
		$char = $_;
		if (isNA($char)) {
		    $NA_count++;
		}
	    }

	    $returnVal[$seqCount] .= $char;
	}
	$desc[$seqCount] = "no name";

	if ($seqCount == 0) {
	    $prevlinelength = $linelength;
	} elsif ($prevlinelength != $linelength) {  # different number of residues, so complain
	    return (undef, undef, undef, 0, $seq);  # 0 for not same length, $seq is name
	}
	$linelength=0;

	$seqCount++;
    }

    # determine whether to use $NA or $AA
    if (!defined $input_kind) {
	if ($NA_count / ($total_residues+1) >= $CONFIDENCE_LIMIT) {
	    $kind = $NA;
	} else { 
	    $kind = $AA;
	}
    }

    return (\@returnVal, \@desc, $prevlinelength, 1, undef);
}

sub getSeqs_CLUSTAL {

    if ($DEBUG) {print STDERR "DOING CLUSTAL\n";}

    my ($input_data_r, $input_kind) = @_;

    my @returnVal;
    my @desc;
    my $seqCount=0;
    my $reslength = 0;
    my ($name, $seq);

#    my $input_kind = pop @_;
#    my $CONFIDENCE_LIMIT = 0.90;
    my $NA_count = 0;
    my $total_residues = 0;
    my ($prevlinelength, $linelength) = (0,0);

#    foreach (@_) {
    foreach (@$input_data_r) {
#	chomp;

	if ($DEBUG) {print STDERR ("line = $_\n")};

	$_ =~ s/\s+$//;

	# skip if it is a comment character -- first character is "#"
	next if (/^\s*\#/);

	# skil if it is a CLUSTAL W header line
	next if (/^\s*CLUSTAL/);

	# if spaces or just "*" and "." and ":"
	if (/^[\*\.\:\s]*$/) {
	    $seqCount=0;
	    $prevlinelength=0;
	    next;
	}

	($name,$seq) = (/^\s*(\S+)\s+(\S+)\s*$/);

	if ($DEBUG) { print STDERR ("name, seq = $name, $seq\n"); }

	# add new entry
	if (!defined $desc[$seqCount]) {
	    $desc[$seqCount] = $name;
	    $returnVal[$seqCount] = "";
	}
    
    if(!defined $seq) {return (undef, undef, undef, undef, undef, 1);} # Something has gone terrible wrong
    
	my @chars = split(//,$seq);
	my $char;
	foreach (@chars) {
	    if ($seqCount == 0) {
		$reslength++;     # all sequences have same residue length, so only count first one
	    }

	    $total_residues++;
	    $linelength++;

	    # set $char
	    if (defined $input_kind) {
		if ($input_kind == $AA) {
		    $char = (isAA($_)) ? $_ : "-";
		} else { # == $NA
		    $char = (isNA($_)) ? $_ : "-";
		}
	    } else {
		$char = $_;
		if (isNA($char)) {
		    $NA_count++;
		}
	    }
	    
	    $returnVal[$seqCount] .= $char;
	}

	if ($seqCount == 0) {
	    $prevlinelength = $linelength;
	} elsif ($prevlinelength != $linelength) {  # different number of residues, so complain
	    return (undef, undef, undef, 0, $name);
	}
	$linelength=0;

	$seqCount++;
    }

    # determine whether to use $NA or $AA
    if (!defined $input_kind ) {
	if ($NA_count / ($total_residues+1) >= $CONFIDENCE_LIMIT) {
	    $kind = $NA;
	} else { 
	    $kind = $AA;
	}
    }

    return (\@returnVal, \@desc, $reslength, 1, undef);
}

# if $input_kind is defined, residues that are not defined are set to space
sub getSeqs_FASTA {

    if ($DEBUG) {print STDERR "DOING FASTA\n";}

    my ($input_data_r, $input_kind) = @_;

    my @returnVal;
    my @desc;
    my $count=-1;
    my $newElem=0;

#    my $input_kind = pop @_;

#    my $CONFIDENCE_LIMIT = 0.90;
    my $NA_count = 0;
    my $total_residues = 0;
    my $reslength = 0;
    my $maxreslength = 0;
    
    my ($goodlength, $currline, $prevline);


#    # skip all lines that are all spaces
#    while ($_[0] =~ /^\s*$/) {
#	shift @_;
#    }

#    foreach (@_) {
    foreach (@$input_data_r) {

	$_ =~ s/\s+$//;

	# skip if it is a comment character -- first character is "#"
	next if (/^\s*\#/);

	# skip all lines that are all spaces
	next if (/^\s*$/);

	$_ =~ s/\s+$//;  # cut trailing white space
	$_ =~ s/^\s+//;  # cut leading white space
	if (/>/) {
	    $currline = $_;
	    ($desc[scalar @desc]) = ($_ =~ />\s*(.+)$/);

	    if (not $newElem) {		
		$count++;
		$newElem = 1;
	    }
	} else {
	    if ($newElem) {
		$maxreslength = $reslength if $maxreslength == 0;
		if (($maxreslength != 0) && ($maxreslength != $reslength)) {
		    return (undef, undef, undef, 0, $prevline);
		}

		$maxreslength = $reslength;
		$reslength = 0;
	    }

	    my @chars = split(//,$_);
	    my $char;
	    foreach (@chars) {
		$reslength++;
		$total_residues++;

		# set $char
		if (defined $input_kind) {
		    if ($input_kind == $AA) {
			$char = (isAA($_)) ? $_ : "-";
		    } else { # == $NA
			$char = (isNA($_)) ? $_ : "-";
		    }
		} else {
		    $char = ($_ =~ /[a-zA-Z]/) ? $_ : "-";  # if not alpha char, use space
		    if (isNA($char) && !isSpace($char)) {
			$NA_count++;
		    }
		}

		if ($newElem) {
		    $returnVal[$count] = $char;
		} else {
		    $returnVal[$count] .= $char;
		}
		$newElem = 0;
	    }
	    $prevline = $currline if $currline =~ />/;
	}
    }

    # check if last is biggest
    if (($maxreslength != 0) && ($maxreslength != $reslength)) {
	return (undef, undef, undef, 0, $prevline);
    }
#    $maxreslength = ($reslength > $maxreslength) ? $reslength : $maxreslength;

    # determine whether to use $NA or $AA
    if (!defined $input_kind) {
	if ($NA_count / ($total_residues+1) >= $CONFIDENCE_LIMIT) {
	    $kind = $NA;
	} else { 
	    $kind = $AA;
	}
    }

    return (\@returnVal, \@desc, $maxreslength || $reslength, 1, undef);  # 1 for good lengths
}

sub isSpace {
    return $_[0] =~ /[ \-]/;
}

sub isAA {
    return (defined $AMINOACIDS{lc $_[0]});
}

sub isNA {
    return (defined $BASES{lc $_[0]});
}

################################################################################
######                       DATA BUILDING FUNCTIONS                      ######
################################################################################


# arguments: takes reference to array and lines aligned sequences of bases or
#            amino acids
# returns: updated reference array to reflect contents of sequences sorted
#          vertically by height as described by (1)
sub buildData {
    
    my $currentx = 0;
    my $h;
    my $count=0;
    my $maxreslength = pop (@_);
    my $stretch = pop(@_);
    my $smallsampletoggle = pop (@_);
    my $totalsize = $#_+1;

    while ($currentx < $maxreslength) {       #(length $_[0])) {
	my $allspaces = 1;
	my $spaceCount=0;

	# get vertical sequence
	my @vert=();
	foreach (@_) {  # foreach sequence
	    my $currentchar;

	    # set currentchar, set to " " if $_ is not long enough
	    if ($currentx >= (length $_)) {
		$currentchar = " ";
	    } else {
		$currentchar = substr($_,$currentx,1);
	    }

	    # in all sequences, "-" is considered as a space
	    # don't count " " and "-"
	    if (($currentchar ne "-") && ($currentchar ne " ")) {
		$vert[scalar @vert] = uc substr($_,$currentx,1);
		$allspaces = 0;
	    } else {
		$spaceCount++;
	    }
	}

	if ($allspaces) {
	    # build @vert
	    @vert = (" 0", ">0");

	    # place in @data
	    $data[scalar @data] = \@vert;

	    $currentx++;
	    next;
	}

	my $error;
	if ($smallsampletoggle) {
	    $error=getError($kind,$totalsize - $spaceCount);
	} else {
	    $error = 0;
	}

	# sort vertical sequence by amino acid or base
	@vert = sort(@vert);
	my $total = $#vert + 1;

	# find H(l) -- must be done before collapsing
	$h = getH(@vert);

	# collect like terms
	@vert = collapse(@vert);

	# get R(l)
	my $r;
	if (!defined $stretch || !$stretch) {
	    $r= getR($kind, $h, $error);
	} else {
	    $r = ($kind == $NA) ? 2 : (log(20) / log(2));
	}

	# place heights
	my $count=0;
	my $height;
	my $elem;
	foreach $elem (@vert) {
	    $height = getHeight(substr($elem, 2) / $total,$r);
	    $vert[$count] = substr($elem,0,1) . $height;
	    $count++;
	}

	# sort by height
	@vert = sort height_sort @vert;

	# put in error bar size
	$vert[$count] = ">$error";

	# place in @data
	$data[scalar @data] = \@vert;

	$currentx++;
    }
}

# uses error approximation given by:
#             e :=  (s-1) / (2 * ln2 * ntrue);
sub getError {
    return ((($_[0] == $NA) ? 4 : 20) - 1) / (2 * log(2) * $_[1]);
}

sub height_sort {
    my ($lettera, $heighta) = ($a =~ /^(.{1})(\S+)$/); #substr($a,1);
    my ($letterb, $heightb) = ($b =~ /^(.{1})(\S+)$/); #substr($b,1);
    
    # compare by height first, then letter
    if ($heighta <=> $heightb) {
	return $heighta <=> $heightb;
    } else {
	return $letterb cmp $lettera;  #reversed for some reason...
    }
}

sub collapse {
    my @returnVal;
    my $current = $_[0];
    my $count=0;
    my $freq;

    foreach (@_) {
        if ($current eq $_) {
            $count++;
        } else {
	    $returnVal[scalar @returnVal] = "$current $count";

            $current = $_;
            $count=1;
        }
    }

    # add last element
    $returnVal[scalar @returnVal] = "$current $count";

    return @returnVal;
}

# arguments : $_[0] : list of bases or amino acids
sub getH {
    my $h = 0;
    my (@vert) = @_;  # vertical sequence (comparing multiple sequences)

    my $current = $vert[0];
    my $count=0;
    my $freq;

    foreach (@vert) {
	if ($current eq $_) {
	    $count++;
	} else {
	    $freq = $count / ($#vert + 1);
	    $h += $freq * (log ($freq) / log(2));

	    $current = $_;
	    $count=1;
	}
    }

    # add last element
    $freq = $count / ($#vert + 1);
    $h += $freq * (log ($freq) / log(2));

    return -$h;
}

# arguments : $_[0] : $AA or $NA
#             $_[1] : uncertainty = H(l)
#             $_[2] : error correction factor for small number of sequences
#                     = e(n) ; assummed to be 0 if not given
sub getR {
    my $max = ($_[0] == $NA) ? 2 : (log(20) / log(2));
    my $e = (defined $_[2]) ? $_[2] : 0;
    return ($max - ($_[1] + $e));
}

# arguments: $_[0] : frequency = f(b,l)
#            $_[1] : R(l)
sub getHeight {
    return $_[0] * $_[1]
}

1;
