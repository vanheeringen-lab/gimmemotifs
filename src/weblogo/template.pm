=head1 NAME

  template.pm - creates logo output in various formats

=head1 SYNOPSIS

  Perl module 

=head1 DESCRIPTION

  logo.cgi and run.pl collect the logo data. They can then enter
  template::create_template to create logo output in the following formats:

   * EPS
   * GIF
   * PDF
   * PNG

  

  If the configuration file "logo.conf" exists in the working directory, then
  it will be parsed for the locations of GhostScript (gs) and convert. The
  following is an
  example of the configuration file ("#" at beginning of line indicates
  comment):

    # Make configuration changes here. Rename this file to logo.conf
    # gs version 5.5 does not work, 6.5 does
    # set the PATHS
    gs=/usr/local/bin/gs
    convert=/usr/X11R6/bin/convert

=cut

package template;
use strict;


################################################################################
########                   STATE FOR FILLINGS                         ##########
################################################################################

#        BARBITS            number of bits in vertical y-axis bar
#        CHARWIDTH          width of characters in logo
#        COLORSCHEME        "a" for amino acid;
#                           "n" for nucleic acid;
#                           "b" for black (no color scheme)
#                           color scheme
#        DATA               string of heights of characters in cm
#        DESC               description of aligned sequences
#        FINEPRINT          enable adverts/credits
#        LOGOLINES          number of lines of logo
#        LOGOHEIGHT         height of EACH LINE in logo (in cm)
#        LOGOWIDTH          width of final logo in cm
#        LOGOHEIGHTPOINTS   height of final logo in points
#        LOGOWIDTHPOINTS    height of final logo in points
#        LOGOSTART          logo output will begin at residue LOGOSTART
#        LOGOEND            logo output will end at residue LOGOEND
#        ERRBAR             1 to include error bar, 0 to exclude
#        ERRORBARFRACTION   percent of error bar to show in range [0,1]
#        KIND               $AA for amino acid, $NA for nucleic acid
#        NUMBERING          1 to show residue numbers, 0 to exclude
#        OUTLINE            1 to print characters in outline form;
#                           0 to print characters in solid form
#        SHOWENDS           "d" to show 5' and 3' ends;
#                           "p" to show N and C termini;
#                           "-" to exclude end markers
#        SHOWINGBOX         "n" to have No boxes around characters;
#                           "s" to have boxes around characters, with Shrinking;
#                           "f" to have Filled boxes around characters
#        SHRINKBOOLEAN      1 to shrink characters; 0 to exclude shrinking
#        SHRINKFACTOR       amount to shrink in range from 1(no shrinking) to
#                           0(full shrinking)
#        START_NUM          start number for very beginning of sequence
#        TITLE              title of logo
#        YAXIS              1 to turn on y axis and its labels

################################################################################
#####                   VARIABLES AND DEFAULT VALUES                       #####
################################################################################

my %defaults = (
		LOGOHEIGHT => 5,
		LOGOWIDTH => 8,

		YAXIS => "false",
		TITLE => "",
		YAXIS_LABEL => "bits",
		XAXIS_LABEL => "",

		BARENDS => "false",
		OUTLINE => "false",
		SHOWINGBOX => "n",
		NUMBERING => "false",
        #FINEPRINT => "Created by: G. E. Crooks, G. Hon, J.-M. Chandonia & S. E. Brenner, (2002) <weblogo.berkeley.edu>",
        FINEPRINT => "weblogo.berkeley.edu",

		ERRORBARFRACTION => "1",
		ERRBAR => "0",
		SHRINKFACTOR => "1",
		START_NUM => "1",

		DEFAULT_COLOR => "black",

		black  => "[0 0 0]",
		red    =>  "[0.8 0 0]",
		green  =>  "[0 0.8 0]",
		blue   =>  "[0 0 0.8]",
		yellow =>  "[1 0.7 1.0]",
		purple =>  "[0.8 0 0.8]",
		orange =>  "[1 0.7 0]"
		);

my $AA = 0;
my $NA = 1;
my $PATH;

################################################################################
#####                          SOME FUNCTIONS                              #####
################################################################################


sub create_template {
    my ($input, $kind, $desc_r, $data_r, $outfile, $path) = @_;

    # set path
    $PATH = $path;

    #Create EPS file
    my %fillings;

    # put parameters in fillings
    makeFillings(\%fillings, $input, $kind, $desc_r, $defaults{FINEPRINT});

    # set default data if not filled
    setDefaults(\%fillings, \%defaults, scalar @$data_r);

    # put color in fillings
    setColors(\%fillings, \%defaults);

    # put data in fillings
    setData(\%fillings, $data_r);

    # make eps output
    my $eps = fillTemplate("$PATH/template.eps", \%fillings);
    my $format = $input->{FORMAT};

    # convert
    my ($gsprog, $convertprog) = getProgs();

#    print STDERR "(gsprog, convertprog) = ($gsprog, $convertprog)\n";

    my $width = $fillings{LOGOWIDTHPOINTS};
    my $height = $fillings{LOGOHEIGHTPOINTS};       # height of entire logo
    my $res = $input->{RES};
    my $antialias = (defined $input->{ANTIALIAS} && $input->{ANTIALIAS}) ? "-dTextAlphaBits=4" : "";

    my $r = $width . "x" . $height;

    if( $format eq "EPS" ) {
	if ($outfile eq "-") { # if standard out
	    print $eps;
	} else {
	    open (WRITEME, ">$outfile");
	    print WRITEME $eps;
	    close (WRITEME);
	}

    } elsif ($format eq "PDF"){
#	print("outfile = $outfile\n");
	my $program = "| $gsprog -sOutputFile=$outfile -sDEVICE=pdfwrite -dPDFSETTINGS=/printer -q -r$res -dDEVICEWIDTHPOINTS=$width -dDEVICEHEIGHTPOINTS=$height -dEmbedAllFonts=true $antialias -dSAFER -dBATCH  -dNOPAUSE -_";
	open(WRITEME, $program);
	print WRITEME "$eps";
	close (WRITEME);

    } elsif ( $format eq "PNG" ) {
	my $program = "| $gsprog -sOutputFile=$outfile -sDEVICE=png16m -q -r$res -dDEVICEWIDTHPOINTS=$width -dDEVICEHEIGHTPOINTS=$height $antialias -dSAFER -dBATCH  -dNOPAUSE -_";
    #print ("$program");
	open(WRITEME, $program);
	print WRITEME "$eps";
        close (WRITEME);

    } elsif ($format eq "GIF") {
	# convert to EPS first, then GIF
        die "Please check logo.conf: convert program does not exist" 
            if (!defined $convertprog || !(-e $convertprog));

	my $program = "| $gsprog -sOutputFile=- -sDEVICE=png16m -q -r$res -dDEVICEWIDTHPOINTS=$width -dDEVICEHEIGHTPOINTS=$height $antialias -dSAFER -dBATCH  -dNOPAUSE -_";

	if ($outfile eq "-") {
	    $program .= " | $convertprog png:- gif:-";
	} else {
	    $program .= " | $convertprog png:- $outfile";
	}

	open(WRITEME, $program);
	print WRITEME "$eps";
        close (WRITEME);
    }
}

#deprecated
sub c {
    return create_template( @_);
}

sub getProgs {
    my ($gsprog, $convertprog) = ("gs", "convert");     


    # No configuration file, then return defaults.
    return ($gsprog, $convertprog) if (! (-e "$PATH/logo.conf"));
    open (CONF, "$PATH/logo.conf");
    
    while (<CONF>) {
	next if (/^\#/);  # skip lines beginning with "#"
	if (m/^gs/i) {    # if looks like gs (case insensitive)
	    ($gsprog) = ($_ =~ /^\S+\=(.+)$/);
	}
	if (m/^convert/i) { # if looks like convert (case insensitive)
	    ($convertprog) = ($_ =~ /^\S+\=(.+)$/);
	}
    }

    # Do these fings exist?
    my ($gsprogname) = ($gsprog =~ /^(\S+)/);

    die "Please check $PATH/logo.conf: gs program ($gsprogname) does not exist" if (!defined $gsprogname || !(-e $gsprogname));
    #die "Please check logo.conf: convert program does not exist" if (!defined $convertprog || !(-e $convertprog));

    return ($gsprog, $convertprog);
}

sub get_eps { 
    my ($input, $kind, $desc_r, $data_r) = @_;
    my %fillings;

    # put parameters in fillings
    makeFillings(\%fillings, $input, $kind, $desc_r);

    # set default data if not filled
    setDefaults(\%fillings, \%defaults, $#$data_r);

    # put data in fillings
    setData(\%fillings, $data_r);

    # make output
    return fillTemplate("$PATH/template.eps", \%fillings);
}

sub fillTemplate {
    my ($filename, $fillings) = @_;

    if (not -e $filename) {
	die "filename ($filename) must exist\n";
    }

    my $text;
    local $/; # slurp mode (undef)
    local *F; # create local filehandle
    open(F, "< $filename\0") || return;
    $text = <F>;
    close(F);

    #replace {$KEYWORDS} with value in %$fillings hash
    $text =~ s{ \{\$ (.*?) \} }
    { exists( $fillings->{$1})
	  ? $fillings->{$1}
                     : ""
		     }gsex;
    return $text;
}


################################################################################
#####                    FILL THE FILLINGS HERE                            #####
################################################################################

sub isChecked {
    return 0 if (!defined $_[0]);
    return $_[0];
}

# negative/positive ints
sub isInt {
    return ($_[0] =~ /^[-\+]?\d+$/) ? 1 : 0;
}

sub makeFillings {

#    my ($fillings, $input, $kind, $desc_r, $data_r, $fineprint) = @_;
    my ($fillings, $input, $kind, $desc_r, $fineprint) = @_;
    $fillings->{KIND} = $kind;
    $fillings->{LOGOHEIGHT} = $input->{LOGO_HEIGHT};
    $fillings->{LOGOWIDTH} = $input->{LOGO_WIDTH};
    $fillings->{OUTLINE} = (isChecked($input->{OUTLINE})) ? "true" : "false";
    $fillings->{NUMBERING} = (isChecked($input->{NUMBERING})) ? "true" : "false";
    $fillings->{FINEPRINT} = (isChecked($input->{FINEPRINT})) ? $fineprint : "";

    $fillings->{LOGOSTART} = $input->{LOGOSTART};
    $fillings->{LOGOEND} = $input->{LOGOEND};
    $fillings->{START_NUM} = $input->{START_NUM};

    $fillings->{YAXIS} = (isChecked($input->{YAXIS}) && !isChecked($input->{STRETCH})) ? "true" : "false";
    
    
    
    $fillings->{TITLE} = $input->{TITLETEXT};
    $fillings->{YAXIS_LABEL} = $input->{YAXIS_LABEL};

    $fillings->{XAXIS_LABEL} = $input->{XAXIS_LABEL};
    $fillings->{ERRBAR} = $input->{ERRBAR};
    $fillings->{SHOWINGBOX} = (isChecked($input->{SHOWINGBOX})) ? "s" : "n";
    $fillings->{SHRINKBOOLEAN} = ($fillings->{SHOWINGBOX} eq "s") ? "true" : "false";
    $fillings->{SHRINKFACTOR} = $input->{BOXSHRINK};

    if ((defined $input->{CHARSPERLINE}) && 
	isInt($input->{CHARSPERLINE}) && 
	($input->{CHARSPERLINE} > 0)) {
	$fillings->{CHARSPERLINE} = $input->{CHARSPERLINE};
    }

    if (defined $input->{BARBITS}) {
	$fillings->{BARBITS} = $input->{BARBITS};
    } else {
	$fillings->{BARBITS} = ($fillings->{KIND} == $AA) ? 4.3 : 2;
    }

    if (defined $input->{TICBITS}) {
	$fillings->{TICBITS} = $input->{TICBITS};
    } else {
	$fillings->{TICBITS} = 1;
    }




#    if (isChecked($input->{NOCOLOR})) {
#        $fillings->{COLORSCHEME} = "b";
#    } else {
#        $fillings->{COLORSCHEME} = ($kind == $AA) ? "a" : "n";
#    }

    #color
    if (defined $input->{DEFAULT_COLOR}) {
	$fillings->{DEFAULT_COLOR} = (isHexColor( $input->{DEFAULT_COLOR})) ? "c" . $input->{DEFAULT_COLOR} :
	                              $input->{DEFAULT_COLOR};
    }

    if (isChecked($input->{SHOWENDS})) {
	$fillings->{SHOWENDS} = ($fillings->{KIND} == $AA) ? "p" : "d";
    } else {
        $fillings->{SHOWENDS} = "-";
    }

    $fillings->{DESC} = getDescription($desc_r, $fillings->{KIND});

    $fillings->{ERRORBARFRACTION} = $input->{ERRORBARFRACTION};
    $fillings->{COLORSCHEME} = $input->{COLORSCHEME};
    $fillings->{COLORS} = $input->{COLORS};
}

sub getDescription {
    my $returnVal = "";

    foreach (@{$_[0]}) {
        if(defined($_)) {
            $returnVal .= "% * $_\n";
        } else {
            $returnVal .= "% * \n";
        }
    }

    if ($_[1] == $AA) {
        $returnVal .= "% * PROTEIN ALIGNMENT";
    } else {
        $returnVal .= "% * NUCLEOTIDE ALIGNMENT";
    }

    return $returnVal;
}


################################################################################
#####                        SETTING DEFAULTS                              #####
################################################################################

sub setDefaults {
    my ($fillings, $defaults, $numchars) = @_;

    $fillings->{LOGOHEIGHT} = $defaults->{LOGOHEIGHT} if !defined $fillings->{LOGOHEIGHT};
    $fillings->{LOGOWIDTH} = $defaults->{LOGOWIDTH} if !defined $fillings->{LOGOWIDTH};

    $fillings->{START_NUM} = $defaults->{START_NUM} if !defined $fillings->{START_NUM};    
    $fillings->{LOGOSTART} = $fillings->{START_NUM} if !defined $fillings->{LOGOSTART};
    $fillings->{LOGOEND} = $numchars + $fillings->{LOGOSTART} - 1 if !defined $fillings->{LOGOEND};

    $fillings->{YAXIS} = $defaults->{YAXIS} if !defined $fillings->{YAXIS};
    $fillings->{TITLE} = $defaults->{TITLE} if !defined $fillings->{TITLE} || $fillings->{TITLE} eq "";
    #$fillings->{YAXIS_LABEL} = $defaults->{YAXIS_LABEL} if !defined $fillings->{YAXIS_LABEL} || $fillings->{YAXIS_LABEL} eq "";
    $fillings->{YAXIS_LABEL} = $defaults->{YAXIS_LABEL} if !defined $fillings->{YAXIS_LABEL} ;
    $fillings->{XAXIS_LABEL} = $defaults->{XAXIS_LABEL} if !defined $fillings->{XAXIS_LABEL} || $fillings->{XAXIS_LABEL} eq "";

    $fillings->{BARENDS} = $defaults->{BARENDS} if !defined $fillings->{BARENDS};
    $fillings->{OUTLINE} = $defaults->{OUTLINE} if !defined $fillings->{OUTLINE};
    $fillings->{SHOWINGBOX} = $defaults->{SHOWINGBOX} if !defined $fillings->{SHOWINGBOX};
    $fillings->{NUMBERING} = $defaults->{NUMBERING} if !defined $fillings->{NUMBERING};

    $fillings->{ERRORBARFRACTION} = $defaults->{ERRORBARFRACTION} if !defined $fillings->{ERRORBARFRACTION};
    $fillings->{SHRINKFACTOR} = $defaults->{SHRINKFACTOR} if !defined $fillings->{SHRINKFACTOR};
    $fillings->{ERRBAR} = $defaults->{ERRBAR} if !defined $fillings->{ERRBAR};

#    printf("logostart = %d, start num = %d, numchars = $numchars, logoend = %d\n", $fillings->{LOGOSTART}, $fillings->{START_NUM},
#	   $fillings->{LOGOEND});

    my $givenrange = $fillings->{LOGOEND} - $fillings->{LOGOSTART} + 1;
    my $possiblerange = $numchars - ($fillings->{LOGOSTART} - $fillings->{START_NUM});

    if (!defined $fillings->{CHARSPERLINE} && ($givenrange > $possiblerange)) {
	$fillings->{CHARSPERLINE} = $numchars - ($fillings->{LOGOSTART} - $fillings->{START_NUM});
    } elsif (!defined $fillings->{CHARSPERLINE}) {
	$fillings->{CHARSPERLINE} = $fillings->{LOGOEND} - $fillings->{LOGOSTART} + 1;
    }

    $fillings->{DEFAULT_COLOR} = $defaults->{DEFAULT_COLOR} if !defined $fillings->{DEFAULT_COLOR} ||
	                                                       $fillings->{DEFAULT_COLOR} eq "";

#    printf("chars per line = %s\n",$fillings->{CHARSPERLINE});    
#    print("givenrange = $givenrange, possiblerange = $possiblerange\n");

    if ($givenrange > $possiblerange) {
	$fillings->{LOGOLINES} = roundup($possiblerange / $fillings->{CHARSPERLINE});
    } else {
	$fillings->{LOGOLINES} = roundup($givenrange / $fillings->{CHARSPERLINE});
    }

    $fillings->{CHARWIDTH} = ($fillings->{LOGOWIDTH} - 1.5) / $fillings->{CHARSPERLINE};

#    # change height if more than 1 line
#    $fillings->{LOGOHEIGHTPOINTS} = int($fillings->{LOGOHEIGHT} * (72 / 2.54)) * $fillings->{LOGOLINES};

    # LOGOHEIGHTPOITNS is the height input by the user
#    $fillings->{LOGOHEIGHTPOINTS} = int($fillings->{LOGOHEIGHT} * (72 / 2.54));                          # user specifies height of entire logo
    $fillings->{LOGOHEIGHTPOINTS} = int($fillings->{LOGOHEIGHT} * (72 / 2.54)) * $fillings->{LOGOLINES};  # user specifies height of logo line
    $fillings->{LOGOWIDTHPOINTS} =  int($fillings->{LOGOWIDTH}  * (72 / 2.54));

    # LOGOLINEHEIGHT is the height of each logo line, in cm
#    $fillings->{LOGOLINEHEIGHT} = $fillings->{LOGOHEIGHT} / $fillings->{LOGOLINES};    # user specifies height of entire logo
    $fillings->{LOGOLINEHEIGHT} = $fillings->{LOGOHEIGHT};                              # user specifies height of logo line
}

sub roundup {
    return ($_[0] - int($_[0]) > 0) ? int($_[0] + 1) : $_[0];
}


################################################################################
#####                         COLORS                                       #####
################################################################################

sub getDefaultColors {
    my ($defaults) = @_;
    my $returnVal = "";
    $returnVal .= "/black " . $defaults->{black} . " def\n";
    $returnVal .= "/red " . $defaults->{red} . " def\n";
    $returnVal .= "/green " . $defaults->{green} . " def\n";
    $returnVal .= "/blue " . $defaults->{blue} . " def\n";
    $returnVal .= "/yellow " . $defaults->{yellow} . " def\n";
    $returnVal .= "/purple " . $defaults->{purple} . " def\n";
    $returnVal .= "/orange " . $defaults->{orange} . " def\n";

    return $returnVal;
}

sub getNAColors {
    my $returnVal = <<END
% Standard DNA/RNA color scheme
/colorDict << 
   (G)  orange
   (T)  red   
   (C)  blue  
   (A)  green 
   (U)  red   
END
    ;

    return $returnVal;
}
    
sub getAAColors {
    my $returnVal = <<END
% Standard Amino Acid colors
/colorDict << 
  (G)  green  
  (S)  green  
  (T)  green  
  (Y)  green  
  (C)  green  
  (N)  purple 
  (Q)  purple 
  (K)  blue   
  (R)  blue   
  (H)  blue   
  (D)  red    
  (E)  red    
  (P)  black  
  (A)  black  
  (W)  black  
  (F)  black  
  (L)  black  
  (I)  black  
  (M)  black  
  (V)  black
END
    ;

    return $returnVal;
}

sub setColors {
    my ($fillings, $defaults, $input) = @_;
    my $colordef = getDefaultColors($defaults);
    my $colordict = "/colorDict <<\n";

    if ($fillings->{COLORSCHEME} eq "DEFAULT") {
	$colordef = getDefaultColors($defaults);

	if ($fillings->{KIND} eq $AA) {
	    $colordict = getAAColors();
	} else {
	    $colordict = getNAColors();
	}
    } elsif ($fillings->{COLORSCHEME} eq "BW") {
	# do nothing for dict
    } else {

	my %colorhash = %{ $fillings->{COLORS} };
	my $colorName = "";

	foreach (keys %colorhash) {  # keys are strings of residues, value = color name or color code (FF0000)
	    # add color to definitions
	    $colorName = $colorhash{$_};

#	    print("color = $_\n");

	    addColorDef(\$colordef, $colorName ) if (isHexColor($colorName));
	    
	    # add have each residue use the color
	    foreach (split(//, $_)) {
		# add color to dictionary
		if (isHexColor($colorName)) {
		    $colordict .= "  ($_)  c$colorName\n" if !($_ =~ /^\s*$/);
		} else {
		    $colordict .= "  ($_)  $colorName\n" if !($_ =~ /^\s*$/);
		}
	    }
	}
    }

    $colordict .= "\n>> def";

    # add to fillings
    $fillings->{COLORDEF} = $colordef;
    $fillings->{COLORDICT} = $colordict;
}

sub addColorDef {
#    print("adding to color def\n");
    my ($colordef_r, $color) = @_;
    my $PSColor = getPSColor($color);
    $$colordef_r .= "/c$color $PSColor def\n";
}

sub isHexColor {
    return ($_[0] =~ /^[0-9a-fA-F]+$/) && (length $_[0] == 6);
}

# know that it is hex color
sub getPSColor {
    return "[" . hex(substr($_[0],0,2)) / 255 . "  " .
          	 hex(substr($_[0],2,2)) / 255 . "  " .
		 hex(substr($_[0],4,2)) / 255 . "]";
}


################################################################################
#####                         SETTING DATA FIELD                           #####
################################################################################

sub setData {
    my ($fillings, $data_r) = @_;

    my @data = @$data_r;
    my ($height, $letter);
    my @slice;
    my $data;
    my $start_num = $fillings->{START_NUM};

    my $start = $fillings->{LOGOSTART} - $start_num;    # where in @data to start
    my $end = $fillings->{LOGOEND} - $start_num;        # where in @data to end
    my $charsperline = $fillings->{CHARSPERLINE};

    my $numlabel = $fillings->{LOGOSTART};

    $end = ($end >= scalar @data) ? (scalar @data - 1) : $end;

    for (my $i=$start ; $i<=$end ; $i++) {

	# if add new lines
#	if ((($i - $start) % $charsperline == 0) &&
#	    ($i != $start) &&    # not first one
#	    ($i != $end)) {      # not last one
	if ((($i - $start) % $charsperline == 0) &&
	    ($i != $start)) {    # not first one
	    $data .= <<END
EndLine
StartLine

END
    ;
	}

        @slice = @{$data[$i]};
        $data .= <<END
($numlabel) startstack
END
    ;

        $numlabel++;

        foreach (@slice) {
            ($letter,$height) = ($_ =~ /^(.{1})(\S+)/);

	    # is space, so leave
	    if ($letter eq " ") {
		last;
	    }
	    
	    # look for ">", which is symbol for error bar, then quit
	    if ($letter eq ">") {
		last;
	    }

#	    # look for negative heights
#	    if ($height < 0) {
#		next;
#	    }

            $letter =  (uc $letter);  # always uppercase
	    $height = ($height < 0) ? 0 : $height;
            $data .= " $height ($letter) numchar\n";
        }

        # put in error bars -- size is in $height as read in before
        if ($fillings->{ERRBAR} && $letter ne " " && $height != 0) {
	    $data .= " $height Ibeam\n";
	}

        $data .= <<END
endstack

END
    ;

    }

    $fillings->{DATA} = $data;
}

################################################################################

1;
