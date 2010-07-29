#!/usr/bin/perl -w

$|++;   # turn off buffering for perl pre-5.6 because of forks

use vars qw($PATH);

BEGIN {
    open(PATH, "pwd|");
    $PATH = <PATH>;
    chomp $PATH;
    close(PATH);
    unshift(@INC, $PATH, "$PATH/lib");
}

# logo.cgi prototype.
use strict;
use POSIX;
use File::Basename;
#use CGI qw(:standard);
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use logo;
use template;


$CGI::POST_MAX=1024 * 1024;  # max 1 megabyte of data

my $MAX_FILESIZE = 10000;

my %error_msg = (
		 BARBITS => "Error: The maximum bits must be a positive real number, but is: ",
		 BADLENGTH => "Error: Number of characters in each logo line is not consistent, starting at: ",
		 BADFORMAT => "Error: Invalid input format does not conform to FASTA, CLUSTAL, or Flat. ",
		 CHARSPERLINE => "Error: Characters per line must be a positive integer, but is: ",
		 FILESIZE => "Error: File size exceeds $MAX_FILESIZE",
		 FIRSTNUM => "Error: First position number must be an integer, but is: ",
		 INVALIDCOLOR => "Error: Invalid RGB color: ",
		 LOGOENDINT => "Error: Ending range must be an integer, but is: ",
		 LOGOENDBOUNDS => "Error: Ending range must be equal to or greater than beginning range, but is: ",
		 LOGOHEIGHT => "Error: Logo height must be a positive real number <= 100 cm, but is (after conversion): ",
		 LOGOSTARTINT => "Error: Beginning range must be an integer, but is: ",
		 LOGOSTARTBOUNDS => "Error: Beginning range must be greater than or equal to first position number, but is: ",
		 LOGOWIDTH => "Error: Logo width must be a positive real number <= 100 cm, but is (after conversion): ",
		 MULTIPLECOLOR => "Error: Residue was found multiple times in the custom color scheme: ",
		 NODATA => "Error: No Data. Enter sequence data in the text field below and hit the Create Logo button",
		 RES => "Error: Resolution must be a positive integer < 1200, but is: ",
		 SHRINK => "Error: Shrink factor must be a positive real number <= 1, but is: ",
		 );

my $query = new CGI;
my $AA = 0;
my $NA = 1;

#FIXME: Error handling: Bad query data
#FIXME: Error handling: Too much form data

my $command = $query->param('command');

clean_cache();


if (!defined $command) {
    # Initial visit. Return the default html page.
    send_html($query);

} elsif ($command eq "Logo") {
    # Example data passed in from the examples page
    # display the standard html page
    send_html($query);

} elsif ($command eq "Edit Logo") {
    # Example data passed in from the examples page
    # display the standard html page
    send_html($query);

} elsif ($command eq "Create Logo") {
    # Create a logo already
    # If anything goes wrong, control is passed to send_html()
    create_logo($query, \%error_msg);

} else {
    # Unkown command. Huh?
    # Display an error message back to the user
    send_html($query, "Internal Form error: Unkown command");
}



# logo.xml is the template html page. Within that
# file the variables that are listed as ${KEYWORD}
# are replaced by form elements.
sub send_html {
    my ($q, $error_msg) = @_;

    print $q->header;

    my %tokens = ();

    if ( defined $error_msg ) {
        $tokens{"ERROR_MESSAGE"} =
            "<p class=\"error\">" . $error_msg ."</p>";
    }

    $tokens{"FORM_TAG"} = 
        $q->start_multipart_form();

    # Bug fix. If the logo form is accessed with an extra slash,
    # http://weblogo.berkeley.edu//logo.cgi        
    # then the action of the form tag is "//logo.cgi". This 
    # causes an error. Fix by replacing double slashes with single slashes
    $tokens{"FORM_TAG"} =~ s/\/\//\//go ;


    $tokens{"SEQUENCE_FILENAME"} = 
        $q->filefield(-name=>'aligned_file');

#<input type="file" name="file" size="40">

    $tokens{"SEQUENCE_TEXTAREA"} = 
        $q->textarea(-name=>'sequence',
                     -default=>'',
                     -rows=>10,
                     -columns=>80);

    $tokens{"FORMAT_SELECT"} = 
        $q->popup_menu(-name=>'format',
                       -values=>['GIF','PNG','EPS','PDF'],  #,'TEST'],
                       -labels=>{'GIF'=>'GIF (bitmap)',
                                 'PNG'=>'PNG (bitmap)',
                                 'EPS'=>'EPS (vector)',
                                 'PDF'=>'PDF (vector)'   #,
#                                 'TEST'=> 'TEST'
				 },
                       -default=>'PNG');

    $tokens{"LOGOWIDTH_TEXTFIELD"} =
        $q->textfield(-name=>'logowidth',
                      -default=>'18',
                      -size=>6,
                      -maxlength=>80);

    $tokens{"LOGOHEIGHT_TEXTFIELD"} =
        $q->textfield(-name=>'logoheight',
                      -default=>'5',
                      -size=>6,
                      -maxlength=>80);

    $tokens{"LOGOUNITS_SELECT"} = 
        $q->popup_menu(-name=>'logounits',
                       -values=>['cm','inches','pixels','points'],
                       -default=>'cm');



    $tokens{"YAXIS_CHECKBOX"} =
        $q->checkbox(-name=>'yaxis',
                     -checked=>'checked',
                     -label=>'');

    $tokens{"YAXIS_LABEL_TEXTFIELD"} =
        $q->textfield(-name=>'yaxis_label',
                      -default=>'bits',
                      -size=>16,
                      -maxlength=>32);


    $tokens{"BARBITS_TEXTFIELD"} =
        $q->textfield(-name=>'barbits',
                      -default=>'',
                      -size=>4,
                      -maxlength=>80);

    $tokens{"TICBITS_TEXTFIELD"} =
        $q->textfield(-name=>'ticbits',
                      -default=>'1',
                      -size=>4,
                      -maxlength=>80);

    $tokens{"XAXIS_CHECKBOX"} =
        $q->checkbox(-name=>'xaxis',
                     -checked=>'checked',
                     -label=>'');


    $tokens{"XAXIS_LABEL_TEXTFIELD"} =
        $q->textfield(-name=>'xaxis_label',
                      -default=>'',
                      -size=>16,
                      -maxlength=>32);



    $tokens{"FIRSTNUM_TEXTFIELD"} =
        $q->textfield(-name=>'firstnum',
                      -default=>'1',
                      -size=>4,
                      -maxlength=>80);


    $tokens{"LOGOSTART_TEXTFIELD"} =
        $q->textfield(-name=>'logostart',
                      -default=>'',
                      -size=>4,
                      -maxlength=>80);

    $tokens{"LOGOEND_TEXTFIELD"} =
        $q->textfield(-name=>'logoend',
                      -default=>'',
                      -size=>4,
                      -maxlength=>80);




    $tokens{"TITLE_TEXTFIELD"} =
        $q->textfield(-name=>'title',
                      -default=>'',
                      -size=>20,
                      -maxlength=>80);




    $tokens{"RES_TEXTFIELD"} = 
        $q->textfield(-name=>'res',
                      -default=>'96',
                      -size=>4,
                      -maxlength=>8);

    $tokens{"RES_UNITS_SELECT"} = 
        $q->popup_menu(-name=>'res_units',
                       -values=>['ppc','ppi','ppp'],
                       -labels=>{'ppc'=>'pixels/cm',
                                 'ppi'=>'pixels/inch (dpi)',
                                 'ppp'=>'pixels/point'},
                       -default=>'ppi');

    $tokens{"ANTIALIAS"} =
        $q->checkbox(-name=>'antialias',
                     -checked=>'checked',
                     -label=>'');

    $tokens{"FINEPRINT_CHECKBOX"} =
	$q->checkbox(-name=>'fineprint',
                     -checked=>'checked',
		     -label=>'');



    $tokens{"KIND_RADIO"} = 
        $q->radio_group(-name=>'kind',
		  -values=>['AA','NA','AUTO'],
                  -labels=>{'AA'=>' amino acid',
                            'NA'=>' DNA / RNA',
                            'AUTO'=>' Automatic Detection'},
                  -default=>'AUTO');


    $tokens{"SHOWENDS_CHECKBOX"} =
        $q->checkbox(-name=>'showends',
                     -checked=>'checked',
                     -label=>'');



    $tokens{"OUTLINE_CHECKBOX"} =
        $q->checkbox(-name=>'outline',
                     -label=>'');


    $tokens{"BOX_CHECKBOX"} =
        $q->checkbox(-name=>'box',
                     -label=>'');

    $tokens{"SHRINK_TEXTFIELD"} =
        $q->textfield(-name=>'shrink',
                      -default=>'0.5',
                      -size=>4,
                      -maxlength=>80);


    
    $tokens{"SSC_CHECKBOX"} =
        $q->checkbox(-name=>'smallsamplecorrection',
                     -checked=>'checked',
                     -label=>'');

    $tokens{"ERRBAR_CHECKBOX"} =
        $q->checkbox(-name=>'errbar',
                     -label=>'');

    $tokens{"STRETCH_CHECKBOX"} =
	$q->checkbox(-name=>'stretch',
		     -label=>'');


    $tokens{"COLOR_RADIO"} = 
        $q->radio_group(-name=>'colorscheme',
                  -values=>['DEFAULT','BW','CUSTOM'],
                  -labels=>{'DEFAULT'=>' Default',
                            'BW'=>' Black & White',
                            'CUSTOM'=>' Custom'},
                  -default=>'DEFAULT');

    addSymbolTextFields($q, \%tokens);
    addColorOptions($q, \%tokens);
    addRGBTextFields($q, \%tokens);
    
    $tokens{"MULTILINE_CHECKBOX"} =
        $q->checkbox(-name=>'multiline',
                     -label=>'');

    $tokens{"SYMBOLS_TEXTFIELD"} =
        $q->textfield(-name=>'symbolsperline',
                      -default=>'32',
                      -size=>4,
                      -maxlength=>4);


    $tokens{"DEFAULTS_BUTTON"} =
        $q->defaults('Reset');

    print template("logo.xml", \%tokens);

    exit(0);
}

sub addSymbolTextFields {
    my ($q, $tokens_r) = @_;
    my @default_symbols = ("GSTYCQN","KRH","DE","AVLIPWFM","","","");

    for (my $i=1 ; $i<=7 ; $i++) {
        $tokens_r->{"SYMBOL" . $i . "_TEXTFIELD"} = 
            $q->textfield(-name=> "symbol" . $i,
                          -default=> $default_symbols[$i],
                          -size=>16,
                          -maxlength=>32);
    }
}

sub addColorOptions {
    my ($q, $tokens_r) = @_;
    my @default_colors = ("black","green","blue","red","black","purple","orange","black");
    
    for (my $i=0 ; $i<=7 ; $i++) {
        $tokens_r->{"COLOR" . $i . "_SELECT"} = 
            $q->popup_menu(-name=> "color" . $i,
                           -values=>['black','red','green','blue','yellow','purple','orange','RGB=>'],
                           -default=> $default_colors[$i]);
    }
}

# RGB0 is the default one
sub addRGBTextFields {
    my ($q, $tokens_r) = @_;
    
    for (my $i=0 ; $i<=7 ; $i++) {
        $tokens_r->{"RGB" . $i . "_TEXTFIELD"} = 
            $q->textfield(-name=> "rgb" . $i,
                          -default=> '',
                          -size=>8,
                          -maxlength=>8);
    }
}

sub template {
    my ($filename, $fillings) = @_;
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

sub getColors {
    my ($q, $colorhash, $default_color, $error_msg) = @_;
    my %symbolhash = ();
    my ($symbol, $color, $rgb);
    my ($inputsymbol, $inputcolor, $inputrgb);
    my $error = 0;

    #unless( $inputsymbol) {
    #    return $error;
    #}
       
    # go through color fields, build up $colorhash
    for (my $i=0 ; $i<=7 ; $i++) {
        $symbol = "symbol" . $i;
        $color = "color" . $i;
        $rgb = "rgb" . $i;

        $inputsymbol = $q->param($symbol);
        if($inputsymbol) {
        
        $inputcolor = $q->param($color);
        $inputrgb = $q->param($rgb);

        #remove leading, trailing spaces
        $inputrgb =~ s/^\s*//;
        $inputrgb =~ s/\s*$//;

        if ($inputcolor eq "RGB=>") {
            if (!template::isHexColor($inputrgb)) {
                send_html($q, $error_msg->{INVALIDCOLOR} . $inputrgb);
                $error = 1;
                return $error;
            }
            $colorhash->{$inputsymbol} = $inputrgb;
        } else {
            $colorhash->{$inputsymbol} = $inputcolor;
        }

        #check for repeats
        foreach (split(//, $inputsymbol)) {
            if (defined $symbolhash{uc $_}) {
                send_html($q, $error_msg->{MULTIPLECOLOR} . $_);
                $error = 1;
                return $error;
            } else {
                $symbolhash{uc $_} = uc $_;
            }
        }
        }
    }

    # add default color
    $colorhash->{""} = $default_color;

    return $error;
}

sub create_logo {
    my ($q, $error_msg) = @_; 

    my %logo_input = ();

    # Process form input

    # things that can't go wrong
    my $title = $q->param('title');
    my $yaxis_label = $q->param('yaxis_label');
    my $xaxis_label = $q->param('xaxis_label');
    my $outline   = ($q->param('outline')   ) ;
    my $box       = ($q->param('box')       ) ;
    my $xaxis     = ($q->param('xaxis')     ) ;
    my $errbar    = ($q->param('errbar')    ) ;
    my $showends  = ($q->param('showends')  ) ;
    my $yaxis     = ($q->param('yaxis')     ) ;
    my $antialias = ($q->param('antialias') ) ;
    my $fineprint = ($q->param('fineprint') ) ;
    my $stretch   = ($q->param('stretch')   ) ;
    my $input_kind = $q->param('kind');
    if ($input_kind eq "AUTO") {
        $input_kind = undef;
    } elsif ($input_kind eq "AA") {
        $input_kind = $AA;
    } else {
        $input_kind = $NA;
    }
    my $colorscheme = $q->param('colorscheme');
    my $smallsamplecorrection = $q->param('smallsamplecorrection') ;    

# This has changed in the HTML
#    my $charsperline = $q->param('charsperline');
#    if ($charsperline ne "" && !isInt($charsperline)) {
#        send_html($q,"Characters per line must be an integer, but is $charsperline.");
#        return;
#    }    

    my $res = $q->param('res');
    if ($res ne "" && !(isInt($res) && $res > 0 && $res < 1200)) {
        send_html($q, $error_msg->{RES} . $res);
        return;
    }
    my $res_units = $q->param('res_units');    # convert to inches
    if ($res_units eq "ppc") {
        $res *= 2.54;                          # 2.54 cm per inch
    } elsif ($res_units eq "ppp") {
        $res *= 72;                            # 72 points per inch
    }                                          # do nothing for ppi
    my $barbits = $q->param('barbits');
    if ( ($barbits ne "") && !isPosReal($barbits) ) {
        send_html($q,  $error_msg->{BARBITS} . $barbits);
        return;
    }
    my $ticbits = $q->param('ticbits');
    if (($ticbits ne "") && !isPosReal($ticbits) ) {
        send_html($q,  $error_msg->{TICBITS} . $ticbits);
        return;
    }

    my $logowidth = $q->param('logowidth');
    if (!isPosReal($logowidth) || ($logowidth = toCM($logowidth, $res, $q)) > 100) {
        send_html($q, $error_msg->{LOGOWIDTH} . $logowidth);
        return;
    }
    my $logoheight = $q->param('logoheight');
    if (!isPosReal($logoheight) || ($logoheight = toCM($logoheight, $res, $q)) > 100) {
        send_html($q, $error_msg->{LOGOHEIGHT} . $logoheight);
        return;
    }

#    print ("w = $logowidth, h = $logoheight\n");
    my $shrink = $q->param('shrink');
    if (! (isPosReal($shrink) && $shrink <= 1) ) {
        send_html($q, $error_msg->{SHRINK} . $shrink);
        return;
    }
    my $firstnum = $q->param('firstnum');
    if (!isInt($firstnum)) {
        send_html($q, $error_msg->{FIRSTNUM} . $firstnum);
        return;
    }
    my $logostart = $q->param('logostart');
    if ($logostart ne "" && !isInt($logostart)) {
        send_html($q, $error_msg->{LOGOSTARTINT} . $logostart);
        return;
    }
    if ( $firstnum ne "" && $logostart ne "" && $logostart < $firstnum ) {
        send_html($q, $error_msg->{LOGOSTARTBOUNDS} . $logostart);
        return;
    }
    my $logoend = $q->param('logoend');
    if ($logoend ne "" && !isInt($logoend)) {
        send_html($q, $error_msg->{LOGOENDINT} . $logoend);
        return;
    }
    if ($logoend ne "" && $logoend < $logostart) {
        send_html($q, $error_msg->{LOGOENDBOUNDS} . $logoend);
        return;
    }
    my $multiline = $q->param('multiline') ;
    my $charsperline = $q->param('symbolsperline');
    if ($multiline) {
        if ( !isInt($charsperline) || $charsperline <= 0) {
            send_html($q, $error_msg->{CHARSPERLINE} . $charsperline);
            return;
        }
    }

    $logo_input{LOGO_HEIGHT} = $logoheight;
    $logo_input{LOGO_WIDTH} = $logowidth;
    $logo_input{SHOWENDS} = $showends;
    $logo_input{OUTLINE} = $outline;
    $logo_input{NUMBERING} = $xaxis;

    $logo_input{START_NUM} = $firstnum;
    $logo_input{LOGOSTART} = ($logostart eq "") ? undef : $logostart;
    $logo_input{LOGOEND} = ($logoend eq "") ? undef : $logoend;

    $logo_input{YAXIS} = $yaxis;
    $logo_input{STRETCH} = $stretch;    
    $logo_input{TITLETEXT} = $title;
    $logo_input{YAXIS_LABEL} = $yaxis_label;
    $logo_input{XAXIS_LABEL} = $xaxis_label;
    $logo_input{BOXSHRINK} = $shrink;
    $logo_input{SHOWINGBOX} = $box;
    $logo_input{BARBITS} = ($barbits eq "") ? undef : $barbits;
    $logo_input{TICBITS} = ($ticbits eq "") ? undef : $ticbits;
    $logo_input{ERRBAR} = $errbar;
    $logo_input{RES} = $res; #($res eq "" || $res > 1000) ? 50 : $res;
    $logo_input{FORMAT} = $q->param('format');
    $logo_input{ANTIALIAS} = $antialias;
    $logo_input{FINEPRINT} = $fineprint;
    $logo_input{CHARSPERLINE} = ($multiline) ? $charsperline : undef;

    #colors
    $logo_input{COLORSCHEME} = $colorscheme;

    if ($colorscheme eq "CUSTOM") {
        $logo_input{DEFAULT_COLOR} = ($q->param('color0') eq "RGB=>") ? $q->param('rgb0') : $q->param('color0');
    }

    # set color hash
    my %colorhash = ();
    my $error =  getColors($q, \%colorhash, $logo_input{DEFAULT_COLOR}, $error_msg);

    if ($error) {
        return;
    }
    
    $logo_input{COLORS} = \%colorhash;

    my @input;
    if ( (not defined $q->param('aligned_file')) ||
	 ($q->param('aligned_file') =~ /^\s*$/) ) {
#	send_html($q, "reading from box");
	#Error handling: No data
	if( $q->param('sequence') eq '') {
	    send_html($q, $error_msg->{NODATA});
	    return;
	}

	# inputs OK.
	# data ok, get height data now
	@input = split(/\n/, $q->param('sequence'));


    # read file
    } else {
#	my $file = $q->param('aligned_file');	
#	my $a;
#	if ($file =~ /^\s*$/) {
#	    $a = "blank";
#	} else {
#	    $a = "not blank";
#	}
#	my $b = defined $file;
#	send_html($q, "reading from file, " . $b . " " . $a . "erg?");
	my $filehandle = $q->upload('aligned_file');
	if (!$filehandle && $q->cgi_error()) {
	    send_html($q, "Error:" . $q->cgi_error());
	    exit 0;
	}

	{
	    no strict 'vars';
	    my $lines;
	    while (<$filehandle>) {
		push (@input, $_);

		if ($lines++ > $MAX_FILESIZE) {
		    send_html($q, $error_msg->{FILESIZE});
		}
	    }
#	    send_html($q, join("::", @input));
	}
    }

    my %heightparams = (
			smallsampletoggle => $smallsamplecorrection,
			input_kind => $input_kind,
			stretch => $stretch
			);

#    print $q->header;
#    print("<pre>----<br/>");

    my ($heightdata_r, $desc_r, $kind, $goodlength, $badline, $validformat) = 
      logo::getHeightData(\@input, \%heightparams);

    # check for errors
    if ((defined $validformat) && ($validformat == 1)) {
	send_html($q, $error_msg->{BADFORMAT} );
    }
    if (!$goodlength) {
	send_html($q, $error_msg->{BADLENGTH} . $badline);
    }

    #FIXME: Error handling: Bad data
    #FIXME: Error handling: Template failed?

    if ($logo_input{FORMAT} eq "GIF" ||
        $logo_input{FORMAT} eq "PNG" ||
        $logo_input{FORMAT} eq "PDF" ||
        $logo_input{FORMAT} eq "EPS") {

        # make temporary filename
        my $handle ="cache/" . basename( tmpnam() ) 
            . "." . (lc $logo_input{FORMAT});

#	print ("handle = $handle\n");

        # prints automatically
        my $text = template::create_template(\%logo_input, $kind,
                                             $desc_r, $heightdata_r, $handle, $PATH);
        print $q->redirect("$handle");

    # for testing
    } else {
        $logo_input{FORMAT} = "EPS";
        print $q->header;
        print "<pre>";

        my $text = template::create_template(\%logo_input, $kind,
                                             $desc_r, $heightdata_r, "-", $PATH);
    }
}

# negative/positive ints
sub isInt {
    return ($_[0] =~ /^\s*[-\+]?\d+\s*$/) ? 1 : 0;
}

# reals >= 0
sub isPosReal {
    return ($_[0] =~ /^\s*\+?\d*\.?\d*?\s*$/) ? ($_[0] > 0) : 0;
}



sub toCM {
    my ($size, $res, $q) = @_;  # res in pixels/inch
    my $unit = $q->param('logounits');

    if ($unit eq "cm") {
	return $size;
    } elsif ($unit eq "inches") {
	return $size * 2.54;
    } elsif ($unit eq "pixels") {
	return $size / $res * 2.54;    # convert to inches, then to cm
    } else { #elsif ($unit eq "points") {
	return $size / 72 * 2.54;
    }
}


#Remove all files in cache more than 1 day old
sub clean_cache {
  my $file;
  my $dir = "./cache" ;
  opendir (DIR, $dir);
  while(defined($file = readdir DIR)) {
    next if $file =~ /^\.\.?$/; # skip . and ..
    if(-M "$dir/$file" >1) {
      unlink("$dir/$file");
   }
  }
  closedir (DIR);
}





