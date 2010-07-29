#!/usr/bin/perl

# $Id: pipeline_trawler_03.pl,v 1.32 2009/05/05 15:01:04 ramialis Exp $

=pod

=head1 NAME

  pipeline_trawler_03.pl

=head1 SYNOPSIS

    pipeline_trawler_03.pl -directory runName -conservation 1/0

=head1 DESCRIPTION

  A script which parses the output files from pipeline_trawler_01.pl and pipeline_trawler_02.pl for a given run:
    1- runName.cluster (list of over-represented instances and SD)
    2- features/STAT.tmp (conservation count)
    3- runName_compare.txt (TFBS database hits)
    The output is an html file.

=head1 OPTIONS

    -directory #name of the run

=head1 CONTACT

  EMBL 2008
  Mirana Ramialison ramialis@embl.de
  Yannick Haudry haudry@embl.de

=cut

#==============================================================================
# Read and Set properties
#==============================================================================

use strict;
use Carp;
use File::Basename;
use Getopt::Long;
use CGI;
use CGI qw(:standard :html3 *table);
use File::Spec;

# Locate Trawler modules
use FindBin ();
use lib "$FindBin::RealBin/../modules";
my $script_name = $FindBin::RealScript;

# Trawler Modules
use Trawler::Constants 1.0 qw(_read_config _tcst);

# START processing
print "\n## Running $script_name\n";

###############################################################################

use constant HTML_INDEX  =>  'index.html';

use constant CLASS_TAB  =>  'tab_panel';
use constant CLASS_SORTER  =>  'tablesorter';
use constant CLASS_TOOLTIP  =>  'tooltip';

use constant LINK_JS_FILE => 'js/base.js';
use constant LINK_CSS_FILE => 'css/style.css';

use constant LINK_TRAWLER => 'http://ani.embl.de/laurence/blog/';
use constant LINK_TRAWLER_WEB => 'http://ani.embl.de/trawler/';

use constant MAILTO_LAURENCE => 'mailto:ettwille@embl.de';
use constant MAILTO_BENEDICT => 'mailto:benedict@soe.ucsc.edu';
use constant MAILTO_MIRANA => 'mailto:ramialis@embl.de';
use constant MAILTO_YANNICK => 'mailto:haudry@embl.de';

###############################################################################

#==============================================================================
# Read and Set properties
#==============================================================================

# Read config file
_read_config($FindBin::RealBin);
my %tcst = _tcst();

# Logging Levels
my $DEBUG = $tcst{DEBUG};
my $INFO  = $tcst{INFO};

#############################################################
# Set parameters                                            #
#############################################################
#global parameters, to uncheck when used done

#PATH
my $RES_PATH = $tcst{RES_PATH};

#==============================================================================

my $directory = undef;
my $conservation = undef;
GetOptions(
  'directory=s'    => \$directory,
  'conservation=s' => \$conservation,
);

# are we in conservation mode ?
if ($conservation) {
  print "HTML ouput with conservation\n" if $DEBUG;
}

# if no directory provided => exit
unless($directory) {
    print STDERR "\nERROR: no directory to run\n";
    exit(1);
}

# TRAWLER VERSION
my $version = $tcst{trawler_version};

# FIXME[YH]: naming conventions..
my ($RES_DIR_NAME, $RES_DIR_PATH) = fileparse($directory);
my $WORKING_DIR = File::Spec->catdir( $directory, $tcst{RES_DIR_NAME} );


my $COMPARE_FILE_NAME = File::Spec->catfile($WORKING_DIR, $RES_DIR_NAME . "_compare.txt");
my $STAT_FILE_NAME = File::Spec->catfile(($WORKING_DIR, "features"), $tcst{STAT_FILE_NAME} );
my $CLUSTER_FILE_NAME = File::Spec->catfile($WORKING_DIR, $RES_DIR_NAME . $tcst{CULSTER_FILE_EXT} );
my $PWM_FILE = File::Spec->catfile($WORKING_DIR, $RES_DIR_NAME . $tcst{PWM_FILE_EXT} );
my $HTML_INDEX = File::Spec->catfile($directory, HTML_INDEX);

# FIXME[YH]: use property
my $JS_SCRIPT = LINK_JS_FILE;
my $CSS_SCRIPT = LINK_CSS_FILE;

#### Downloadable files [not system dependent !!!, hence do not use cat() function]
my $input_link_dir = $tcst{INPUT_DIR_NAME} . "/";
my $result_link_dir = $tcst{RES_DIR_NAME} . "/";

my $input_file_link = $input_link_dir . $tcst{INPUT_FILE};
my $readme_file_link = $input_link_dir . $tcst{README_FILE};
my $license_file_link = $input_link_dir . $tcst{LICENSE_FILE};
my $trawler_raw_file_link = $result_link_dir . $RES_DIR_NAME . $tcst{TRAWLER_FILE_EXT};
my $trawler_sorted_file_link = $result_link_dir . $RES_DIR_NAME . $tcst{TRAWLER_SHORT_FILE_EXT};
my $cluster_file_link = $result_link_dir . $RES_DIR_NAME . $tcst{CULSTER_FILE_EXT};
my $pwm_file_link = $result_link_dir . $RES_DIR_NAME . $tcst{PWM_FILE_EXT};
####

my $html_download_dir = File::Spec->catdir( $directory, $tcst{HTML_DOWNLOAD} );
my $html_input_dir = File::Spec->catdir( $directory, $tcst{HTML_INPUT} );

if ($DEBUG) {
    print "abs path[RES_DIR_PATH]: $RES_DIR_PATH \n";
    print "dir name[RES_DIR_NAME]: $RES_DIR_NAME \n";
    print "working directory[WORKING_DIR]: $WORKING_DIR \n";
    print "cluster file path[CLUSTER_FILE_NAME]: $CLUSTER_FILE_NAME \n";
    print "compare file path[COMPARE_FILE_NAME]: $COMPARE_FILE_NAME \n";
    print "stat file path[STAT_FILE_NAME]: $STAT_FILE_NAME \n";
    print "pwm file [PWM_FILE]: $PWM_FILE \n";
    print "result path $RES_PATH \n";
    print "$JS_SCRIPT \n";
}

#############################################################
# HTML / JavaScript                                         #
#############################################################

# NOTES
# - external links configuration: rel="external"

# declare scripts, styles
# init Tabs
sub printHTMLhead() {
 return <<Head;
<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
     "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
<head>
  <title>TRAWLER: over-represented motifs</title>
  <meta http-equiv="content-type" content="text/html;charset=utf-8" />
  <meta http-equiv="Content-Style-Type" content="text/css" />
  <meta name="keywords" content="motif discovery" />
  <meta name="copyright" content="copyright 2008 EMBL" />
  <link rel="stylesheet" type="text/css" href="$CSS_SCRIPT" />
  <!--[if IE]><script language="javascript" type="text/javascript" src="js/excanvas.js"></script><![endif]-->
  <script src="$JS_SCRIPT" type="text/javascript"></script>
  <script type="text/javascript">
  // <![CDATA[
    \$(function(){
      \$("a[rel*='external']").click(function(){this.target="_blank";});
      \$("#container").tabs({selected:2,fx:{opacity:"toggle"}});
      sortIndex();
      \$(".delete").click(function(){\$(this).parent("div").hide("slow");return false;});
      \$(".tooltip a").hover(function(){\$(this).next("em").animate({opacity:"show",top:"-80"},"slow");},
                             function(){\$(this).next("em").animate({opacity:"hide",top:"-85"},"fast");});
    });
  // ]]>
  </script>
</head>
<body>

Head
}

#############################################################
# Run Main                                                  #
#############################################################
msg_pipeline(); # console output

#load cluster
print "Loading Cluster\n" if $DEBUG;
my %cluster = loadFileToHash($CLUSTER_FILE_NAME);
print "Cluster loaded from $CLUSTER_FILE_NAME\n" if $INFO;

#load compare
print "Loading Compare\n" if $DEBUG;
my %compare = loadFileToHash($COMPARE_FILE_NAME);
print "Compare loaded from $COMPARE_FILE_NAME\n" if $INFO;

#load stat
print "Loading Stat\n" if $DEBUG;
my %stat = loadFileToHash($STAT_FILE_NAME);

#load stat data for graphic display
my @stat_data = parse_stat($STAT_FILE_NAME);

print "Stat loaded from $STAT_FILE_NAME\n" if $INFO;

#create main list hash
my %families;
my %families_instances;
foreach my $h (keys %cluster){
  if ($h>0) { #avoids first line
    my @line=@{$cluster{$h}};
    my $fname=$line[1];
    my $new_SD = $line[2];
    $families_instances{$fname."_".$line[0]}=$new_SD;
    my $old_SD=0;
    if (exists $families{$fname}) { #gets old SD for this family
      $old_SD=$families{$fname};
    }
    if ($new_SD>$old_SD){ #check whether best SD
      $families{$fname}=$new_SD;
    }
  }
}

my %main_list;
my $tab_count = 2; # 2 is the Result tab [input | download | results | famailies...]
foreach my $family(keys %families) {
  $tab_count++;

  # relative link to image directory
  my $image = File::Spec->catfile($tcst{HTML_IMG}, $family . $tcst{MOTIF_PNG_EXT});

  # tab link (image column)
  my $image_scr = a( { href => '#', onclick => "\$('#container').tabs('select', $tab_count);" },
                       img( { -src => $image, -class => 'small', -alt=>$image } ) );

  # tab link (family column)
  my $html_link = a( {href=>'#', onclick=>"\$('#container').tabs('select', $tab_count);"}, $family );

  #extract compare Table1
  my $compare_list;
  my $compare_family;
  $compare_family = start_table( { class => 'tablesorter tooltip', id => $family."Tf" } );
#	foreach my $hits_num_comp (sort keys %compare){
#        if ($compare{$hits_num_comp}[0] eq 'Query_ID'){ #headline
#            $compare_family.=Tr(th([
#                                    @{$compare{$hits_num_comp}}
#                                    ]));
#        }
#        if ($compare{$hits_num_comp}[0] eq $family){
#            $compare_family.=Tr(td([@{$compare{$hits_num_comp}}]));
#            if (length($compare{$hits_num_comp}[2])>0){
#                $compare_list.=$compare{$hits_num_comp}[2].";";
#            }
#        }
#    }
  my $compare_family_thead;
    my $compare_family_tbody;
    my $Query_ID="Name of the input family matrix (query)";
  my $Query_Consensus="Consensus of the query matrix";
  my $Subject_Name="Name of the hit matrix with divergence smaller than the given cutoff (subject)";
  my $Source_DB="Source database of the subject matrix";
  my $Subject_ID="Database ID of the subject matrix";
  my $Length="Length of the consensus sequence of the subject matrix";
  my $Orientation="Orientation between query and subject matrices ";
  my $Offset="Shift between query and subject matrices";
  my $Divergence="Dissimilarity score between query and subject matrices";
  my $Overlap="Overlap between query and subject matrices";
  my $Subject_Consensus="Consensus sequence of the subject matrix";

  foreach my $hits_num_comp (sort keys %compare){
      # headline
    if ($compare{$hits_num_comp}[0] eq 'Query_ID'){
      #$compare_family_thead .= Tr( th( [ @{$compare{$hits_num_comp}} ] ) );
      $compare_family_thead .= Tr( th( [
                   div( a({href=>"#"}, $compare{$hits_num_comp}[0]), em($Query_ID)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[1]), em($Query_Consensus)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[2]), em($Subject_Name)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[3]), em($Source_DB)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[4]), em($Subject_ID)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[5]), em($Length)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[6]), em($Orientation)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[7]), em($Offset)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[8]), em($Divergence)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[9]), em($Overlap)),
                   div( a({href=>"#"}, $compare{$hits_num_comp}[10]), em($Subject_Consensus)),
                   ] ) );
    }
    # tbody TF db
    if ($compare{$hits_num_comp}[0] eq $family) {
        # Add TRANSFAC / JASPAR links
        if ($compare{$hits_num_comp}[3] eq 'JASPAR') {
           my $id = $compare{$hits_num_comp}[4];
           $compare{$hits_num_comp}[3] = a( {
            href => "http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl?ID=".$id."&rm=present&db=0",
            title => "JASPAR:".$id,
            rel => "external" },
            'JASPAR' );
        } elsif ($compare{$hits_num_comp}[3] eq 'TRANSFAC') {
            my $id = $compare{$hits_num_comp}[4];
            $compare{$hits_num_comp}[3] = a( {
                href => "http://www.biobase.de/cgi-bin/biobase/transfac/8.4/bin/getTFProf.cgi?".$id,
                title => "TRANSFAC:".$id,
                rel => "external" },
                'TRANSFAC' );
        }
        # Create table rows
      $compare_family_tbody .= Tr( td( [ @{$compare{$hits_num_comp}} ] ) );
      if (length($compare{$hits_num_comp}[2])>0){
        $compare_list .= $compare{$hits_num_comp}[2].";";
      }
    }
  } # end build TF db Table
  #print "****THEAD\n";
  #print $compare_family_thead."\n";
  #print "****TBODY\n";
  #print $compare_family_tbody."\n";

  $compare_family .= thead( $compare_family_thead );
  $compare_family .= tbody( $compare_family_tbody );
  $compare_family .= end_table();

  #extract sequences
  my %unique_positions;
  my %family_bestCS;
  foreach my $hits_num_seq (sort keys %stat) {
    if ($stat{$hits_num_seq}[1] eq $family) { # treats one family at a time
      #gets family best CS
      my $new_CS = $stat{$hits_num_seq}[5];
      my $old_CS = -1;
      if (exists $family_bestCS{$family}) { #gets old_CS for this family
         $old_CS = $family_bestCS{$family};
      }
      if ($new_CS>$old_CS) { #check whether $new_CS
         $family_bestCS{$family} = $new_CS;
      }

      #checks if same position
      my $instance_position = $stat{$hits_num_seq}[0]."_".$stat{$hits_num_seq}[3];
      # FIXME[YH]: use href as id
      #my $jalview_link = createJalview($stat{$hits_num_seq}[0], $family.'_'.$stat{$hits_num_seq}[0]);
      my $jalview_link = $stat{$hits_num_seq}[0]; # seq_name

      if (exists $unique_positions{$instance_position}) {

        my $bestCS_at_this_position = $unique_positions{$instance_position}[4];
        my $instance_with_bestCS_at_this_position = $unique_positions{$instance_position}[5];
        if ($stat{$hits_num_seq}[5] > $unique_positions{$instance_position}[4]) { #best CS
          $bestCS_at_this_position = $stat{$hits_num_seq}[5];
          $instance_with_bestCS_at_this_position = $stat{$hits_num_seq}[2];
        }

        my $bestSD_at_this_position=$unique_positions{$instance_position}[5];
        my $instance_with_bestSD_at_this_position = $unique_positions{$instance_position}[6];
        if ($families_instances{$family."_".$stat{$hits_num_seq}[2]} > $unique_positions{$instance_position}[5]) { #best SD
          $bestSD_at_this_position = $families_instances{$family."_".$stat{$hits_num_seq}[2]};
          $instance_with_bestSD_at_this_position = $stat{$hits_num_seq}[2];
        }

        if ($conservation) {
          @{$unique_positions{$instance_position}} = ( a( { -href => "#", title => $jalview_link }, $stat{$hits_num_seq}[0] ),
                                                       $stat{$hits_num_seq}[3],
                                                       $stat{$hits_num_seq}[8], # start position from end
                                                       $stat{$hits_num_seq}[4],
                                                       $bestCS_at_this_position,
                                                       $instance_with_bestCS_at_this_position,
                                                       $bestSD_at_this_position,
                                                       $instance_with_bestSD_at_this_position,
                                                       $stat{$hits_num_seq}[7]); # strand
        }
        else {
          @{$unique_positions{$instance_position}} = ( a( { -href => "#", title => $jalview_link }, $stat{$hits_num_seq}[0] ),
                                                       $stat{$hits_num_seq}[3],
                                                       $stat{$hits_num_seq}[8], # start position from end
                                                       $bestSD_at_this_position,
                                                       $instance_with_bestSD_at_this_position,
                                                       $stat{$hits_num_seq}[7]); # strand

        }

      # if ($stat{$hits_num_seq}[4]>$unique_positions{$instance_position}[2]){ #best CS
      #	@{$unique_positions{$instance_position}}=(a({-href=>$jalview_link},$stat{$hits_num_seq}[0]),$stat{$hits_num_seq}[3],$stat{$hits_num_seq}[4],$stat{$hits_num_seq}[5],$stat{$hits_num_seq}[2],$unique_positions{$instance_position}[4],$unique_positions{$instance_position}[5]);
      # }
      # if ($families_instances{$family."_".$stat{$hits_num_seq}[2]}>$unique_positions{$instance_position}[5]){ #best SD
      #	@{$unique_positions{$instance_position}}=(a({-href=>$jalview_link},$stat{$hits_num_seq}[0]),$stat{$hits_num_seq}[3],$stat{$hits_num_seq}[4],$stat{$hits_num_seq}[5],$stat{$hits_num_seq}[2],$families_instances{$family."_".$stat{$hits_num_seq}[2]},$stat{$hits_num_seq}[2]);
      # }

      }
      else { #new position
         if ($conservation) {
           @{$unique_positions{$instance_position}} = ( a( { -href=>"#",title=>$jalview_link }, $stat{$hits_num_seq}[0] ),
                                                        $stat{$hits_num_seq}[3],
                                                        $stat{$hits_num_seq}[8], # start from end
                                                        $stat{$hits_num_seq}[4],
                                                        $stat{$hits_num_seq}[5],
                                                        $stat{$hits_num_seq}[2],
                                                        $families_instances{$family."_".$stat{$hits_num_seq}[2]},
                                                        $stat{$hits_num_seq}[2],
                                                        $stat{$hits_num_seq}[7]); # strand

         }
         else {
           @{$unique_positions{$instance_position}} = ( a( { -href=>"#",title=>$jalview_link }, $stat{$hits_num_seq}[0] ),
                                                        $stat{$hits_num_seq}[3],
                                                        $stat{$hits_num_seq}[8], # start position from end
                                                        $families_instances{$family."_".$stat{$hits_num_seq}[2]},
                                                        $stat{$hits_num_seq}[2],
                                                        $stat{$hits_num_seq}[7]); # strand


         }
      }
    }
  }
  my $seq_family = start_table( { class => 'tablesorter tooltip', id => $family."Occ" } );
  # Tooltip text (on table header)
  my $seq_txt='Click on the sequence name to view the motif in the sequence (Jalview)';
  my $sta_txt='Start position of the motif within the sequence (from start)';
  my $end_txt='Start position of the motif within the sequence (from end)';
  my $avc_txt='Average conservation of the whole sequence';
  my $bcs_txt='Best conservation score of the motif within this sequence';
  my $ibc_txt='Consensus sequence of the instance of the motif with best conservation score';
  my $bzs_txt='Best Z-score of the motif within this sequence';
  my $ibz_txt='Consensus sequence of the instance of the motif with best Z-score';
  my $strand_txt='Strand: 1=forward, -1=reverse';
  if ($conservation) { # with alignments
    $seq_family.= thead( Tr( th [
            #'Sequence', 'Start_position (from start)', 'Start_position (from end)', 'Average_conservation',
            #'Best_conservation_score', 'Instance_with_best_CS', 'Best_Z-score', 'Instance_with_best_ZS'
            div( a({href=>"#"}, "Sequence"), em($seq_txt) ),
            div( a({href=>"#"}, "Start_position (from start)"), em($sta_txt) ),
            div( a({href=>"#"}, "Start_position (from end)"), em($end_txt) ),
            div( a({href=>"#"}, "Average_conservation"), em($avc_txt) ),
            div( a({href=>"#"}, "Best_conservation_score"), em($bcs_txt) ),
            div( a({href=>"#"}, "Instance_with_best_CS"), em($ibc_txt) ),
            div( a({href=>"#"}, "Best_Z-score"), em($bzs_txt) ),
            div( a({href=>"#"}, "Instance_with_best_ZS"), em($ibz_txt) ),
            div( a({href=>"#"}, "Strand"), em($strand_txt) )
            ]
          ) );
  }
  else { # without alignments
    $seq_family.= thead( Tr( th [
            #'Sequence', 'Start_position (from start)', 'Start_position (from end)',
            #'Best_Z-score', 'Instance_with_best_ZS'
            div( a({href=>"#"}, "Sequence"), em($seq_txt) ),
            div( a({href=>"#"}, "Start_position (from start)"), em($sta_txt) ),
            div( a({href=>"#"}, "Start_position (from end)"), em($end_txt) ),
            div( a({href=>"#"}, "Best_Z-score"), em($bzs_txt) ),
            div( a({href=>"#"}, "Instance_with_best_ZS"), em($ibz_txt) ),
            div( a({href=>"#"}, "Strand"), em($strand_txt) )
            ]
          ) );
  }

  my $seq_family_tbody;
  for my $kf (keys %unique_positions) {
      $seq_family_tbody .= Tr( td( [ @{$unique_positions{$kf}} ] ) );
  }
  $seq_family .= tbody( $seq_family_tbody );
  $seq_family .= end_table();

  createFamilyHTML($family, $compare_family, $seq_family);
  #populate the main list
  # FIXME[YH]: tablink
    #push( @{$main_list{$family}},a({href=>$html_link},$family));
    push( @{$main_list{$family}}, $html_link );
    push( @{$main_list{$family}}, $image_scr );
    push( @{$main_list{$family}}, $compare_list );
    push( @{$main_list{$family}}, $families{$family} );
    if ($conservation) {
       push( @{$main_list{$family}}, $family_bestCS{$family});
    }
}
print "Main list loaded\n" if $DEBUG;

################################################################################
# Create HTML index

open(FINDEX, '>'.$HTML_INDEX) or croak "Can not create index.html file ($HTML_INDEX): $!";
print FINDEX printHTMLhead(); # header and scripts
#### --- HTML header
print FINDEX "<div id=\"tr-siteContain\">";
print FINDEX "<div id=\"tr-header\">";
print FINDEX div( {id=>'tr-helpNav'},
             ul(
             li( a( { href => '#', onclick => 'showDiv(\'h_links\');' }, "Links"), ' | ' ),
             li( a( { href => '#', onclick => 'showDiv(\'h_readme\');' }, "Readme"), ' | ' ),
             li( a( { href => '#', onclick => 'showDiv(\'h_license\');' }, "License"), ' | ' ),
             li( a( { href => '#', onclick => 'showDiv(\'h_contact\');' }, "Contact") ),
             ) );
print FINDEX "</div>"; # tr-header
# TITLE
print FINDEX div( h2("Over-represented motifs discovery results") );
# HELP containers: links
print FINDEX div( { id => 'h_links', class => 'help hide-first' },
                  a( { href => '#', class => 'delete' }, "Hide" ),
                  a( { href => LINK_TRAWLER, rel => 'external' }, "Trawler Home Page"), br(),
                  a( { href => LINK_TRAWLER_WEB, rel => 'external' }, "Trawler Web"), br(),
                  );
print FINDEX div( { id=>'h_contact', class=>'help hide-first' },
                  a( { href => '#', class=>'delete' }, "Hide"),
                  "Contacts: ", br(),
                  a( { href => MAILTO_LAURENCE }, "Laurence Ettwiller"), br(),
                  a( { href => MAILTO_BENEDICT }, "Benedict Paten"), br(),
                  a( { href => MAILTO_MIRANA }, "Mirana Ramialison"), br(),
                  a( { href => MAILTO_YANNICK }, "Yannick Haudry"), br()
                );
print FINDEX div( { id=>'h_license', class => 'help hide-first' },
                  a( { href => $license_file_link,
                       rel => 'external',
                       onmouseover => "showfile('$license_file_link', 'lfile');" },
                       "LICENSE" ),
                  a( { href => '#', class => 'delete'}, "Hide" ),
                  pre( { id => 'lfile' } )
                );
print FINDEX div( { id=>'h_readme', class=>'help hide-first' },
                  a( { href => $readme_file_link,
                       rel => 'external',
                       onmouseover => "showfile('$readme_file_link', 'rfile');" },
                       "README" ),
                  a( { href => '#', class => 'delete'}, "Hide" ),
                  pre( { id => 'rfile' } )
                );

####
print FINDEX "<div id=\"container\">";

#### --- tabbed menu ---
my $family_tab_index;
foreach my $family(keys %families) {
  $family_tab_index .= li( a( { href => $family.".html", title => $family.'_tab' }, span($family) ) );
}
print FINDEX ul(
               li( a( { href => '#uinput' }, span('Input') ) ),
               li( a( { href => '#udownload' }, span('Download') ) ),
               li( a( { href => '#trmain' }, span('Results') ) ),
               $family_tab_index
);

#### --- Input tab ---
my $input_tab = div( { id => 'uinput', class => CLASS_TAB },
                b("Input:"), br(),
                ul(
                  li( a( { href => $input_file_link,
                           rel => 'external',
                           onmouseover => "showfile('$input_file_link', 'ifile');"},
                           "Trawler options")
                     ) # li input file
                  ), # ul
                  # div container input file
                  div( { class=>'pane hide-first' },
                          a( { href=>'#', class=>'delete'}, "Hide" ),
                          pre( {id=>'ifile'} )
                     )
                ); # uinput div
$input_tab .= input_file_script();
# print input tab
print FINDEX $input_tab;

#### --- Download tab ---

#Legend: Trawler raw data and trawler sorted
my $legend_trawler = "<b>column legend:</b><br/>" .
                     "1- Motif occurrence in the sample <br/>" .
                     "2- Motif occurrence in the background <br/>" .
                     "3- Z score <br/>" .
                     "4- Motif";
#Legend: Clustered motif
my $legend_cluster = "<b>column legend:</b><br/>" .
                     "1- Motif <br/>" .
                     "2- Family name <br/>" .
                     "3- Z score <br/>" .
                     "4- Occurrence in the sample <br/>" .
                     "5- Occurrence in the background <br/>" .
                     "6- Strand";

my $download_tab = div( { id => 'udownload', class => CLASS_TAB },
                   b("Download:"), br(),
                   ul( # downloadable files
                     li( a( { href => $trawler_raw_file_link,
                             rel => 'external',
                              onmouseover => "showfilelegend('$trawler_raw_file_link', 'dlfile', '$legend_trawler');" },
                              "Trawler raw data" ) ),
                     li( a( { href => $trawler_sorted_file_link,
                             rel => 'external',
                              onmouseover => "showfilelegend('$trawler_sorted_file_link', 'dlfile', '$legend_trawler');" },
                              "Trawler sorted data" ) ),
                     li( a( { href => $cluster_file_link,
                             rel => 'external',
                              onmouseover => "showfilelegend('$cluster_file_link', 'dlfile', '$legend_cluster');" },
                              "clustered motifs" ) ),
                     li( a( { href => $pwm_file_link,
                             rel => 'external',
                              onmouseover => "showfile('$pwm_file_link', 'dlfile');" },
                              "PWMs" ) ),
                   ),
                   # download file container
                   div( { class => 'pane hide-first' },
                     a( { href => '#', class => 'delete' }, "Hide" ),
                     pre( { id => 'dlfile' } )
                   )
); # end download tab HTML
print FINDEX $download_tab;

##### --- Results tab ---
# Page title
my $main_tab_title = h1('TRAWLER\'S OVER-REPRESENTED MOTIFS');
# Tooltip text (on table header)
my $fam_txt = 'Motif family';
my $pwm_txt = 'Position Weight Matrix';
my $hit_txt = 'From TRANSFAC (Knueppel et al., 1994), JASPAR (Sandelin et al., 2004), homeobox matrices (Noyes et al., 2008; Wilson et al., 2008)';
my $zsc_txt = 'Over-representation score of this motif compared to its background occurrence';
my $csc_txt = 'Number of species where the motif is conserved';
# print THEAD
my $main_tab_table_thead;
if ($conservation) {
    $main_tab_table_thead = thead( Tr( th [
          div( a( { href => "#" }, "Motif"), em($fam_txt) ),
          div( a( { href => "#" }, "PWM"), em($pwm_txt) ),
          div( a( { href => "#" }, "Hits against known TFBS databases"), em($hit_txt) ),
          div( a( { href => "#" }, "Z-score"), em($zsc_txt) ),
          div( a( { href => "#" }, "Conservation_score"), em($csc_txt) )
          #'Motif', 'PWM', 'hits against known TFBS databases', 'Z-score', 'Conservation_score'
          ] ) );
}
else {
    $main_tab_table_thead = thead( Tr(th [
          div( a( { href => "#" }, "Motif"), em($fam_txt) ),
          div( a( { href => "#" }, "PWM"), em($pwm_txt) ),
          div( a( { href => "#" }, "Hits against known TFBS databases"), em($hit_txt) ),
          div( a( { href => "#" }, "Z-score"), em($zsc_txt) )
          #'Motif', 'PWM', 'hits against known TFBS databases', 'Z-score'
          ] ) );
}
# print TBODY
my $main_tab_table_tbody;
for my $k (keys %main_list) {
    $main_tab_table_tbody .= Tr( td( [ @{$main_list{$k}} ] ) );
}
# print main DIV
print FINDEX div( { id => 'trmain', class => CLASS_TAB},
             $main_tab_title,
             start_table( { id => 'idx', class => ''.CLASS_SORTER.' '.CLASS_TOOLTIP.'' } ),
             $main_tab_table_thead,
             $main_tab_table_tbody,
             end_table() );


#### --- Family tabs (Ajax mode) ---
my $family_tabs;
# NOTE: div ID must match title attribute (container)
foreach my $family(keys %families) {
  $family_tabs .= div( { id => $family.'_tab', class => CLASS_TAB }, '' );
}
print FINDEX $family_tabs;

#### --- Finalize tab container ---
print FINDEX "</div>"; # container div

#### --- HTML Footer ---
print FINDEX div( { id => 'footer' }, "Trawler " . $version);

#### --- Finalize page ---
print FINDEX "</div>"; # site-contain
print FINDEX end_html;

#### --- Close index file ---
close(FINDEX) or croak "Can't close file '$HTML_INDEX': $!";
print "html index created\n" if $DEBUG;


#############################################################
# Sub routines                                              #
#############################################################

sub loadFileToHash {
  #------------------------------------------------------------------
  #General script to load the lines of the file into a hash of arrays
  #------------------------------------------------------------------
  my ($file_to_load) = @_;

  my %returned_hash;
  open(F, $file_to_load) or croak "Can not open file ($file_to_load): $!";

  my $hitno = 0;
  while (my $ligne = <F>) {
    $hitno++;
    chomp $ligne;
    my @c = split/\t/, $ligne;
    @{$returned_hash{$hitno}} = @c;
  }
  close(F) or croak "Can't close file '$file_to_load': $!";
  return %returned_hash;

} # end loadFileToHash()

# $family.Occ $family.Tf $family.form
sub print_family_script {
 my ($family_name) = @_;
 my $tid = $family_name . "Occ";

 return <<Script;
  <script type="text/javascript">
  // <![CDATA[
    \$("a[rel*='external']").click(function(){this.target="_blank";});
    tablefunc("$family_name");
    \$("#$tid td a").click(function(){runJalview("$family_name", this.title, "jalview_$family_name");});
    \$(".tooltip a").hover(function(){\$(this).next("em").animate({opacity:"show",top:"-80"},"slow");},
                           function(){\$(this).next("em").animate({opacity:"hide",top:"-85"},"fast");});
  // ]]>
  </script>
Script
}

sub input_file_script {
 return <<Script;
  <script type="text/javascript">
  // <![CDATA[
    showfile('$input_file_link', 'ifile');
  // ]]>
  </script>
Script
}

sub createFamilyHTML {
  my ($family_name, $compare_family_name, $seq_family_name) = @_;

  #create Family HTML
  my $FAMHTML_INDEX = File::Spec->catfile($directory, $family_name.".html");
  open(FAMINDEX, '>'.$FAMHTML_INDEX) or croak "Can not create index.html file ($FAMHTML_INDEX): $!";
  # FIXME[YH]: remove header declaration
  #print FAMINDEX printHTMLhead();

  # Jalview container
  print FAMINDEX div( { id => 'jalview_'.$family_name } );

  #header
  # FIXME[YH]: relative link to images (images folder)
   #my $image_link = $WORKING_DIR."/".$family_name."_all_motif.png";
   my $image_link = File::Spec->catfile($tcst{HTML_IMG}, $family_name . $tcst{MOTIF_PNG_EXT});

   my $PWMFAM = createPWMfile($family_name); # create PWM file for this family
   my $CLUSTFAM = createClusterfile($family_name); # create the motifs file for this family
   # create HTML link for these two files
   my $PWMlink = File::Spec->catfile($tcst{HTML_DOWNLOAD}, $PWMFAM);
   my $CLUSTERlink = File::Spec->catfile($tcst{HTML_DOWNLOAD}, $CLUSTFAM);
   print FAMINDEX table( { class => 'pwm' },
          thead(
          #Tr( [
          Tr( th( h2($family_name) ) ),
            #th([h2($family_name)]),
            #th([img({-src=>$image_link,-alt=>$image_link})]),
            #th([a({-href=>$PWMlink},'Download PWM')])
          ), # end thead
          tbody( Tr( td( img( { -src => $image_link, -alt => $image_link } ) ) ), # PWM logo
                 Tr( td( a( { href => $PWMlink, rel => 'external' }, 'Download PWM' ) ) ), # PWM file
                 Tr( td( a( { href => $CLUSTERlink, rel => 'external' }, 'Download instances (motifs)' ) ) ), # motifs file
                 Tr( td( a( { href => "#", onclick => "showPlot();" }, 'Motifs distribution' ) ) ), # distribution plot
          ), # end tbody
        );

  # motifs distribution plot container / script
  my @fam_graph = family_stat_graph($family_name, @stat_data);
  print FAMINDEX create_graph_javascript($family_name, @fam_graph)."<br/>";

  #database hits Table
  print FAMINDEX h1('Hits against transcription factor databases');
  print FAMINDEX "<br>".$compare_family_name."<br/>";

  #sequences and conservation Table
  print FAMINDEX h1('Occurrences of the motif in the input sequences').br();

  print FAMINDEX create_filter_form($family_name)."<br/>";
  print FAMINDEX $seq_family_name."<br/>";

  #print FAMINDEX "</ul><br />";

  print FAMINDEX print_family_script($family_name);

  # FIXME[YH]: remove html structure
  # print FAMINDEX end_html;
  close(FAMINDEX) or croak "Can't close file '$FAMHTML_INDEX': $!";

  print "html $family_name created\n" if $DEBUG;

}

sub create_graph_javascript {
   my ($fname, @fam_graph) = @_;

#<a href="#" onclick="return showPlot();">Motifs distribution</a>
   my $data_gr = join ",", @fam_graph;

  return <<GraphJS;
<div id="plot-$fname" style="width:600px;height:300px;display:none;"></div><br />
<div id="plot-$fname-legend" style="display:none;">Y-axis: number of motifs<br />X-axis: % of sequence length</div>
<script language="javascript" type="text/javascript">
function showPlot() {
    var fdata = [$data_gr];
    \$.plot(\$("#plot-$fname"),[{data:fdata,bars:{show: true}}],{xaxis:{tickFormatter:function(v,axis){return v * 10}}});
    \$("#plot-$fname").show();
    \$("#plot-$fname-legend").show();
}
</script>
GraphJS
}

sub create_filter_form {
   my ($family_name) = @_;

  return <<FilterForm;
<form id="$family_name-form">Search:
  <input name="$family_name-filter" id="$family_name-filter" value="" maxlength="30" size="30" type="text">
</form><br />
FilterForm
}

# Creates the PWM file
sub createPWMfile {
  my ($PWM_name) = @_;

  open(PWM, $PWM_FILE) or croak "Can not open .pwm file ($PWM_FILE): $!";

  # Retrieve pwm files from download directory
  my $PWM_FAMFILE_NAME = $PWM_name."_pwm.txt";
  my $PWM_FAMFILE = File::Spec->catfile($html_download_dir, $PWM_FAMFILE_NAME);

  open(FAMPWM, '>'.$PWM_FAMFILE) or croak "Can not create family.pwm file ($PWM_FAMFILE): $!";
  my $started='false';
  while (my $ligne = <PWM>) {
    if ($ligne=~/^>/){
      $started='false';
    }
    if ($ligne=~/^>$PWM_name/) {
      $started='true';
    }
    if ($started eq 'true') {
      print FAMPWM $ligne;
    }
  }
  close(PWM) or croak "Can't close file '$PWM_FILE': $!";
  close(FAMPWM) or croak "Can't close file '$PWM_FAMFILE': $!";

  return $PWM_FAMFILE_NAME;
}

# Creates the clustered motifs file
sub createClusterfile {
  my ($Clutser_name) = @_;

  open(CLUSTER, $CLUSTER_FILE_NAME) or croak "Can not open .pwm file ($CLUSTER_FILE_NAME): $!";

  # Retrieve cluster file from download directory
  my $CLUSTER_FAMFILE_NAME = $Clutser_name."_cluster.txt";
  my $CLUSTER_FAMFILE = File::Spec->catfile($html_download_dir, $CLUSTER_FAMFILE_NAME);
  open(FAMCLUSTER, '>'.$CLUSTER_FAMFILE) or croak "Can not create family.pwm file ($CLUSTER_FAMFILE): $!";
  print FAMCLUSTER "#Motif\tFamily name\tZ score\tOccurrence in the sample\tOccurrence in the background\tStrand\n";
  while (my $ligne = <CLUSTER>) {
    if ($ligne=~/$Clutser_name/) {
      print FAMCLUSTER $ligne;
    }
  }
  close(CLUSTER) or croak "Can't close file '$CLUSTER_FILE_NAME': $!";
  close(FAMCLUSTER) or croak "Can't close file '$CLUSTER_FAMFILE': $!";

  return $CLUSTER_FAMFILE_NAME;
}

sub parse_stat {

    my ($stat_file) = @_;

    open(FILE, $stat_file) or croak "Can not open STAT file ($stat_file): $!";
    my @extract;
    my $line = <FILE>; # skip header line
    while ($line = <FILE>) {
      # 0: id | 1: family | 3: start | 6: seq length
        my @array = split(/\t/,$line);
        push @extract, "$array[1];$array[0];$array[3];$array[6]";
    }
    close(FILE) or croak "Can't close STAT file '$stat_file': $!";
    # uniq and sorted
    my %seen = map { $_, 1 } @extract;
    my @uniqed = sort keys %seen;

    return @uniqed;
}

sub family_stat_graph {

    my ($family_id, @stat_array) = @_;

    my @scoring;
    foreach my $elem ( @stat_array ) {
        my @uniq_extract = split /;/, $elem;
        if ($uniq_extract[0] eq $family_id) {
            my $calc = sprintf("%.2f", ($uniq_extract[2]/$uniq_extract[3])*10);
            push @scoring, $calc;
        }
    }

    my %matrix = ();
    # generate beans
    foreach my $elem (@scoring) {
        $elem = int($elem); # absc
        if (exists $matrix{$elem}) {
          my $tmp_val = $matrix{$elem};
            $matrix{$elem} = $tmp_val + 1;
        }
        else {
          $matrix{$elem} = 1;
        }
    }

    # check all beans
    foreach(0..9) {
        if (!exists $matrix{$_}) {
            $matrix{$_} = 0;
        }
    }

    # generate data graph for javascript
    my @data_graph;
    for my $key ( sort keys %matrix ) {
        my $value = $matrix{$key};
        push @data_graph, "[$key,$value]";
    }

    return @data_graph;
}

sub msg_pipeline {
  my $msg_pipeline = <<"MSG";
  ==========
  Creating HTML output...
  ==========
MSG

  print $msg_pipeline . "\n";
}

1;
