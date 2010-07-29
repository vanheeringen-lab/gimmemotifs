#!@WHICHPERL@ -w

# Author: Timothy Bailey

# This script copies the list of supported databases to
# the client. It is intended to support clients using 
# XMLHttpRequest to update a selection list of supported databases
# on a web page.
use strict;
use warnings;
use lib qw(@PERLLIBDIR@);
require "Utils.pm";                     # must come before "whines"
use CGI qw/:standard/;
use CGI::Carp qw(fatalsToBrowser set_message);

# Set standard message for error reporting.
set_message('Please report this message to <a href="mailto:@contact@">@contact@</a>.');

my $csv_file = param('db_names');
my $short_only = param('short_only');	# only include databases with "short_seqs"
my $doc = param('doc');			# print documentation (depreciated, use mode=doc)
my $mode = param('mode');               # output mode, doc, xml, list (default)
$mode = "doc" if (!$mode && $doc);      # for backwards compatiblity still support doc
my $categoryIndex = param('catid');     # filter to category at index (xml mode)
my $is_motif_db = param('is_motif_db'); # used by doc mode to provide correct titles

my $dir = "@MEME_DIR@";	# Installation directory
my $cgi = new CGI;

# Open list of supported databases and read them in.
my $result = open(DB_NAMES, "<$dir/etc/$csv_file"); 
if (not $result) {
  print $cgi->header(-status=>'500 Internal Server Error');
  die "Unable to open the list of supported databases.\n";
}
my @DB_NAMES = <DB_NAMES>;
close(DB_NAMES);
my $db_list_empty = (@DB_NAMES == 0);

# create documentation if the databases haven't been installed
if ($db_list_empty) {
  if ($mode eq "doc") {
    $DB_NAMES[0] = ",,,,</a>Click <a href=\"../doc/meme-install.html#$csv_file\">here</a> for database installation instructions.<a>,,,";
  } else {
    $DB_NAMES[0] = ",,,,Click on \"supported database\" for database installation instructions.,,,";
  }
}

# Print desired contents (names or documentation)
if ($mode eq "doc") {
  doc_print();
} elsif ($mode eq "xml") {
  xml_print();
} elsif ($mode eq "categories_xml") {
  cat_xml_print();
} elsif ($mode eq "categories") {
  cat_list_print();
} else {
  list_print();
}

exit(0);

#
# Ouput the list of categories
#
sub cat_list_print {
  my @categories = ();
  for $_ (@DB_NAMES) {
    next if (/^#/ || !/\S/);		# skip comments and blank lines
    chomp;
    my @fields = split(/,/);
    #see if it's a category
    if ($fields[0] !~ /\S/ && $fields[4] =~ /\S/) {
      #trim the category name to remove the unneeded ---
      $fields[4] =~ s/^-+//g; 
      $fields[4] =~ s/-+$//g;
      push(@categories, $fields[4]);
    }
  }
  my $output = join("\n", @categories);

  print "Content-type: text/plain\n\n";
  print $output;
}

#
# Output the xml form of the list of categories
#
sub cat_xml_print {

  print "content-type:text/xml", "\n\n";
  print xml_dtd();
  print "<dblist>\n";
  for $_ (@DB_NAMES) {
    next if (/^#/ || !/\S/);		# skip comments and blank lines
    chomp;
    my @fields = split(/,/);
    #see if it's a category
    if ($fields[0] !~ /\S/ && $fields[4] =~ /\S/) {
      #trim the category name to remove the unneeded ---
      $fields[4] =~ s/^-+//g; 
      $fields[4] =~ s/-+$//g;
      print "<category name=\"" . $fields[4] . "\"/>\n";
    }
  }
  print "</dblist>\n";
}

#
# Output the list of databases
# Fields in db file are:
# <db_base_name>,<prot?>,<nuc?>,<short_seqs?>,<db_menu_name>,<db_long_name>,<prom_start> or <prot_url>,<prom_stop> or <nucl_url>[,comp_db_name,db_long_name]
# 0) The file name of the database
# 1) Is it a protein database
# 2) Is it a nuc database
# 3) Is it a short sequence
# 4) The name that appears in the menu
# 5) The description that appears on the help page
# 6) The start of the sequence relative to the TSS or the protein url
# 7) The stop of the sequence relative to the TSS or the nuc url
# 8) The name of a comparative database.
# 9) The description for the comparative database
#
sub list_print {
  print "Content-type: text/plain\n\n";
  foreach (@DB_NAMES) {
    next if (/^#/ || /^\s*$/);		# skip comments and whitespace
    chomp;
    my @fields = split /,/;
    # remove invisible whitespace from name and yes/no fields
    my $i;
    for ($i=0; $i<@fields; $i++) {
      next if $i >= 4 && $i < 8 || $i > 8 && $i % 2;
      $fields[$i] =~ s/\s+//g;
    }
    $fields[2] = 'no' if ($short_only && $fields[3] eq 'no');	# turn off long DNA if need be
    # create a list of types and add to name if not a header
    if ($fields[0] =~ /\S/ && $fields[4] =~ /\S/) {
      my $db_types;
      if ($is_motif_db) {
        $db_types = "";
      } elsif ($fields[1] eq 'yes' && $fields[2] eq 'yes') {
        $db_types = "(peptide and nucleotide)"; 
      } elsif ($fields[1] eq 'yes') {
        $db_types = "(peptide only)"; 
      } elsif ($fields[2] eq 'yes') {
        $db_types = "(nucleotide only)"; 
      } else {
        $db_types = "(unknown database type)"; 
      }
      $fields[4] .= " $db_types";
    }
    print join(",", @fields) . "\n";
  }
} # list_print 

#
# Output the documentation on supported databases.
# 
sub doc_print {
  # get the number of separator lines
  my @categories = ();
  my @descriptions = ();
  for $_ (@DB_NAMES) {
    next if (/^#/ || !/\S/);		# skip comments and blank lines
    chomp;
    my @fields = split(/,/);
    #next if ($short_only && $fields[3] eq 'no');	# skip categories with only long seqs
    $fields[2] = 'no' if ($short_only && $fields[3] eq 'no');	# turn off long DNA if need be
    if ($fields[0] !~ /\S/ && $fields[4] =~ /\S/) { # check for category
      $fields[4] =~ s/^-+//g; 
      $fields[4] =~ s/-+$//g;
      push(@categories, $fields[4]);    #db menu name
      push(@descriptions, $fields[5]);  #db long name
    }
  } 
  # set globals
  my $db_type1 = "Sequence";
  my $db_type2 = "sequence";
  if ($is_motif_db) {
    $db_type1 = "Motif";
    $db_type2 = "motif";
  }
  my $tab = " ";
  my $form = "";
  my $tabx3 = $tab x 3;
  my $css = 
      $tabx3."h4 {\n".
      $tabx3.$tab."margin-bottom: 0;\n".
      $tabx3."}\n".
      $tabx3."h5 {\n".
      $tabx3.$tab."font-size: 14px;\n".
      $tabx3.$tab."font-family: Arial;\n".
      $tabx3.$tab."font-weight: bold;\n".
      $tabx3."}\n".
      $tabx3.".db {\n".
      $tabx3.$tab."margin-left: 1.5em;\n".
      $tabx3.$tab."margin-top: 0px;\n".
      $tabx3.$tab."margin-bottom: 0px;\n".
      $tabx3."}\n".
      $tabx3."#main ul {\n".
      $tabx3.$tab."margin-top: 0;\n".
      $tabx3.$tab."margin-bottom: 0;\n".
      $tabx3."}\n".
      $tabx3."#main li {\n".
      $tabx3.$tab."font-size: 12px;\n".
      $tabx3."}\n";

  $form .= make_form_header("$db_type1 Databases", "Documentation", $css, 0, $tab);
  $form .= doc_start_body($db_type1, $db_type2, \@categories, \@descriptions, 1, $tab);
  $form .= doc_add_dbs(\@DB_NAMES, 8, $tab);
  $form .= doc_finish_body(1, $tab);
  $form .= "</html>\n";

  print "Content-type: text/html\n\n$form";
} # doc_print

sub doc_start_body {
  my ($db_type1, $db_type2, $categories_p, $descriptions_p, $indent_lvl, $tab) = @_;

  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;

  my $indent = $tab x $indent_lvl;

  # create list of categories
  my ($i, $cat_descr);
  my $ncats = @{$categories_p};
  my $desc_indent = $indent . ($tab x 6);
  for ($i=1; $i<=$ncats; $i++) {
    my $cat = $categories_p->[$i-1];
    my $descr = $descriptions_p->[$i-1];
    $cat_descr .= $desc_indent."<li><h3><a href=#cat$i>$cat</a></h3>\n";
    $cat_descr .= $desc_indent.$descr."\n";
  }
  my $groups_desc = "$ncats categories";
  if ($ncats == 1) {
    $groups_desc = "1 category";
  }

  my $content = 
      $indent."<body class=\"body\">\n".
      $indent.$tab."<table class=\"maintable\">\n".
      $indent.$tab.$tab."<tr>\n".
      $indent.$tab.$tab.$tab."<td class=\"maintablewidth\">\n".
      $indent.$tab.$tab.$tab.$tab."<div id=\"main\">\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<script src=\"../doc/meme-suite-logo.js\" type=\"text/javascript\"></script>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<noscript>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."<h1>MEME Suite</h1>\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."The MEME Suite web application requires the use of JavaScript<br />\n".
      $indent.$tab.$tab.$tab.$tab.$tab.$tab."Javascript doesn't seem to be available on your browser.\n".
      $indent.$tab.$tab.$tab.$tab.$tab."</noscript>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<hr />\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<h1 align=\"center\">$db_type1 databases available for search</h1>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<hr />\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<h3>The $db_type2 databases that can be searched are grouped into $groups_desc:</h3>\n".
      $indent.$tab.$tab.$tab.$tab.$tab."<ul>\n".
      $cat_descr.
      $indent.$tab.$tab.$tab.$tab.$tab."</ul>\n";
} # doc_start_body

sub doc_add_dbs {
  my ($lines, $indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;

  my $i = 1;
  my $content;

  foreach $_ (@$lines) {
    next if (/^#/ || !/\S/);		# skip comments and blank lines
    chomp;
    my @fields = split(/,/);
    $fields[2] = 'no' if ($short_only && $fields[3] eq 'no');	# remove DNA if need be
    my @extra_dbs = ();
    my $ind;
    for ($ind = 9; $ind < @fields; $ind += 2) {
      push(@extra_dbs, $fields[$ind]);
    }
    if ($fields[0] !~ /\S/) {
      $fields[4] =~ s/^-+//g; $fields[4] =~ s/-+$//g;
      $content .= 
          $indent."<hr />\n".
          $indent."<a name=\"cat$i\"></a>\n".
          $indent."<h3>$fields[4]</h3>\n".
          $indent."<hr />\n";
      $i++;
    } else {
      my $db_types;
      if ($is_motif_db) {
        $db_types = "";
      } else {
        if ($fields[1] eq 'yes' && $fields[2] eq 'yes') {
          $db_types = "peptide and nucleotide"; 
        } elsif ($fields[1] eq 'yes') {
          $db_types = "peptide only"; 
        } elsif ($fields[2] eq 'yes') {
          $db_types = "nucleotide only"; 
        } else {
          $db_types = "unknown database type"; 
        }
      }
      $content .= 
          $indent.$tab."<h4>$fields[4]</h4>\n".
          $indent.$tab."<p class=\"db\">$fields[5]</p>\n";
      $content .= $indent.$tab."<p class=\"db\" >($db_types)</p>\n" unless ($is_motif_db); 
      if (@extra_dbs > 0) {
        $content .=
            $indent.$tab."<h5 class=\"db\">Comparative database(s):</h5>\n".  
            $indent.$tab."<ul>\n";
        foreach my $db (@extra_dbs) {
          $content .=
              $indent.$tab.$tab."<li>$db</li>\n";
        }
        $content .=
            $indent.$tab."</ul>\n";
      }
    }
  }
  $content .= 
      $indent."<hr />\n";

  return($content);
} # doc_add_dbs

sub doc_finish_body {
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;

  my $content .= 
      $indent.$tab.$tab.$tab.$tab.$tab."<script src=\"../template-footer.js\" type=\"text/javascript\"></script>\n".
      $indent.$tab.$tab.$tab.$tab."</div>\n".
      $indent.$tab.$tab.$tab."</td>\n".
      $indent.$tab.$tab."</tr>\n".
      $indent.$tab."</table>\n".
      $indent."</body>\n";

  return($content);
} # doc_finish_body


sub xml_print {
  my $catid = 0;
  my $prev = 0;
  print "content-type:text/xml", "\n\n";
  print xml_dtd();
  print "<dblist>\n";
  for $_ (@DB_NAMES) {
    next if (/^#/ || !/\S/);		# skip comments and blank lines
    chomp;
    my @fields = split(/,/);
    #why did the skipping get removed?
    $fields[2] = 'no' if ($short_only && $fields[3] eq 'no');	# turn off long DNA if need be
    #see if it's a category
    if ($fields[0] !~ /\S/ && $fields[4] =~ /\S/) {
      $catid = $catid + 1; #increment the category id
      #optionally filter to a specific category
      next if (defined($categoryIndex) && $categoryIndex != $catid);
      #trim the category name to remove the unneeded ---
      $fields[4] =~ s/^-+//g; 
      $fields[4] =~ s/-+$//g;
      if ($prev) {
        print "</category>\n";
      }
      print "<category name=\"", $fields[4], "\" >\n";
      $prev = 1;
    } else {
      #optionally filter to a specific category
      next if (defined($categoryIndex) && $categoryIndex != $catid);
      print "<db base=\"", fix($fields[0]), "\" prot=\"", yn($fields[1]), "\" nucl=\"", yn($fields[2]);
      print "\" short=\"", yn($fields[3]), "\" menu=\"", fix($fields[4]), "\" desc=\"", fix($fields[5]);
      print "\" ";
      #see if the field looks like the promoter start and end from gomo dbs
      if ($fields[6] =~ m/^[+-]?\d+$/) {
      	print "start=\"", fix($fields[6]), "\" end=\"", fix($fields[7]), "\" ";
      } else {
      	print "prot_url=\"", fix($fields[6]), "\" " if ($fields[6]);
        print "nucl_url=\"", fix($fields[7]), "\" " if ($fields[7]);
      }
      print ">\n";
      for (my $i = 8; ($i + 1) < @fields; $i += 2) {
        print "<cmp_db base=\"", fix($fields[$i]), "\" desc=\"", fix($fields[$i+1]), "\" />\n";
      }
      print "</db>\n";
    }
  } 
  if ($prev) {
    print "</category>\n";
  }
  print "</dblist>\n";
}

sub yn {
  my $test = shift;
  if ($test eq 'yes') {
    return 'y';
  } else {
    return 'n';
  }
}

sub fix {
  my $text = shift;
  $text =~ s/&/&amp;/g;
  $text =~ s/</&lt;/g;
  $text =~ s/>/&gt;/g;
  $text =~ s/"/&quot;/g;
  return $text;
}

sub xml_dtd {
  my ($indent_lvl, $tab) = @_;
  $indent_lvl = 0 unless $indent_lvl;
  $tab = "\t" unless $tab;
  my $indent = $tab x $indent_lvl;
  my $dtd =
      $indent."<?xml version=\"1.0\"?>\n".
      $indent."<!DOCTYPE dblist [\n".
      $indent.$tab."<!ELEMENT dblist (category*)>\n".
      $indent.$tab."<!ELEMENT category (db*)>\n".
      $indent.$tab."<!ATTLIST category name CDATA #REQUIRED>\n".
      $indent.$tab."<!ELEMENT db (cmp_db*)>\n".
      $indent.$tab."<!ATTLIST db base CDATA #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db prot (y|n) #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db nucl (y|n) #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db short (y|n) #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db menu CDATA #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db desc CDATA #REQUIRED>\n".
      $indent.$tab."<!ATTLIST db start CDATA #IMPLIED>\n".
      $indent.$tab."<!ATTLIST db end CDATA #IMPLIED>\n".
      $indent.$tab."<!ATTLIST db prot_url CDATA #IMPLIED>\n".
      $indent.$tab."<!ATTLIST db nucl_url CDATA #IMPLIED>\n".
      $indent.$tab."<!ELEMENT cmp_db (EMPTY)>\n".
      $indent.$tab."<!ATTLIST cmp_db base CDATA #REQUIRED>\n".
      $indent.$tab."<!ATTLIST cmp_db desc CDATA #REQUIRED>\n".
      $indent."]>\n";

  return($dtd);
}
