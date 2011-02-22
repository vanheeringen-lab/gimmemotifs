#!@WHICHPERL@ -w

use strict;
use warnings;

use lib qw(@PERLLIBDIR@);
use CatList qw(load_categories load_entry open_entries_iterator DB_PROT DB_NUCL DB_SHORT DB_NAME DB_DESC G_COMP_DESC);
use CGI;
use CGI::Carp qw(fatalsToBrowser set_message);
use HTML::Template;

# Set standard message for error reporting.
set_message('Please report this message to <a href="mailto:@contact@">@contact@</a>.');

process_request();
exit(0);

#
# process_request
#
# Process the request.
#
# Parameters
# db_names    - one of gomo_db.csv, motif_db.csv or fasta_db.csv
# mode        - one of doc, xml, categories_xml, categories or list
# doc         - for backwards compatibility, enable doc mode
# short_only  - filter the DNA sequences to only allow short versions
# catid       - the category id (for xml mode)
#
sub process_request {
  my $q = CGI->new;

  # Get the name of the database and select the associated file.
  # This used to just take the parameter and append the path but 
  # that had the nasty effect of exposing all files accessable 
  # to the server.
  my $dbname = $q->param('db_names'); #get the database name
  $dbname = "fasta_db.csv" unless (defined $dbname);
  my $index_file = "@MEME_DIR@/etc/";
  my $append_types = 0;
  my $title = "";
  if ($dbname eq "gomo_db.csv") {
    $index_file .= "gomo_db.index";
    $title = "GOMO Databases";
  } elsif ($dbname eq "motif_db.csv") {
    $index_file .= "motif_db.index";
    $title = "Motif Databases";
  } else { # assume fasta_db.csv
    $index_file .= "fasta_db.index";
    $append_types = 1;
    $title = "Sequence Databases";
  }

  # get the output mode
  my $mode = $q->param('mode');                   # output mode: doc, xml, categories_xml, categories, list (default)
  $mode = "list" unless (defined $mode);
  $mode = "doc" if (!$mode && $q->param('doc'));  # for backwards compatiblity still support doc flag

  #read some filtering parameters
  my $short_only = $q->param('short_only');   # filter to only include dna databases with "short_seqs"
  my $category_index = $q->param('catid');     # filter to category at index (xml mode)

  # Print desired contents (names or documentation)
  if ($mode eq "doc") {
    doc_print($q, $index_file, $append_types, $short_only, $title);
  } elsif ($mode eq "xml") {
    xml_print($q, $index_file, $append_types, $short_only, $category_index);
  } elsif ($mode eq "categories_xml") {
    cat_xml_print($q, $index_file);
  } elsif ($mode eq "categories") {
    cat_list_print($q, $index_file); #this may be removed in future as categories_xml does its job
  } else {
    list_print($q, $index_file, $append_types, $short_only); #used by update_db.pl which is possibly on another server
  }
}

#
# cat_list_print
#
# Print out a list of categories with one on each line.
#
sub cat_list_print {
  my ($q, $index_filename) = @_;
  my @categories = load_categories($index_filename);
  print $q->header('text/plain'), join("\n", @categories);
}

#
# list_print
#
# Print out the content of the csv file but with:
# - yes/no converted into 1/0
# - types appended based on $append_types
# - filtered for short DNA sequences based on $short_only
#
sub list_print {
  my ($q, $index_filename, $append_types, $short_only) = @_;
  my $iter = open_entries_iterator($index_filename);
  print $q->header('text/plain');
  while ($iter->load_next()) { #will auto close
    if (!($iter->get_seen())) { # is this the first sighting of the category
      print ',,,,-----', $iter->get_category(), '-----', "\n";
    }
    my @entry = modify_entry([$iter->get_entry()], $append_types, $short_only);
    if (@entry) {
      print join(',', @entry), "\n";
    }
  }
}

#
# doc_print
#
# Print out the documentation. Uses the template db_list_view.tmpl .
#
sub doc_print {
  my ($q, $index_filename, $append_types, $short_only, $title) = @_;
  #open the template
  my $template = HTML::Template->new(filename => 'db_list_view.tmpl', global_vars => 1);
  $template->param(title => $title);
  $template->param(append_types => $append_types);
  
  #fill in parameters
  my $iter = open_entries_iterator($index_filename);
  my @cat_loop_links = ();
  my @cat_loop_main = ();
  my $entries;
  while ($iter->load_next()) { #will auto close
    if (!$iter->get_seen()) {
      my $i = scalar(@cat_loop_main);
      my $category = $iter->get_category();
      push(@cat_loop_links, {i => $i, name => $category});
      $entries = {i => $i, name => $category, db_loop => []};
      push(@cat_loop_main, $entries);
    }
    my @entry = modify_entry([$iter->get_entry()], $append_types, $short_only);
    my @comp_dbs = ();
    for (my $i = G_COMP_DESC; $i < scalar(@entry); $i += 2) {
      push(@comp_dbs, {name => $entry[$i]});
    }
    my $db = {name => $entry[DB_NAME], desc => $entry[DB_DESC], is_nucl => $entry[DB_NUCL], is_prot => $entry[DB_PROT], comp_dbs => \@comp_dbs};
    push(@{$entries->{db_loop}}, $db);
  }
  $template->param(cat_loop_links => \@cat_loop_links);
  $template->param(cat_loop_main => \@cat_loop_main);

  #output the html
  print $q->header;
  $template->output(print_to => *STDOUT);
}

#
# xml_print
#
# Output an xml formatted list of database entries names for a particular category.
#
sub xml_print {
  my ($q, $index_filename, $append_types, $short_only, $cat_id) = @_;
  my $iter = open_entries_iterator($index_filename, $cat_id);
  print $q->header('text/xml'), list_dtd(), "\n<list>\n";
  while ($iter->load_next()) { #will auto close
    my @entry = modify_entry([$iter->get_entry()], $append_types, $short_only);
    if (@entry) {
      my $i = $iter->get_index();
      my $name = $entry[DB_NAME];
      print "  <item index=\"$i\" name=\"$name\"/>\n"
    }
  }
  print "</list>";
}

#
# cat_xml_print
#
# Output an xml formatted list of database categories.
#
sub cat_xml_print {
  my ($q, $index_filename) = @_;
  my @categories = load_categories($index_filename);
  print $q->header('text/xml'), list_dtd(), "\n<list>\n";
  for (my $i = 0; $i < scalar(@categories); $i++) {
    my $name = $categories[$i];
    print "  <item index=\"$i\" name=\"$name\"/>\n"
  }
  print "</list>";
}

#
# list_dtd
#
# Returns the XML DTD for a simple list.
#
sub list_dtd {
  return
      '<?xml version="1.0"?>'."\n".
      '<!DOCTYPE list ['."\n".
      '  <!ELEMENT list (item)>'."\n".
      '  <!ELEMENT item EMPTY>'."\n".
      '  <!ATTLIST item index CDATA #REQUIRED name CDATA #REQUIRED>'."\n".
      ']>';
}

#
# modify_entry
#
# Returns a modified entry with the types appended 
# if $append_types is enabled and filters for short DNA
# sequences if $short_only is enabled.
#
sub modify_entry {
  my ($entry_ref, $append_types, $short_only) = @_;
  my @fields = @{$entry_ref};
  if ($short_only) { # turn off long DNA
    $fields[DB_NUCL] = '0' if (!$fields[DB_SHORT]);
  }
  if ($append_types) {
    #add type to name
    if ($fields[DB_PROT] && $fields[DB_NUCL]) {
      $fields[DB_NAME] .= " (peptide and nucleotide)"; 
    } elsif ($fields[DB_PROT]) {
      $fields[DB_NAME] .= " (peptide only)"; 
    } elsif ($fields[DB_NUCL]) {
      $fields[DB_NAME] .= " (nucleotide only)"; 
    } else {
      return (); #non-valid entry that needs skipping
    }
  }
  return @fields;
}
