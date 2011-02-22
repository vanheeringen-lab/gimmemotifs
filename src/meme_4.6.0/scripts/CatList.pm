# File: CatList.pm
# Project: website
# Author: James Johnson
# Created: June 2010
#
# Description:
#   This module loads the list of lists used in our database selection system. 
#   It uses a binary format created from the csv files and has a file size 
#   limitation of 4Gb (due to use of 32 bit integers to store file offsets).
#   
#
# Format Notes: 
#     1) It's a binary format - there are no spaces, gaps, or anything 
#        other than the fields in the <> brackets.
#     2) Unless otherwise specified to be utf-8 then everything is a 
#        bigendian unsigned 32-bit integer.
#     3) Strings (utf-8) are not null terminated but instead have their 
#        length in the preceeding 4 bytes.
#
#
# Format:
#   ---------- version -------------
#   <version number>
#   ----------- counts -------------
#   <total num categories> <total num options>
#   -------- category index --------
#   <category 1 entries offset> 
#   <category 2 entries offset>
#   ... more category index items ...
#   <category n entries offset>
#   --------- entry index ----------
#   <entry fields 1 offset>
#   <entry fields 2 offset>
#   ... more entry index items ...
#   <entry fields n offset>
#   -------- category names --------
#   <category 1 name length> <category 1 name (utf-8)>
#   <category 2 name length> <category 2 name (utf-8)>
#   ... more categories ...
#   <category n name length> <category n name (utf-8)>
#   ------- entry fields ---------
#   <category 1 entries index> <category 1 entries count>
#   <category 1 name length> <category 1 name (utf-8)>
#   <entry 1 fields count> <entry 1 field 1 length> <entry 1 field 1 value (utf-8)> ... more ... <entry 1 field x length> <entry 1 field x value (utf-8)>
#   <entry 2 fields count> <entry 2 field 1 length> <entry 2 field 1 value (utf-8)> ... more ... <entry 2 field y length> <entry 2 field y value (utf-8)>
#   <category 2 entries index> <category 2 entries count>
#   <category 2 name length> <category 2 name (utf-8)>
#   ... more entry and category fields ...
#   <entry n fields count> <entry n field 1 length> <entry n field 1 value (utf-8)> ... more ... <entry n field z length> <entry n field z value (utf-8)>
#   

package CatList;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = ('load_categories', 'load_entry', 'open_entries_iterator', 'index_catlist_csv',
  'DB_ROOT', 'DB_PROT', 'DB_NUCL', 'DB_SHORT', 'DB_NAME', 'DB_DESC', 
  'F_PROT_URL', 'F_NUCL_URL', 'G_PROM_START', 'G_PROM_STOP', 'G_COMP_ROOT', 'G_COMP_DESC');
our %EXPORT_TAGS = (index_constants => [qw(DB_ROOT DB_PROT DB_NUCL DB_SHORT DB_NAME DB_DESC)], 
  fasta_index_constants => [qw(F_PROT_URL F_NUCL_URL)], gomo_index_constants => [qw(G_PROM_START G_PROM_STOP G_COMP_ROOT G_COMP_DESC)]);

use Fcntl qw(:DEFAULT SEEK_CUR SEEK_SET);
use Encode;

use constant VERSION => 1;
use constant FPOS_SIZE => 4;
use constant INT_SIZE => 4;
use constant CINDEX_SIZE => FPOS_SIZE;
use constant EINDEX_SIZE => FPOS_SIZE;

#indexes common to all db types
use constant DB_ROOT => 0; # used for generating the names of the files
use constant DB_PROT => 1; # is there a amino-acid form avaliable (yes/no/<empty>)
use constant DB_NUCL => 2; # is there a nucleic-acid form avaliable (yes/no/<empty>)
use constant DB_SHORT => 3; # is the nucleotide sequence short (yes/no/<empty>)
use constant DB_NAME => 4; # the name shown in menus
use constant DB_DESC => 5; # the description

#indexes for fasta db columns
use constant F_PROT_URL => 6; # the source url for protein sequences?
use constant F_NUCL_URL => 7; # the source url for nucleotide sequences?

#indexes for gomo db columns
use constant G_PROM_START => 6; # the promoter start
use constant G_PROM_STOP => 7; # the promoter stop
use constant G_COMP_ROOT => 8; # the base used to generate the files for a comparative database
use constant G_COMP_DESC => 9; # the descriptive name for the comparative database


##############################################################################
# PUBLIC FUNCTIONS
##############################################################################

#
# load_categories (
#   $database
#   );
#
#   Given a database filename load the category names.
#
#   Returns:
#     1) The names of the categories.
#
sub load_categories {
  my $db_filename = shift;
  die("Missing db file.") unless ((defined $db_filename) && (-e $db_filename));
  
  my $db_fp;
  sysopen($db_fp, $db_filename, O_RDONLY) or die("sysopen failed!");
  binmode($db_fp) or die("binmode failed!");

  my ($total_cats, $total_entries) = read_header($db_fp);

  seek($db_fp, ($total_cats * CINDEX_SIZE + $total_entries * EINDEX_SIZE), SEEK_CUR);

  my @categories = ();
  for (my $i = 0; $i < $total_cats; $i++) {
    my $category = read_string($db_fp);
    push(@categories, $category);
  }
  close($db_fp) or die("close failed!");
  return @categories;
}


#
# load_entry (
#   $database,
#   $entry_index
#   );
#
#   Given the database filename and the index of an entry returns the fields of that entry.
#
#   Returns
#     1) An array containing the fields of the requested entry.
#
sub load_entry {
  my ($db_filename, $entry_index) = @_;
  die("Missing db file.") unless ((defined $db_filename) && (-e $db_filename));
  die("Entry index must be a number larger than or equal to zero.") unless ((defined $entry_index) && ($entry_index =~ m/^\d+$/));

  my $db_fp;
  sysopen($db_fp, $db_filename, O_RDONLY) or die("sysopen failed!");
  binmode($db_fp) or die("binmode failed!");

  my ($total_cats, $total_entries) = read_header($db_fp);
  
  die("Entry index must be a number smaller than the number of entries.") unless ($entry_index < $total_entries);
  
  #lookup the index for the entry
  seek($db_fp, ($total_cats * CINDEX_SIZE) + (EINDEX_SIZE * $entry_index), SEEK_CUR);
  
  #read the file position of the entry
  my $entry_pos = read_fpos($db_fp);

  #seek to the file position of the entry
  seek($db_fp, $entry_pos, SEEK_SET);

  my @entry = read_entry($db_fp);

  close($db_fp) or die("close failed!");

  return @entry;
}

#
# open entries
#
# Create a new entry iterator.
#
sub open_entries_iterator {
  my $db_filename = shift;
  die("Missing db file.") unless ((defined $db_filename) && (-e $db_filename));
  my $category_index = shift;
  die("Category index must either be left undefined or be a number larger than or equal to zero.") 
      unless ((not defined $category_index) || ($category_index =~ m/^\d+$/));
  my $db_fp;
  sysopen($db_fp, $db_filename, O_RDONLY) or die("sysopen failed!");
  binmode($db_fp) or die("binmode failed!");

  my ($total_categories, $total_entries) = read_header($db_fp);

  my $self;
  if ((not defined $category_index) || ($category_index < $total_categories)) {

    if (defined $category_index) {
      seek($db_fp, CINDEX_SIZE * $category_index, SEEK_CUR);
    }

    if ($total_categories > 0) {
      #read the file position for the start of the entries
      my $entries_fpos = read_fpos($db_fp);
      #seek to the start of the entries
      seek($db_fp, $entries_fpos, SEEK_SET);
    }
  } else {
    warn("Category index $category_index is not smaller than the total number of categories $total_categories.");
  }
  $self = {
    _database_file => $db_fp,
    _remaining_category_count => (defined $category_index ? ($category_index < $total_categories ? 1 : 0) : $total_categories),
    _category_entry_count => 0,
    _category_name => undef,
    _current_entry => undef,
    _entry_index => undef,
    _category_seen => -1
  };

  bless $self;
  return $self;
}

#
# load_next
#
# Load the next entry in the database
#
# Returns true on success
#
sub load_next {
  my $self = shift;
  my $remaining = $self->{_remaining_category_count};
  if ($self->{_category_seen} < ($self->{_category_entry_count} - 1)) {
    #reading existing category
  } elsif ($remaining > 0) {
    #read next category
    my ($entry_index_start, $entry_count, $category_name);
    do {
      ($entry_index_start, $entry_count, $category_name)= read_category_header($self->{_database_file});
      $remaining--;
    } while ($entry_count == 0 && $remaining > 0); #handle empty categories
    return 0 if ($entry_count == 0);#return if there were no non-empty categories
    $self->{_remaining_category_count} = $remaining;
    $self->{_category_entry_count} = $entry_count;
    $self->{_entry_index} = $entry_index_start - 1;
    $self->{_category_name} = $category_name;
    $self->{_category_seen} = -1;
  } else {
    #nothing else
    #automaticlly close
    $self->close_entries_iterator();
    return 0;
  }
  #read the next entry
  my @entry = read_entry($self->{_database_file});
  $self->{_current_entry} = \@entry;
  $self->{_entry_index} = $self->{_entry_index} + 1;
  $self->{_category_seen} = $self->{_category_seen} + 1;
  return 1;
}

#
# get_seen
#
# Get the number of times an entry of this category has previous been seen. Zero if this is the first entry of the category.
#
sub get_seen {
  my $self = shift;
  return $self->{_category_seen};
}

#
# get_category
#
# Get the category of the current entry.
#
sub get_category {
  my $self = shift;
  return $self->{_category_name};
}

#
# get_entry
#
# Get the fields of the current entry
#
sub get_entry {
  my $self = shift;
  return @{$self->{_current_entry}};
}

#
# get_index
#
# Get the index of the current entry
#
sub get_index {
  my $self = shift;
  return $self->{_entry_index};
}

#
# close_entries_iterator
#
# Close the file associated with the entries iterator.
#
sub close_entries_iterator {
  my $self = shift;
  #close but ignore error because the iterator can autoclose
  close($self->{_database_file});
}

#
# index_catlist_csv(
#   $csv,
#   $database
#   );
#
#   Given an input csv filename and the target database filename creates a binary 
#   file that has the content of the csv file but with indexes.
#
#   Returns true
#
sub index_catlist_csv {
  my ($csv_filename, $db_filename) = @_;
  die("Missing csv file.") unless ((defined $csv_filename) && (-e $csv_filename));
  die("Missing db file.") unless (defined $db_filename);

  #read the csv file into memory

  my ($csv_fp, $line);
  sysopen($csv_fp, $csv_filename, O_RDONLY) or die("sysopen failed!");

  my @categories = ();
  my %lookup = ();
  my @options = ();
  my $total_opt_count = 0;
  my $prev_cat;
  FILE_LOOP: 
  while ($line = <$csv_fp>) {
    #skip comment lines and empty lines
    next if ($line =~ m/^#/ || $line !~ m/\S/);
    chomp $line;
    my @fields = split(/,/, $line);
    #skip if we don't have the minimum number of fields
    if (scalar(@fields) < DB_NAME) {
      warn("Skipping line because it had inadequate field count: \n". join(",", @fields));
      next FILE_LOOP;
    }
    if (scalar(@fields) > 50) {
      warn("A Line had a massive field count which normally indicates a file format error.");
    }
    # check if this is a category
    if ($fields[DB_ROOT] =~ /^\s*$/ and $fields[DB_NAME] =~ /\S/) {
      #check that all fields but DB_NAME are empty as expected
      for (my $i = 0; $i < scalar(@fields); $i++) {
        next if ($i == DB_NAME);
        if ($fields[DB_ROOT] !~ /^\s*$/) {
          warn("Skipping category because all fields except DB_NAME should be empty: \n". join(",", @fields));
          next FILE_LOOP;
        }
      }
      if ($prev_cat) {
        #for a previous category store the options we saved
        my @options_copy = @options;
        $lookup{$prev_cat} = \@options_copy;
      }
      my $category = $fields[DB_NAME];
      $category =~ s/^-+//g; #remove dashes 
      $category =~ s/-+$//g;
      #add the category
      push(@categories, $category);
      #clear the options list
      @options = ();
      #update the previous category
      $prev_cat = $category;
    } else {
      #check that the fields we know how to parse are correct, 
      #collapse 'no' to '0' and 'yes' to '1' so they can be directly used as perl booleans
      if ($fields[DB_ROOT] !~ m/\S/) {
        #not right, so skip
        warn("Skipping record because DB_ROOT is empty: \n". join(",", @fields));
        next;
      }
      for (my $i = DB_PROT; $i <= DB_SHORT; $i++) {
        if ($fields[$i] =~ m/^\s*(yes|1)\s*$/) { # match yes or 1
          $fields[$i] = '1';
        } elsif ($fields[$i] =~ m/^\s*((no|0)\s*)?$/) { #match no, 0 or empty
          $fields[$i] = '0';
        } else {
          #not right, so skip
          warn("Skipping record because expected (y/1/n/0/<empty>): \n". join(",", @fields));
          next FILE_LOOP;
        }
      }
      if ($fields[DB_NAME] !~ m/\S/) {
        #not right, so skip
        warn("Skipping record because DB_NAME is empty: \n". join(",", @fields));
        next;
      }

      #store the fields
      push(@options, \@fields);
      $total_opt_count++;
    }
  }
  if ($prev_cat) {
    my @options_copy = @options;
    $lookup{$prev_cat} = \@options_copy;
  }
  close($csv_fp) or die("close failed!");

  # DEBUG check that I haven't messed up (too badly) by checking the loaded options
  my $test_count = 0;
  foreach (@categories)  {
    my @cat_opts = @{$lookup{$_}};
    $test_count += scalar(@cat_opts);
  }
  die("Test count got: $test_count but real count was: $total_opt_count") unless ($test_count == $total_opt_count);

  #write the categories out in a quick to parse format

  my $db_fp;
  sysopen($db_fp, $db_filename, O_WRONLY | O_CREAT) or die("sysopen failed!");
  binmode($db_fp) or die("binmode failed!");

  #write the version and the two counts
  syswrite($db_fp, pack('NNN', VERSION, scalar(@categories), $total_opt_count)) or die("syswrite failed!");

  #categories index into options
  my $opt_offset = 0;
  for (my $i = 0; $i < scalar(@categories); $i++)  {
    #as we don't know the actual pointer to the entries yet we write a temporary value
    syswrite($db_fp, pack('N', 0)) or die("syswrite failed!");
  }

  #for each of the entries, write out the pointer to the fields for that entry
  foreach (@categories)  {
    my @cat_opts = @{$lookup{$_}};
    for (my $i = 0; $i < scalar(@cat_opts); $i++) {
      #as we don't know the actual pointer to the entry yet we write a temporary value
      syswrite($db_fp, pack('N', 0)) or die("syswrite failed!");
    }
  }

  #write out the categories
  foreach (@categories) {
    syswrite_string($db_fp, $_);
  }

  #write out the entries
  my @cat_pt = ();
  my @opt_pt = ();
  my $entries_index = 0;
  foreach (@categories) {
    my @cat_opts = @{$lookup{$_}};
    push(@cat_pt, systell($db_fp));
    #write out index of the entries and the number of entries for this category
    syswrite($db_fp, pack('NN', $entries_index, scalar(@cat_opts))) or die("syswrite failed!");
    #write out the category name
    syswrite_string($db_fp, $_);
    #write out the entries for this category
    for (my $i = 0; $i < scalar(@cat_opts); $i++) {
      push(@opt_pt, systell($db_fp));
      my @fields = @{$cat_opts[$i]};
      syswrite($db_fp, pack('N', scalar(@fields))) or die("syswrite failed!");
      my $field;
      foreach $field (@fields) {
        syswrite_string($db_fp, $field);
      }
    }
    #increment the entries index for the next category
    $entries_index += scalar(@cat_opts);
  }

  #go back to the start of the file and fix the pointers
  sysseek($db_fp, 3*INT_SIZE, SEEK_SET) or die("sysseek failed");

  #for each of the categories, write out the pointer to the entries
  foreach (@cat_pt)  {
    syswrite($db_fp, pack('N', $_)) or die("syswrite failed!");
  }

  #for each of the entries, write out the pointer to the fields for that entry
  foreach (@opt_pt) {
    syswrite($db_fp, pack('N', $_)) or die("syswrite failed!");
  }

  close($db_fp) or die("close failed!");
  return 1;
}

##############################################################################
# PRIVATE FUNCTIONS
##############################################################################

#
# Internal Only Function
#
# systell(
#   $file_handle_ref
#   );
#
#   Given a file returns the current position in the file.
#
#   Since I'm using syswrite then I must use sysseek
#   to get the file position instead of tell.
#
#   Returns:
#     1) The current file position or undefined on failure.
#
sub systell {
  return sysseek($_[0], 0, SEEK_CUR) or die("sysseek (for systell implementation) failed!");
}

#
# Internal Only Function
#
# syswrite_string(
#   $file_handle_ref,
#   $string
#   );
#
#   Given a file and a string writes the string to the file.
#
#   As I couldn't figure out the docs on write I used syswrite.
#   This converts the string to utf 8 and writes out its length
#   in bigendian 32-bit int followed by the utf-8 octlets.
#
#   Returns nothing, though I should change this to check for error states.
#
sub syswrite_string {
  my $octlets = encode_utf8($_[1]);
  my $data = pack('Na*', length($octlets), $octlets);
  syswrite($_[0], $data) or die("syswrite failed!");
}

#
# Internal Only Function
#
# read_int (
#   $file_handle_ref
#   );
#
#   Given a file, read an unsigned 4-byte bigendian integer.
#
#   Returns:
#     1) the integer
#
sub read_int {
  my $data;
  my $read_len = read($_[0], $data, INT_SIZE);
  die("read failed!") if ($read_len != INT_SIZE);
  my $pos = unpack('N', $data);
  return $pos;
}

#
# Internal Only Function
#
# read_string(
#   $file_handle_ref
#   );
#
#   Given a file, read a string.
#
#   Reads a string written using syswrite_string.
#   A string consists of:
#   <string length (4 byte bigendian integer)><utf-8 encoded string>
#
#   Returns:
#     1) The string.
#
sub read_string {
  my $len = read_int(@_);
  if ($len > 0) {
    my $octlets;
    my $read_len = read($_[0], $octlets, $len, 0);
    die("read failed!") if ($read_len != $len);
    my $str = decode_utf8($octlets);
    return $str;
  }
  return '';
}

#
# Internal Only Function
#
# read_header (
#   $file_handle_ref
#   );
#
#   Given a file read the header. Check that the version is ok.
#   The header consists of (all are 4 byte bigendian integers):
#   <version><count of categories><count of entries>
#
#   Returns:
#     1) count of categories
#     2) count of entries
#
sub read_header {
  my $version = read_int(@_);
  die("Database has version $version but this program parses version " . VERSION . ".") unless ($version == VERSION);
  my $total_cats = read_int(@_);
  my $total_entries = read_int(@_);
  return ($total_cats, $total_entries);
}

#
# Internal Only Function
#
# read_category_header (
#   $file_handle_ref
#   );
#
#   Given a file read a category header.
#
#   Returns:
#     1) first entry index
#     2) entry count
#     3) category name
#
sub read_category_header {
  my $entry_index = read_int(@_);
  my $entry_count = read_int(@_);
  my $category = read_string(@_);
  return ($entry_index, $entry_count, $category);
}

#
# Internal Only Function
#
# read_entry (
#   $file_handle_ref
#   );
#
#   Given a file read an entry.
#
#   Returns:
#     1) the list of fields
#
sub read_entry {
  my $field_count = read_int(@_);
  my @fields = ();
  for (my $i = 0; $i < $field_count; $i++) {
    push(@fields, read_string(@_));
  }
  return @fields;
}

#
# Internal Only Function
#
# read_fpos (
#   $file_handle_ref
#   );
#
#   Given a file, read a number that is large enough to fit a file position offset.
#   The file position offset is a unsigned 4-byte bigendian integer.
#
#   Returns:
#     1) the read file offset
#
sub read_fpos {
  return read_int(@_);
}
