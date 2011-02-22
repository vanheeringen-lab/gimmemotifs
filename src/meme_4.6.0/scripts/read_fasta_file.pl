################################################################################
#	read_fasta_sequence(fh)
#
#	File handle passed by reference.
#	Returns (name, seq, comment) where comment is anything after the name
#   on the same line; can call ignore the comment value on return. To start
#	a new file, on the first call add a parameter of any defined value at end
################################################################################
sub read_fasta_sequence {
  my ($fh, $newfile) = @_;
  my ($name, $seq, $comment);

  if (defined $newfile) {
  	$__fasta_next_seq_name__ = $__fasta_next_comment__ = undef;
  	$__fasta_eof__ = 0;
  } else {
  	$__fasta_eof__ = 0 unless defined $__fasta_eof__;
  	if ($__fasta_eof__ == 1) { return(()); }	# EOF reached already
  	$name = $__fasta_next_seq_name__;		# use saved name from last call
  	$comment = $__fasta_next_comment__;   # and rest of name line
  }

  $seq = "";

  while (<$fh>) {
    if (/^>(\S+)\s+(.*$)/) {				# found start of next sequence
      $__fasta_next_seq_name__ = $1;		# save next sequence name
      $__fasta_next_comment__ = $2;		# save next sequence name
      if (defined $name) {
	return($name, $seq, $comment);			# return the sequence just read
      } else {
        $name = $__fasta_next_seq_name__;	# first sequence
        $comment = $__fasta_next_comment__; # and rest of name line
      }
    } else {
      s/\s//g;					# remove space from sequence
      s/\.//g;					# remove "."s from sequence
      next if ($_ eq "");
      $seq .= $_;
    }
  }
  
  $__fasta_eof__ = 1;
  return($name, $seq, $comment);				# last sequence
} # read_fasta_sequence

################################################################################
#	read_fasta_file
#
#	Read a FASTA sequence file and return the sequences in an
#	associative array indexed by sequence name.
#
#	$file		open file pointer
#
################################################################################
sub read_fasta_file {
  local ($file) = @_;

  open(F, "<$file") || die("Can't open file $file\n");

  # read the file
  while (<F>) {
    if (/^>(\S+)\s/) {				# found start of sequence
      if ($name) {				# save previous sequence in array
        $seqs{$name} = $seq;
      }
      $name = $1;				# save sequence name
      $seq = "";
    } else {
      s/\s//g;					# remove space from sequence
      s/\.//g;					# remove "."s from sequence
      next if ($_ eq "");
      $seq .= $_;
    }
  } # read file
  if ($name && $seq) {				# save very last sequence in array
    $seqs{$name} = $seq;
  }

  %seqs;					# return array of sequences
} # read_fasta_file

1;
