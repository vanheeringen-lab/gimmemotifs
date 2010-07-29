#!@WHICHPERL@ -w
# FILE: metatest
# AUTHOR: William Stafford Noble
# CREATE DATE: 7-13-97
# PROJECT: MHMM
# COPYRIGHT: 2008, University of Washington
# DESCRIPTION: Run a complete set of tests on the Meta-MEME software.

my $usage = "
USAGE: metatest [options]

  This script tests the programs ama, fimo, gomo, mhmms, motiph, psp-gen,
  qvalue, shadow, tomtom.

  Options:

    -update :   Turns on update mode, replaces expected output files
                with observed output files when a difference occurs.

    -complete : Test additional scaffolding programs:
                clustalw2fasta, fasta-io, log-hmm, meme-io, mhmm-io

    -continue : Allow the tests to continue after a program test fails.

";


# Set defaults for global variables.
my $complete = 0;
my $continue = 0;
my $update = 0;

# Parse the command line.
while (scalar(@ARGV) > 0) {
  $param = shift(@ARGV);
  if ($param eq "-complete") {
    $complete = 1;
  } elsif ($param eq "-continue") {
    $continue = 1;
  } elsif ($param eq "-update") {
    $update = 1;
  } else {
    print(STDERR "Unrecognized option ($param).\n");
    print(STDERR $usage);
    exit(1);
  }
}
if (scalar(@ARGV) != 0) {
  print(STDERR $usage);
  exit(1);
}

# Check if the env tell us
# where to run (needed by 'make distcheck' in build process)
my $testbindir = `pwd`;
chomp $testbindir;
# directory containing scripts that may need to be tested
$testscriptdir = $testbindir . "/../scripts";
$testbindir = $testbindir . "/../src";
if ($ENV{testdir}) {
	chdir $ENV{testdir};
	print("testdir = $ENV{testdir}\n");
}

if ($complete) {
  &test_scaffold("crp0", 1);
  &test_scaffold("lipo", 0);
}
#&test_core("crp0", 1);
#&test_core("lipo", 0);

# Set the location of the stylesheets
$ENV{'MEME_ETC_DIR'} = '../etc';
# Test mhmms with --maxseqs option (r1635)
&test_program($update, $continue, "mhmms", 
  "--maxseqs 3 --quiet"
    . " mhmm/crp0.linear.mhmm"
    . " common/crp0.fasta",
    "mhmms/r1635.maxseqs.mhmms");

# Test mhmms with --global option (r1635)
&test_program($update, $continue, "mhmms", 
  "--global --quiet"
    . " mhmm/lipo.linear.mhmm"
    . " common/lipo.fasta",
    "mhmms/r1635.global.mhmms");

# Test mhmms for bug where runs of 
# spacer states are scored as matches (r1635).
# With theses input files mhmms should find no matches.
&test_program($update, $continue, "mhmms", 
  "--quiet --motif-scoring"
    . " mhmm/lipo.linear.mhmm"
    . " mhmms/r1635.fasta",
    "mhmms/r1635.mhmms");

# Test mhmms for segmentation fault when --motif-scoring option 
# is used. (r1584)
&test_program($update, $continue, "mhmms", 
  "--motif-scoring --quiet"
    . " mhmms/r1584.hmm"
    . " mhmms/r1584.fasta",
    "mhmms/r1584.mhmms");

# Test motiph
&test_program($update, $continue, "motiph", 
  "--bg 2.0 --pseudocount 0.01 --text "
    . " motiph/spiked.aln"
    . " motiph/yeast.tree"
    . " motiph/MCM1.meme.html",
    "motiph/motiph.gff");

# Test motiph with --motif option
&test_program($update, $continue, "motiph", 
  "--bg 2.0 --pseudocount 0.01 --text  "
    . " --output-pthresh 1.0"
    . " --motif 2 --motif 3"
    . " motiph/spiked.aln"
    . " motiph/yeast.tree"
    . " common/crp0.meme.html",
    "motiph/motiph-motif23.gff");

# Test psp-gen on DNA
&test_program($update, $continue, "psp-gen", 
  "-revcomp "
    . " -pos psp-gen/one-peak-dna.fasta"
    . " -neg psp-gen/all-A.fasta",
    "psp-gen/one-peak-dna-revcomp.psp",
    $testscriptdir);

# Test psp-gen on protein
&test_program($update, $continue, "psp-gen", 
  "-alpha prot -maxrange -triples "
    . " -pos psp-gen/one-peak-protein.fasta"
    . " -neg psp-gen/all-A.fasta",
    "psp-gen/one-peak-protein.psp",
    $testscriptdir);

# Test fimo
# The MEME_ETC_DIR contains the needed style sheets.
&test_program($update, $continue, "fimo", 
  "--text --motif-pseudo 0.01"
    . " motiph/MCM1.meme.html"
    . " motiph/spiked.fasta",
    "motiph/fimo.txt");

# Test fimo with --motif option
# The MEME_ETC_DIR contains the needed style sheets.
&test_program($update, $continue, "fimo", 
  "--text --motif 2 --motif 3"
    . " --output-pthresh 0.01"
    . " common/crp0.meme.html"
    . " motiph/spiked.fasta",
    "motiph/fimo-motif23.txt");

# Test ama without --z-scoring
&test_program($update, $continue, "ama",
  " --pvalues --verbosity 1" 
    . " gomo/motif.meme"
    . " gomo/seqs.fasta"
    . " gomo/seqs.norc.bg",
    "gomo/ama.nozscoring.xml");
    
# Test ama with maxodds scoring
&test_program($update, $continue, "ama",
  " --verbosity 1 --scoring max-odds"
    . " gomo/motif.meme"
    . " gomo/seqs.fasta"
    . " gomo/seqs.bg",
    "gomo/ama.withMaxodds.xml");

# Test gomo on single species
&test_program($update, $continue, "gomo", 
  " --nostatus --verbosity 1 --text" 
    . " gomo/GO2Gene.map.csv"
    . " gomo/ama.nozscoring.xml",
    "gomo/gomo.smallthreshold.txt");

# Test gomo on multiple species
&test_program($update, $continue, "gomo", 
  " --nostatus --verbosity 1 --text" 
    . " gomo/GO2Gene.map.csv"
    . " gomo/ama.nozscoring.xml"
    . " gomo/ama.nozscoring.xml",
    "gomo/gomo.multipeSpecies.txt");

# Test shadow 
&test_program($update, $continue, "shadow", 
  "--output-pthresh 0.1 --bg 2.0 --text"
    . " motiph/spiked.aln"
    . " motiph/yeast.tree",
    "motiph/shadow.gff");

# tomtom test is failing on some platforms. Issue with platfrom variation
# in random number generators? FIXME
# Test tomtom (7 distance measures)
# foreach $score (
#     "allr",
#     "ed",
#     "kullback",
#     "pearson", 
#     "sandelin",
#     "blic1",
#     "blic5"
#     ) {
#   &test_program($update, $continue, "tomtom", 
#     " -query common/sample.meme -target common/sample.meme"
#     . " -dist " . $score . " -text ",
#     "tomtom/tomtom.out." . $score);
# }

# qvalue test is failing on some platforms. Issue with platfrom variation
# in random number generators? FIXME
#&test_program($update, $continue, "qvalue",
#	      "--header 1 --append --column 2 --seed 7718 qvalue/uniform.txt", 
#	      "qvalue/uniform.out");

if (0) { # TLB; broken test
&test_program($update, $continue, "qvalue",
	      "--null qvalue/null.txt qvalue/observed.txt", 
	      "qvalue/observed.out");
}
exit(0);

###########################################################################
sub test_scaffold {
  my($fileroot, $is_dna) = @_;

  # Create input filenames.
  my $train_file = "common/$fileroot.fasta";
  my $meme_file = "common/$fileroot.meme.html";
  my $test_file = "common/$fileroot-test.fasta";
  my $linear_model = "mhmm/$fileroot.linear.mhmm";
  my $complete_model = "mhmm/$fileroot.complete.mhmm";
  my $test_aln = "common/test.aln";
  my $test_aln_out = "clustalw2fasta/test.fasta";
  my $test_aln_nogap_out = "clustalw2fasta/test.nogap.fasta";
  my $test_aln_consensus_out = "clustalw2fasta/test.consensus.fasta";

  # All files are in the scaffold directory.
  my $scaffold_root = "scaffold/$fileroot";

  my $dna_switch;
  if ($is_dna == 1) {
    $dna_switch = "--dna";
  } else {
    $dna_switch = "";
  }

  # Run a test on each program in succession.
  &test_program($update, $continue, "meme-io", " --verbosity 1 $meme_file ",
		"$scaffold_root.meme-io");
  &test_program($update, $continue, "fasta-io",
		" --verbosity 1 $dna_switch $train_file ", $train_file);
  &test_program($update, $continue, "fasta-io", 
		" --verbosity 1  --many $dna_switch $train_file", $train_file);
  &test_program($update, $continue, "fasta-io", 
		" --verbosity 1 --blocksize 100 $train_file", "$scaffold_root.fasta");

  # Test model I/O.
  &test_program($update, $continue, "mhmm-io", " --verbosity 1 $linear_model", 
		"$linear_model");
  &test_program($update, $continue, "mhmm-io", " --verbosity 1 $complete_model",
		"$complete_model");

  # Test model log conversion.
  &test_program($update, $continue, "log-hmm", " --verbosity 1 $linear_model", 
		"$linear_model");
  &test_program($update, $continue, "log-hmm", " --verbosity 1 $complete_model",
		"$complete_model");

  # Test clustalw2fasta
  &test_program($update, $continue, "clustalw2fasta", " $test_aln", 
    $test_aln_out);
  &test_program($update, $continue, "clustalw2fasta", " -nogap $test_aln", 
    $test_aln_nogap_out);
  &test_program($update, $continue, "clustalw2fasta", 
    " -consensus 100 $test_aln", $test_aln_consensus_out); 
}

###########################################################################
sub test_core {
  my($fileroot, $is_dna) = @_;

  # Create input filenames.
  my $train_file = "common/$fileroot.fasta";
  my $meme_file = "common/$fileroot.meme.html";
  my $test_file = "common/$fileroot-test.fasta";

  # Check different topologies.
  foreach $topology ("linear", "complete", "star") {

    # Create various types of spacer models.
    foreach $spacer ("1", "3", "fim") {

      # FIXME: There is a bug with star topology plus fims.
      if (($spacer eq "fim") && ($topology eq "star")) {
				next;
      }

      $model = "$fileroot.$topology";
      $mhmm_params = " --noheader --noparams --type $topology ";

      if ($spacer eq "fim") {
				$model .= ".fim";
				$mhmm_params .= "--fim ";
      } elsif ($spacer eq "3") {
				$model .= ".spacer";
				$mhmm_params .= "--nspacer 3";
      }
			$mhmm_params .= " --verbosity 1 $meme_file";
      
      # Create the model.
      &test_program($update, $continue, "mhmm", $mhmm_params, 
		    "mhmm/$model.mhmm");

      # Draw the model.
      &test_program($update, $continue, "draw-mhmm",
		    " --verbosity 1 --consensus mhmm/$model.mhmm ", 
		    "draw-mhmm/$model.gvz");

      # Set parameters for search routines.
      $search_params = "--ethresh 99999 --quiet --width 79 ";

      # Use different scoring schemes.
      foreach $paths ("all", "single") {

				# Search with the model.	    
				$search_file = "$model.$paths";
				if ($paths ne "all") {
	  			$search_params .= "--fancy ";
				}
				&test_program($update, $continue, "mhmms", 
		      " --verbosity 1 --paths $paths $search_params mhmm/$model.mhmm $test_file",
		      "mhmms/$search_file.mhmms");
      }
    }
  }
}

############################################################################
#
# programs that are run from src/ need not specify $bindir
#
############################################################################

sub test_program {
  my($update, $continue, $program, $arguments, $good_file, $bindir) = @_;
  my $output = "/var/tmp/$program.$$.tmp";
  
  $bindir = $testbindir unless $bindir;
  # Tell the user what's happening.
  printf("************************************************************\n");
  printf("Testing $program . . . \n");

  # Run the program and collect the error level.
  printf("$bindir/$program $arguments \n");
  system("$bindir/$program $arguments > $output");
  my $result = $? >> 8;

  # Eliminate any lines that might differ, but are irrelevant
  system("sed -e '/xml-stylesheet/d' $output > $output.sed");

  # Diff the results with the desired results.
  my $diff = `diff $good_file $output.sed`;

  # Check for problems and die if any happened.
  if (($diff ne "") || ($result != 0) || ($? >> 8 != 0)) {
    printf("FAIL\n");
    if ($result != 0) {
      printf("$program had a non-zero exit status ($result)\n");
    } else {
      print("$good_file\n");
      print($diff);
      if ($update) {
	print(STDERR "Updating $good_file.\n");
	`cp $output.sed $good_file`;
      }
    }
    if (!$continue) {
      exit(1);
    }
  } else {

    # Return success.
    printf("SUCCESS\n");
  }

  # Delete the temporary file.
  unlink($output);
  unlink("$output.sed");
}
