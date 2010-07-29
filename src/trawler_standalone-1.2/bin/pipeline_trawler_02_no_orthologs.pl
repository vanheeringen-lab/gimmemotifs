#!/usr/bin/perl

# $Id: pipeline_trawler_02_no_orthologs.pl,v 1.13 2009/04/27 11:49:10 haudry Exp $

=head1 NAME

    pipeline_trawler_02_no_orthologs.pl

=head1 SYNOPSIS

    If you do not have orthologous sequences this program takes the sample sequences file (in FASTA format)
    and the motif cluster file (from pipeline_01) and locate the motifs onto the sequences.
    Because no orthologous sequences are provided, no conservation scores will be calculate.
    Instead the sequences will be sorted according to the number of different sites found.

=head1 DESCRIPTION

    please use the help function : ./pipeline_trawler_02_no_orthologs.pl -help


=head1 OPTIONS

    -motif (motif cluster from trawler clustering .clsuter)
    -sequences(sample file containing the sequences in FASTA format)
    -help

=head1 CONTACT

Contact laurence Ettwiller (EMBL) ettwille@embl.de

=cut

use strict;
use Carp;
use File::Basename;
use File::Spec::Functions qw[catfile catdir];
use Getopt::Long;

# Locate Trawler modules
use FindBin ();
use lib "$FindBin::RealBin/../modules";
my $script_name = $FindBin::RealScript;

# Internal Modules
use Trawler::Constants 1.0 qw(_read_config _tcst);
use Trawler::Utils 1.0 qw(_parse_seq_file _parse_cluster_file _reg_exp);

# START processing
print "\n## Running $script_name\n";

#==============================================================================
# Read and Set properties
#==============================================================================

# Read config file
_read_config($FindBin::RealBin);
my %tcst = _tcst();

# Logging Levels
my $DEBUG = $tcst{DEBUG};
my $INFO  = $tcst{INFO};

# get the different options ===================================================
my $file_motif;
my $file_sequences;
my $help = 1;
my $DIR = undef;

GetOptions(
  'motif:s'     => \$file_motif,
  'sequences:s' => \$file_sequences,
  'directory=s' => \$DIR,
  'help:s'      => \$help,
);

if (!$help || !$file_motif || !$file_sequences) {

    print STDERR "\nUSAGE : \n \n perl pipeline_trawler_02_no_orthologs.pl  -motif [the complete path to the XX.cluster file from pipeline_trawler1.pl ]  -sequences [complete path to the sequences sequence in fasta format (sample file in pipeline trawler1 ] ]\n\n for example \n\n
perl pipeline_trawler_02_no_orthologs.pl -motif /path/ip-test/ip-test.cluster -sequences /path/sequences/sample.fasta\n\n

===========================================\n===========================================\n \n";
    exit(1);
}

#==============================================================================
# Directories handling

# if no directory provided => exit [must set up in trawler.pl]
unless($DIR) {
    print STDERR "\nERROR: no directory to run\n";
    exit(1);
}
# extract tmp_xxx
my $tmp_dir_name = fileparse($DIR);
# default result directory is like $TRAWLER_HOME/tmp_YYYY-MM-DD_HHhmm:ss/result
my $tmp_result_dir = File::Spec->catdir( $DIR, $tcst{RES_DIR_NAME} );
my $tmp_fasta_dir = File::Spec->catdir( $tmp_result_dir, $tcst{FASTA_DIR_NAME} );
my $tmp_features_dir = File::Spec->catdir( $tmp_result_dir, $tcst{FEATURES_DIR_NAME} );

if ($DEBUG) {
  print "the data will be stored into $tmp_result_dir \n";
  print "file_seq: $file_sequences \n";
  print "file_motif: $file_motif \n";
}

#==============================================================================
# Run
msg_pipeline();

my $id2sequences = _parse_seq_file($file_sequences);

#======================end getting options====================================

my $id2motif;

if ($file_motif) {
    $id2motif = _parse_cluster_file($file_motif);
}

#==============================================================================

# STAT file: $TRAWLWER_RESULT/features/STAT.tmp
my $file_stat = File::Spec->catfile($tmp_features_dir, $tcst{STAT_FILE_NAME});
open(STAT, ">$file_stat") or croak "Cannot open file $file_stat: $!";
print STAT "#id\tmotif_id\tmotif\tstart\taverage_conservation\tconserved_count\tsequence_length\tstrand\tstart_from_end\n";

my %FEATURES;

foreach my $id (keys %$id2sequences) {
    # fasta file: $TRAWLER_RESULT/fasta/<id>.fasta
    my $file = File::Spec->catfile($tmp_fasta_dir, $id . $tcst{FASTA_FILE_EXT});

    open(OUT, ">$file") or croak "Cannot open file $file: $!";
    my $ref_seq = $$id2sequences{$id};
    my $length_seq = length($ref_seq);
    print OUT ">$id\n$ref_seq\n";
    foreach my $idmotif (keys %$id2motif) {
      my @motifs = @{$$id2motif{$idmotif}};

        foreach my $motif (@motifs) {
            #[YH:remove] get_stat_motif($ref_seq, $motif, $id, $idmotif);
            my ($mot, $mstarts, $mends, $strands) = get_motif_loc($ref_seq, $motif);
            if($mot) {
                my $size_mot = @$mot;
                for(my $incr = 0; $incr<$size_mot; $incr++) {
                    my $m_start = $$mstarts[$incr];
                    my $m_end = $$mends[$incr];
		            my $m_strand = $$strands[$incr];
                    my $line_feature = "motif\t".$id."\t1\t".$m_start."\t".$m_end."\tmotif\n";

                    #[YH] Replace get_stat_motif function
                    my $start_from_end = ($length_seq - $m_start) + 1;
                    print STAT "$id\t$idmotif\t$motif\t$m_start\t0\t0\t$length_seq\t$m_strand\t$start_from_end\n";

                    push @{$FEATURES{$idmotif}->{$id}->{$motif}->{"start"}}, $m_start;
                    push @{$FEATURES{$idmotif}->{$id}->{$motif}->{"end"}}, $m_end;
                }
            }
        }
    }
    close(OUT) or croak "Can't close file '$file': $!";
}
close(STAT) or croak "Can't close file '$file_stat': $!";

#==========================================================
#deal with features
#==========================================================

foreach my $file1 (keys %FEATURES) {
    my $tmp = $FEATURES{$file1};

    foreach my $file2 (keys %$tmp) {
        # feature file: $TRAWLER_RESULT/features/<file1>_<file2>.txt
        my $feature_file = File::Spec->catfile($tmp_fasta_dir, $file1 . "_" . $file2. ".txt");
        open(FEATURE, ">$feature_file") or croak "Cannot open file $feature_file: $!";
        print FEATURE "$file1\tff00ff\nmotif\t009ba5\n";

        foreach my $feature (keys %FEATURES) {
            my $id2motif = $FEATURES{$feature};

            foreach my $id (keys %$id2motif) {
                my $motif2loc = $$id2motif{$id};

                foreach my $motif (keys %$motif2loc) {
                    my @starts = @{$$motif2loc{$motif}->{"start"}};
                    my @ends = @{$$motif2loc{$motif}->{"end"}};
                    my $l = @starts;

                    for (my $i=0; $i<$l; $i++) {
                        my $s = $starts[$i];
                        my $e = $ends[$i];

                        if($id eq $file2 && $feature ne $file1) {
                            print FEATURE "$feature\t$id\t1\t$s\t$e\tmotif\n";
                        }
                        elsif ($feature eq $file1 && $id eq $file2) {#same family therefore same color
                            print FEATURE "$motif\t$id\t1\t$s\t$e\t$file1\n";
                        }
                    }
                }
            }
        }
        close(FEATURE) or croak "Can't close file '$feature_file': $!";
    }
}


#==============================================================================
# Subroutines
#==============================================================================

#called by get_motif_position
sub get_motif_loc {
    my($DNA, $motif) = @_;
    my $motif_rc = reverse_complement($motif);
    my $pattern_matching1 = _reg_exp($motif);
    my $pattern_matching2 = _reg_exp($motif_rc);
    #simply a pattern matching of the motif on the sequence
    my $matches = 0;
    my $motif_length = length($motif);
    my @motifs; my @start; my @end; my @strand;
    while ($DNA =~ m/($pattern_matching1)/g) { #the motif
        my $m = $1;
        my $motif_length = length($m);
        my $pos = pos $DNA;
        push @start, ($pos-$motif_length);
        push @end, $pos;
        push @motifs, $m;
	    push @strand, "1";
    }
    while ($DNA =~ m/($pattern_matching2)/g) { #the reverse complement motif
        my $m = $1;
        my $motif_length = length($m);
        my $pos = pos $DNA;
        push @start, ($pos-$motif_length);
        push @end, $pos;
        push @motifs, $m;
	    push @strand, "-1";
    }
    return (\@motifs, \@start, \@end, \@strand) ;
}

#get the statistique about the conservation of motifs.
#sub get_stat_motif {
#    my ($seq,  $motif, $id, $motif_id) = @_;
#    #print "*** CALL get_stat_motif\n";
#
#    my ($motifs, $starts, $ends, $strands) = get_motif_loc($seq, $motif);
#    my $size = @$motifs;
#
#    for(my $i=0; $i< $size; $i++) {
#        my $start = $$starts[$i];
#        my $end = $$ends[$i];
#        my $dif = $end - $start;
#
#        print STAT "$id\t$motif_id\t$motif\t$start\t0\t0\n";
#        #print "$id\t$motif_id\t$motif\t$start\t0\t0 ### Func \n";
#
#    }
#}

sub create_array {
    my ($seq) = @_;
    print "caS: $seq\n" if $DEBUG;
    my @result;
    my $length = length($seq);
    print "caL: $seq\n";

    for (my $i=0; $i<$length; $i++) {
        $result[$i] = "-";
    }
    print "caR: @result\n" if $DEBUG;
    return \@result;
}

sub reverse_complement {
    my ($patt) = @_;
    my @array = split //, $patt;
    my $size = @array;
    my @reverse;

    for (my $i=$size-1; $i>=0; $i--) {
        my $nt = $array[$i];
        my $comp_nt = complement($nt);
        push @reverse, $comp_nt;
    }
    my $result = join "", @reverse;
    return $result;
}

sub complement {
    my ($nt) = @_;

    my $result;

    if($nt eq "A") { $result = "T"; }
    elsif($nt eq "T") { $result = "A"; }
    elsif($nt eq "C") { $result = "G"; }
    elsif($nt eq "G") { $result = "C"; }
    elsif($nt eq "M") { $result = "K"; }
    elsif($nt eq "K") { $result = "M"; }
    elsif($nt eq "R") { $result = "Y"; }
    elsif($nt eq "Y") { $result = "R"; }
    elsif($nt eq "V") { $result = "B"; }
    elsif($nt eq "B") { $result = "V"; }
    elsif($nt eq "D") { $result = "H"; }
    elsif($nt eq "H") { $result = "D"; }
    elsif($nt eq "W") { $result = "W"; }
    elsif($nt eq "S") { $result = "S"; }
    elsif($nt eq ".") { $result = "."; }
    elsif($nt eq "N") { $result = "N"; }
    else { $result = "n"; }

    return $result;
}

sub msg_pipeline {
  my $msg_pipeline = <<"MSG";
  ==========
  Processing motifs...
  ==========
MSG

  print $msg_pipeline . "\n";
}

# TODO[YH]: method not used, should be removed ?
#sub get_motif_positions {
#    my ($motif, $conserved_seq) = @_;
#    my $array_null = create_array($conserved_seq);
#    #first look at the real sequences iterate through all the sequences
#    my $total_real = 0;
#    my ($motifs, $starts, $ends) = get_motif_loc($conserved_seq, $motif);
#    my $size = @$motifs;
#    for(my $i=0; $i< $size; $i++) {
#        my $start = $$starts[$i];
#        my $end = $$ends[$i];
#
#        for (my $j=$start; $j<$end; $j++) {
#            my $nt = substr($conserved_seq, ($j), 1);
#
#            $$array_null[$j]=$nt;
#        }
#    }
#    my $seq = join "", @$array_null;
#    if ($seq =~/[A-Z]/) {
#      return $seq;
#    }
#    else {return}
#}

# TODO[YH]: method not used, should be removed ?
#sub get_stat_sequence {
#    my ($ref_seq, $species2seq) = @_;
#
#    my $size = length($ref_seq);
#    my $stat = 0;
#    for (my $i=0; $i<$size; $i++) {
#        my $conservation_count =0;
#        my $nt = substr($ref_seq, $i, 1);
#        foreach my $specie (keys %$species2seq) {
#            my $seq = $$species2seq{$specie};
#            my $nt_to_compare = substr($seq, $i, 1);
#            if ($nt_to_compare eq $nt && $nt =~/[ATCG]/) {
#                $conservation_count ++;
#            }
#        }
#
#        $stat = $stat+$conservation_count;
#    }
#
#    $stat = $stat / $size;
#    return $stat;
#}

1;
