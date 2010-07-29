#!/usr/bin/perl

# $Id: pipeline_trawler_02_single_strand.pl,v 1.4 2009/04/27 11:49:10 haudry Exp $

=head1 NAME

     pipeline_trawler_02.pl

=head1 SYNOPSIS


=head1 DESCRIPTION

    takes the result of the motif clustering (from trawler_01.pl),
    a list of alignment files and the reference species and locate the motif on the reference species.

=head1 OPTIONS

    -motif (motif cluster from trawler clustering .clsuter)
    -list (file contaning the path to all the alignment in FASTA format. The name of the file is the subsequent ID)
    -dir (name of the directory, default tmp/)
    -ref_species #name of the ref species (id of the sequences) >species_is\nseq
    -help

=head1 CONTACT

Contact laurence Ettwiller (EMBL) ettwille@embl.de


=head1 DEPENDENCY

No dependency.

=cut


####
# print STDERR "arrived a specie $specie for the alignment\n"; [470]
# print STDERR "$id\t$motif_id\t$motif\t$start\t$average_conservation\t$conserved_count\n"; [353]
####


use strict;
use Carp;
use File::Basename;
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

#==========================get the different options

my $DIR         = undef;
my $file_motif  = undef;
my $file_list   = undef;
my $ref_species = undef;
my $help        = 1;

GetOptions(
  'motif:s'       => \$file_motif,
  'list:s'        => \$file_list,
  'ref_species:s' => \$ref_species,
  'directory=s'   => \$DIR,
  'help:s'        => \$help,
);

if (!$help || !$file_motif || !$file_list || !$ref_species ) {
     print STDERR "\nUSAGE : \n \n perl pipeline_trawler_02.pl  -motif [the complete path to the XX.cluster file from pipeline_trawler1.pl ]  -list [list of the files contaning the alignment. Each file is a locus and the file name is the id of the locus. ] -ref_species [exact name of the reference species]\n\n for example \n\n

perl pipeline_trawler_02.pl -motif /path/ip-test/ip-test.cluster -list /path/orthologous_sequences/list -ref_species dr6\n\n

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

#==============================================================================
# Run
msg_pipeline();

my $id2specie2sequences = parse_list_file($file_list);

#======================end getting options====================================

my $id2motif;
if ($file_motif) {
    $id2motif = parse_cluster_file($file_motif);
}

#==============================================================================

# features/STAT.tmp file
my $file_stat = File::Spec->catfile($tmp_features_dir, $tcst{STAT_FILE_NAME});
open(STAT, ">$file_stat") or croak "Cannot open file $file_stat: $!";
print STAT "#id\tmotif_id\tmotif\tstart\taverage_conservation\tconserved_count\tsequence_length\tstrand\tstart_from_end\n";

my %FEATURES;
my %ref_id2seq;
foreach my $id (keys %$id2specie2sequences) {

    #GET THE REF SEQUENCE
    if ($$id2specie2sequences{$id}->{$ref_species}) { #put back to "reference"
        my $ref_seq = $$id2specie2sequences{$id}->{$ref_species};#put back to "reference"
        delete  $$id2specie2sequences{$id}->{$ref_species};#put back to "reference"
        my $species2seq = $$id2specie2sequences{$id};
        $ref_id2seq{$id}=$ref_seq;
        # fasta file = $TRAWLER_RESULT/fasta/<id>.fasta
        my $file = File::Spec->catfile($tmp_fasta_dir, $id . $tcst{FASTA_FILE_EXT});
        open(OUT, ">$file") or croak "Cannot open file $file: $!";

        my $average_conservation = get_stat_sequence($ref_seq, $species2seq);

        foreach my $specie (keys %$species2seq) {
            my $seq_out = $$species2seq{$specie};
            if ( $seq_out) {
                print OUT ">$specie\n$seq_out\n";
            }
        }
        print OUT ">$ref_species\n$ref_seq\n";

        foreach my $idmotif (keys %$id2motif) {

            my @motifs = @{$$id2motif{$idmotif}};

            foreach my $motif (@motifs) {

                get_stat_motif($ref_seq, $species2seq, $motif, $id, $idmotif, $average_conservation);
                my ($mot, $mstarts, $mends) = get_motif_loc($ref_seq,$motif);
                if($mot) {
                    my $size_mot = @$mot;
                    for(my $incr =0; $incr<$size_mot; $incr++) {
                        my $m_start = $$mstarts[$incr];my $m_end = $$mends[$incr];
                        my $line_feature = "motif\t".$id."\t1\t".$m_start."\t".$m_end."\tmotif\n";
                        #print FEATURE "motif\t$id\t1\t$m_start\t$m_end\tmotif\n";
                        push @{$FEATURES{$idmotif}->{$id}->{$motif}->{"start"}}, $m_start;
                        push @{$FEATURES{$idmotif}->{$id}->{$motif}->{"end"}}, $m_end;
                    }
                }
            }
        }
        close(OUT) or croak "Can't close file '$file': $!";
    }

}
close(STAT) or croak "Can't close file '$file_stat': $!";

#==========================================================
#deal with features
#==========================================================

foreach my $file1 (keys %FEATURES) {
    print "$file1\n" if $DEBUG;
    my $tmp = $FEATURES{$file1};
    foreach my $file2 (keys %$tmp) {
        my $ref_seq = $ref_id2seq{$file2};
        print "$ref_seq\n" if $DEBUG;
        # FIXME[YH]: to be removed
        # my $feature_file = $DIR."/fasta/".$file1."_".$file2.".txt";
        my $feature_file = File::Spec->catfile($tmp_fasta_dir, $file1."_".$file2.".txt");
        open(FEATURE, ">$feature_file") or croak "Cannot open file $feature_file: $!";
        print FEATURE "$file1\tff00ff\nmotif\t009ba5\n";
        print "$file1\tff00ff\nmotif\t009ba5\n" if $DEBUG;

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

                        my $sub_seq = substr($ref_seq, 0, $s);
                        $sub_seq =~ s/\w+//g;
                        my $diff = length($sub_seq);
                        #print "$sub_seq\n" if $DEBUG;
                        #print "$diff\n" if $DEBUG;
                        $s = $s - $diff;
                        $e = $e - $diff;

                        if($id eq $file2 && $feature ne $file1) {
                            print FEATURE "$feature\t$ref_species\t1\t$s\t$e\tmotif\n";
                        }
                        elsif ($feature eq $file1 && $id eq $file2) { #same family
                            print FEATURE "$motif\t$ref_species\t1\t$s\t$e\t$file1\n";
                        }
                    }
                }
            }
        }
        close FEATURE;
    }
}


#subroutines=============================================

sub parse_cluster_file {
    my ($file) = @_;
    my %result;

    open(FILE, $file) or croak "Cannot open cluster file $file: $!";
    foreach my $line (<FILE>) {
        chomp $line;
        if ($line !~/^java/) {
            my @tmp = split /\s+/, $line;
            my $size = @tmp;

            #if ($size ==5) { #it is from a cluster of more than one motif
            my $cluster_name = $tmp[1];
            my $motif = $tmp[0];

            if ($cluster_name) {
                push @{$result{$cluster_name}}, $motif;
            }
            #}
        } # end escape java command line
    } # end processing lines
    close(FILE) or croak "Cannot close cluster file $file: $!";

    return \%result;
}

sub get_motif_positions {
    my ($motif, $conserved_seq) = @_;
    my $array_null = create_array($conserved_seq);
    #first look at the real sequences iterate through all the sequences
    my $total_real = 0;
    my ($motifs, $starts, $ends) = get_motif_loc($conserved_seq, $motif);
    my $size = @$motifs;
    for(my $i=0; $i< $size; $i++) {
        my $start = $$starts[$i];
        my $end = $$ends[$i];

        for (my $j=$start; $j<$end; $j++) {
         my $nt = substr($conserved_seq, ($j), 1);
         $$array_null[$j]=$nt;
        }
    }
    my $seq = join "", @$array_null;
    if ($seq =~/[A-Z]/) {
        return $seq;
    }
    else {return}
}

sub get_motif_loc {
    my($DNA, $motif) = @_;
    print "motif $motif\n" if $DEBUG;

    my $pattern_matching1 = _reg_exp($motif);

    print "motif $motif $pattern_matching1 \n" if $DEBUG;
    my $matches = 0;
    my $motif_length = length($motif);
    my @motifs; my @start; my @end;

    while ($DNA =~ m/($pattern_matching1)/g) {
        my $m = $1;
        print "$m\n" if $DEBUG;
        my $motif_length = length($m);
        my $pos = pos $DNA;
        push @start, ($pos-$motif_length);
        push @end, $pos;
        push @motifs, $m;
    }


    print "@motifs @start @end\n" if $DEBUG;
    return (\@motifs, \@start, \@end) ;

}

sub get_stat_sequence {
    my ($ref_seq, $species2seq) = @_;

    my $size = length($ref_seq);
    my $stat = 0;
    for (my $i=0; $i<$size; $i++) {
        my $conservation_count =0;
      my $nt = substr($ref_seq, $i, 1);
      foreach my $specie (keys %$species2seq) {
            my $seq = $$species2seq{$specie};
            my $nt_to_compare = substr($seq, $i, 1);
            if ($nt_to_compare eq $nt && $nt =~/[ATCG]/) {
                $conservation_count ++;
            }
        }
        $stat = $stat+$conservation_count;
    }

    $stat = $stat / $size;
    return $stat;
}

#get the statistique about the conservation of motifs.
sub get_stat_motif {
    my ($seq, $id2seq, $motif, $id, $motif_id, $average_conservation) = @_;

    my ($motifs, $starts, $ends) = get_motif_loc($seq, $motif);
    my $size_seq = length($seq);
    my $size = @$motifs;
    my $motif_rc = reverse_complement($motif);
    my $pattern_matching1 = _reg_exp($motif);
    my $pattern_matching2 = _reg_exp($motif_rc);
    for(my $i=0; $i< $size; $i++) {
        my $conserved_count =0;
        my $start = $$starts[$i];
        my $end = $$ends[$i];
        my $dif = $end - $start;
        foreach my $id_species (keys %$id2seq) {
            my $conserved_seq = $$id2seq{$id_species};
            my $sub_seq = substr($conserved_seq, $start, $dif);

            if ($sub_seq =~/$pattern_matching1/ || $sub_seq =~/$pattern_matching2/) {
                $conserved_count++;
            }
        }
        my $start_from_end = ($size_seq - $start) + 1;
        print STAT "$id\t$motif_id\t$motif\t$start\t$average_conservation\t$conserved_count\t$size_seq\t1\t$start_from_end\n";
        print "$id\t$motif_id\t$motif\t$start\t$average_conservation\t$conserved_count\t$size_seq\n" if $DEBUG;
    }

} # /get_stat_motif

sub parse_jaspar_list {
    my ($file) = @_;
    my %id2jasp;
    open(FILE_DIR, $file) or croak "Cannot open file $file: $!";
    foreach my $line (<FILE_DIR>) {
        chomp $line;
        my @result;
        open(FILE, $line) or croak "Cannot open file $line: $!";
        my $count = 0;
        $line =~ s/\/home\/ettwille\/mef2_project\/JASPAR\///;

        $line =~ s/\.pfm//;
        foreach my $line(<FILE>) {
            chomp $line;

            my @tmp = split /\s+/, $line;
            my $size = @tmp;
            for (my $i=0; $i<$size; $i++) {
                $result[$count][$i] = $tmp[$i];
            }
            $count++;
        }
        close(FILE) or croak "Can't close file '$line': $!";
        $id2jasp{$line} = \@result;
    }
    close(FILE_DIR) or croak "Can't close file '$file': $!";
    return \%id2jasp;
}

sub  parse_cutoff_file {
    my ($file) = @_;

    my %result;
    open(FILE, $file) or croak "Cannot open cutoff file $file_stat: $!";
    foreach my $line (<FILE>) {
        chomp $line;
        my ($jaspar, $cutoff) = split /\s+/, $line;
        $result{$jaspar} = $cutoff;
    }
    close(FILE) or croak "Can't close file '$file_stat': $!";
    return \%result;
}

sub get_pwm_positions {
    my($id,$conserved_seq, $pwm, $cutoff) = @_;

    my $array_null = create_array($conserved_seq);
    #first look at the real sequences iterate through all the sequences
    my $total_real = 0;
    my ($motifs, $starts, $ends) = run_pwm_return_hits($conserved_seq, $pwm, $cutoff);
    my $size = @$motifs;
    for(my $i=0; $i< $size; $i++) {
        my $start = $$starts[$i];
        my $end = $$ends[$i];

        for (my $j=$start-1; $j<$end; $j++) {
            my $nt = substr($conserved_seq, ($j), 1);
            $$array_null[$j]=$nt;
        }
    }
    my $seq = join "", @$array_null;
    if ($seq =~/[A-Z]/) {
        return $seq;
    }
    else {return}
}

sub parse_transfac_file {
    my ($file) = @_;

    my %result; my @tmp;
    open(FILE, $file) or croak "Cannot open cutoff file $file: $!";
    foreach my $line (<FILE>) {
        chomp $line;
        my ($transfac, $cutoff) = split /\s+/, $line;
        push @tmp, $transfac;
        $result{$transfac} = $cutoff;
    }
    close(FILE) or croak "Can't close file '$file': $!";
    return (\%result, \@tmp);

}

sub run_pwm_return_hits {
    my ($DNA, $pwm, $cut_off) = @_;

    my $seq = Bio::Seq->new( -display_id => 'my_id',
           -seq => $DNA);

    my $siteset = $pwm->search_seq(-seqobj=>$seq, -threshold=>$cut_off);
    my $siteset_iterator = $siteset->Iterator(-sort_by =>'score');
    my @result1;my @result2;my @result3;
    while (my $site = $siteset_iterator->next) {
        my $start = $site->start;
        my $end = $site->end;
        my $motif = $site->seq()->seq;

        push @result1, $motif;
        push @result2, $start;
        push @result3, $end;
    }
    return (\@result1, \@result2, \@result3);
}

# FIXME[YH]: make use of File::Spec
sub parse_list_file {
    my ($list_file) = @_;

    my %result;
    open(LIST_FILE, $list_file) or croak "Cannot open list file $list_file: $!";
    foreach my $file (<LIST_FILE>) {
        chomp $file;
        my $id = $file;
        $id =~ s/\S+\///g; #remove the directory path !
        $id =~ s/\.\S+//g; #remove extentions !

        my $species2seq = parse_fasta_file($file);
        foreach my $species (keys %$species2seq) {
            my $seq = $$species2seq{$species};
           $result{$id}->{$species}=$seq;
        }
    }
    close(LIST_FILE) or croak "Can't close file '$list_file': $!";
    return (\%result)

}

sub parse_fasta_file {
    my ($fasta_file) = @_;

    my %result;
    open (FASTA_FILE, $fasta_file)  or croak "Cannot open list file $fasta_file: $!";
    local $/ = "\n>";
    foreach my $entry (<FASTA_FILE>) {

        my ($top,$sequence) = $entry =~ /^(.+?)\n([^>]*)/s
         or die("Can't parse entry fasta from $fasta_file [$entry]");
        my ($tmp,$fulldesc) = $top =~ /^\s*(\S+)\s*(.*)/
         or die("Can't parse fasta heade");
        $sequence =~ s/\n+//g;

        my @tmp_array = split /\|/, $tmp;
        my $id = $tmp_array[0];

        #remove characters that can't be used in the file name!
        #$id =~ s/\_//g;
        $id =~ s/\>//g;
        $id =~ s/\///g;
        $id =~ s/\\//g;

        #$sequence =~ s/[a-z]/N/g; #exons masked - we need to deal with that later.
        $result{$id}=$sequence;

    }
    close(FASTA_FILE) or croak "Can't close file '$fasta_file': $!";
    return \%result;
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
    my $result = join '', @reverse;
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

1;
