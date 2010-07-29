# $Id: Utils.pm,v 1.7 2009/05/07 09:44:28 haudry Exp $

=head1 NAME

Trawler::Utils - Trawler Utility methods

=head1 DESCRIPTION

This module ...

=head1 CONTACT

Yannick Haudry - EMBL (haudry@embl.de)

=cut

package Trawler::Utils;

use strict;
use Carp;

require Exporter;
our $VERSION = '1.00';
our @ISA = qw(Exporter);

our @EXPORT = (); # Symbols to autoexport (:DEFAULT tag)
our @EXPORT_OK = qw(_parse_DNA_file _reg_exp _reverse_complement _get_motif_loc _calculate_overlap
                    _parse_seq_file _parse_cluster_file); # Symbols to export on request()
our %EXPORT_TAGS = (); # Define names for sets of symbols - eg: TAG => [ qw(!name1 name2) ]

#==============================================================================
# subroutines extracted from pipeline_trawler_01_som.pl
#==============================================================================

# Used in pipeline_trawler_01.pl and pipeline_trawler_01_som.pl
#this is the parser for fasta !
sub _parse_DNA_file {
    my ($fasta_file) = @_;

    open(FASTA_FILE, $fasta_file) or croak "Cannot open fasta file $fasta_file: $!";

    local $/ = "\n>";
    my $DNA;
    foreach my $entry (<FASTA_FILE>) {

        my ($top, $sequence) = $entry =~ /^(.+?)\n([^>]*)/s
           or croak("Can't parse entry fasta from $fasta_file. You files do not appear to be in fasta format.");

        my ($id, $fulldesc) = $top =~ /^\s*(\S+)\s*(.*)/
           or croak("Can't parse fasta header ");
        $sequence =~ s/\n+//g; #remove the \n at each lines
        if ($sequence =~ /[^ATCGatcgnN]/) {
        	croak("this is not a nucleotid alphabets - check if there is other characters than ATCG or N (upper or lower case)\n")
        }
        $DNA = $DNA."X".$sequence;

    }
    close(FASTA_FILE) or croak "Can't close file '$fasta_file': $!";

    return $DNA;
}

#------------------------------------------------------------------------------

# revised version of reg_exp with elsif loops, is this OK ?
sub _reg_exp {
    my ($motif) = @_;

    my $regexp;
    my @tmp = split //, $motif;

    my $reg;
    foreach my $letter (@tmp) {
        LETTER: {
            if ($letter eq 'A') { $reg='A'; last LETTER; }
            elsif ($letter eq 'C') { $reg='C'; last LETTER; }
            elsif ($letter eq 'G') { $reg='G'; last LETTER; }
            elsif ($letter eq 'T') { $reg='T'; last LETTER; }
            elsif ($letter eq 'M') { $reg='[AC]'; last LETTER; }
            elsif ($letter eq 'R') { $reg='[AG]'; last LETTER; }
            elsif ($letter eq 'W') { $reg='[AT]'; last LETTER; }
            elsif ($letter eq 'S') { $reg='[CG]'; last LETTER; }
            elsif ($letter eq 'Y') { $reg='[CT]'; last LETTER; }
            elsif ($letter eq 'K') { $reg='[GT]'; last LETTER; }
            elsif ($letter eq 'V') { $reg='[ACG]'; last LETTER; }
            elsif ($letter eq 'H') { $reg='[ACT]'; last LETTER; }
            elsif ($letter eq 'D') { $reg='[AGT]'; last LETTER; }
            elsif ($letter eq 'B') { $reg='[CGT]'; last LETTER; }
            elsif ($letter eq 'N') { $reg='[ACGT]'; last LETTER; }
            elsif ($letter eq '.') { $reg='[ACGT]'; last LETTER; }
            else { $reg = "\."; }
        }
        $regexp .= "$reg";
    }
    return $regexp;
}

# used in pipeline_trawler_01.pl and pipeline_trawler_01_som.pl
sub _reverse_complement {
    my ($patt) = @_;

    my $rv = reverse $patt;
    $rv =~ tr/ACTGMKRYVBHD/TGACKMYRBVDH/;

    return $rv;
}

#------------------------------------------------------------------------------

# used in pipeline_trawler_01.pl and pipeline_trawler_01_som.pl
# WARNING different implementation in 02 pipelines !!

#IF YOU CHANGE THIS METHOD, ALWAYS MAKE SURE THAT RETURNED POSITION ARRAY ARE ORDERED ASC!!!
sub _get_motif_loc {
    my($DNA, $motif) = @_;

    my $matches = 0;
    my $motif_length = length($motif);
    my @motifs;
    my @start;
    #CHARLES: removed extension of match
    while ($DNA =~ m/($motif)/g) {
        my $m = $1;
        my $motif_length = length($m);
        my $pos = pos $DNA;
        $matches++;
        push @start, ($pos-$motif_length);
        push @motifs, $m;
    }
    return (\@motifs, \@start);
}

#------------------------------------------------------------------------------


#ARRAY1 AND ARRAY2 MUST BE SORTED ASC
sub _calculate_overlap {
    my ($array1, $array2, $motif1, $motif2) = @_;

    my $overlap = 0;
    my $size1 = @$array1;
    my $size2 = @$array2;
    my $m1len = length($motif1);
    my $m2len = length($motif2);

    #CHARLES: MODIF TO EVALUATE REAL OVERLAP
    foreach my $loc1 (@$array1) {
        foreach my $loc2 (@$array2) {

            next if ($loc2 + $m2len) < $loc1;
            last if $loc2 > $loc1+$m1len;

            if ( ($loc1 >= $loc2 && $loc1 < $loc2+$m2len) || ($loc2>=$loc1 && $loc2 < $loc1+$m1len)) {
                my $olapstart = $loc1>=$loc2 ? $loc1 : $loc2;
                my $olapstop = ($loc1+$m1len) <= ($loc2+$m2len) ? ($loc1+$m1len) : ($loc2+$m2len);
                #CHARLES: should we consider a min olap percentage to consider this olap real?
                #my $olaplen = $olapstop-$olapstart;
                $overlap++;
            }
        }
    }
    if ($size1 > 0 && $size2 > 0) {
        my $percent1 = ($overlap/$size1)*100;
        my $percent2 = ($overlap/$size2)*100;
        my $percent;
        if ($percent1 > $percent2) {
            $percent = $percent1;
        }
        else {
            $percent = $percent2;
        }
        return $percent;
    }
    else {
        return 0;
    }
}

#------------------------------------------------------------------------------

#==============================================================================
# subroutines extracted from pipeline_trawler_02_no_orthologs.pl
#==============================================================================

sub _parse_seq_file {
    my ($fasta_file) = @_;

    my %result;
    open(FASTA_FILE, $fasta_file) or croak "Cannot open fasta file $fasta_file: $!";

    local $/ = "\n>";
    foreach my $entry (<FASTA_FILE>) {

        my ($top,$sequence) = $entry =~ /^(.+?)\n([^>]*)/s
           or die("Can't parse entry fasta from $fasta_file [$entry]");
        my ($id,$fulldesc) = $top =~ /^\s*(\S+)\s*(.*)/
           or die("Can't parse fasta header");
        $sequence =~ s/\n+//g;
        #remove characters that can't be used in the file name!
        $id =~ s/\_//g;
        $id =~ s/\>//g;
        $id =~ s/\///g;
        $id =~ s/\\//g;

        $result{$id} = $sequence;

    }
    close(FASTA_FILE) or croak "Can't close file '$fasta_file': $!";

    return \%result;
}

#------------------------------------------------------------------------------

sub _parse_cluster_file {
    my ($file) = @_;

    my %result;

    open(FILE, $file) or croak "can't open motif file $file: $!";
    foreach my $line (<FILE>) {
        chomp $line;
        if ($line !~/\#/) { #this is not the first line of the file
            my @tmp = split /\t/, $line;
            my $size = @tmp;
            my $cluster_name = $tmp[1];
            my $motif = $tmp[0];
            if ($cluster_name) {
                push @{$result{$cluster_name}}, $motif;
            }
        }
    }
    close(FILE) or croak "Can't close file '$file': $!";

    return \%result;
}

#------------------------------------------------------------------------------

#warn "Trawler::Utils successfully loaded!\n";
1;
__END__
