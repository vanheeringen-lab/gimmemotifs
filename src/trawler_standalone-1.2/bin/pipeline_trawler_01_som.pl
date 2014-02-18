#!/usr/bin/perl

# $Id: pipeline_trawler_01_som.pl,v 1.28 2009/04/23 15:15:42 haudry Exp $

=head1 NAME

  pipeline_trawler_01.pl

=head1 SYNOPSIS

    this program takes 2 files [1] the sample file (in FASTA format) and [2] the background file (FASTA format)
    and will [1] run trawler and cluster all the motifs.
    the parameters used are described in the OPTIONS.

=head1 DESCRIPTION

This program is the main program runing Trawler.

=head1 OPTIONS

    -sample (FASTA format) better to be repeat-masked.
    -background (FASTA format)

    OPTIONAL PARAMETERS
    ===================

    [MOTIF DISCOVERY]
    -occurrence (optional) minimum occurrence in the sample sequences. [DEFAULT = 10]
    -mlength (optional) minimum motif length. [DEFAULT = 6]
    -wildcard (optional) number of wild card in motifs. [DEFAULT = 2]

    [CLUSTERING]
    -overlap (optional) in percentage. [DEFAULT = 70]
    -motif_number (optional) total number of motifs to be clustered. [DEFAULT = 200]
    -nb_of_cluster (optional) fixed number of cluster. if this option is used, the k-mean clustering algorithm with fixed k will be used instead of the self organizing map (SOM). [DEFAULT = NULL]

    [VISUALIZATION]
    -directory (optional) output directory. [DEFAULT = "TRAWLER_HOME/myResults"]
    -xtralen (optional) add bases around the motifs for the logo. [DEFAULT = 0]

=head1 CONTACT

Contact laurence Ettwiller (EMBL) ettwille@embl.de

=head1 DEPENDENCY

    perl
    java (http://java.sun.com/javase/downloads/index.jsp)
    trawler (from pecan)
    weblogo (http://weblogo.berkeley.edu/)
    treg_comparator.pl (http://treg.molgen.mpg.de)

=cut

use strict;
use Carp;
use File::Basename;
use File::Spec::Functions qw[catfile catdir];
use Getopt::Long;
use Time::Local;

# Locate Trawler modules
use FindBin ();
use lib "$FindBin::RealBin/../modules";
my $script_name = $FindBin::RealScript;
use Algorithm::Cluster;

# Internal Modules
use Trawler::Constants 1.0 qw(trawler_usage _read_config _tcst);
use Trawler::FileUtils 1.0 qw(_get_tstamp);
use Trawler::Utils 1.0 qw(_parse_DNA_file _reg_exp _reverse_complement _get_motif_loc _calculate_overlap);

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


# TRAWLER input ===============================================================
my $DNA_file = undef; # that 's the sample sequence file in FASTA format
my $background_DNA_file = undef; # that the background FASTA file for trawler
my $K_cutoff = 0; # SD from trawler default 0 (no filtering based on the standard deviation)
my $occurrence = 3; # minimum occurrence of the motif
my $WC = 2; # number of wild cards in motifs

# TRAWLER CLUSTERING input ====================================================
my $overlap_percentage = 70; # minimum overlap position of two motifs
my $seqlogo = 1; # if 1 it will create a png file for each cluster and a pwm file (FASTA like to use with t-reg comparator)
my $motif_number = 200;
my $EXTRA_SIDE_LEN = 0; # to build seq logos, we add $EXTRA_SIDE_LEN on each side
my $avoid; # motif to avoid ex : CA..TG
my $MOT_LENGTH = 0; # minimum motif length, used when parsing trawler
my $low_content = 0; # low occurrence are not filtered
my $nb_of_cluster = 0; # fixed number of cluster.
my $DIR = undef;
my $treg_threshold = 0.9;


#=================================================================

GetOptions(
  # Required
  'sample=s'          => \$DNA_file,
  'background=s'      => \$background_DNA_file,
  # motif discovery
  'occurrence:s'      => \$occurrence,
  'mlength:s'         => \$MOT_LENGTH,
  'wildcard:s'        => \$WC,
  # clustering
  'overlap:s'         => \$overlap_percentage, # default 70
  'motif_number:s'    => \$motif_number,
  'nb_of_cluster:s'   => \$nb_of_cluster, # default undef.
  # visulaization
  'directory=s'       => \$DIR,
  'xtralen:s'         => \$EXTRA_SIDE_LEN, # default 0
  # excluded for now
  'seqlogo:s'         => \$seqlogo, #default 1
  'avoid:s'           => \$avoid, # default NULL
  'low_content:s'     => \$low_content,
  'treg_threshold:s'  => \$treg_threshold, # default 0.9
  'k:s'               => \$K_cutoff,
);

# Whack out some help if lost...
if(!defined $DNA_file || !defined $background_DNA_file) {
    print STDERR "===================\n";
    print STDERR "Error: Invalid command line arguments\n";
    trawler_usage();
}

# if no directory provided => exit [initialized in trawler.pl]
unless($DIR) {
    print STDERR "\nERROR: no directory to run\n";
    exit(1);
}
# extract tmp_xxx
my $tmp_dir_name = fileparse($DIR);
# default result directory is like $TRAWLER_HOME/tmp_YYYY-MM-DD_HHhmm:ss/result
my $tmp_result_dir = File::Spec->catdir( $DIR, $tcst{RES_DIR_NAME} );
my $html_imgs_dir = File::Spec->catdir( $DIR, $tcst{HTML_IMG} );

# .cluster file
my $cluster_file = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{CULSTER_FILE_EXT});
# .trawler file
my $trawler_file = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{TRAWLER_FILE_EXT});

if ($DEBUG) {
    print "Cluster file: $cluster_file\n";
    print "Trawler file: $trawler_file\n";
}

#==============================================
# Run
msg_pipeline();

open(CLUST, ">$cluster_file") or croak "Cannot open file $cluster_file: $!";

#==============================================
#now that everthing is read, start to run trawler with the 2 FASTA files

my ($DNA) = _parse_DNA_file($DNA_file);

run_trawler($background_DNA_file, $DNA_file, $trawler_file, $occurrence, $WC); #background and sample file and output file


#==============================================
#print all the parrameters used in case you need them and you forgot to write them down
#print CLUST "#sequence : $DNA_file \t trawler : $trawler_file \t K cutoff = $K_cutoff\t minimum overlap = $overlap_percentage directory $DIR with name $DIR_NAME xtra length around the motif $EXTRA_SIDE_LEN avoid -$avoid- this motifs the number of motifs clustered are $motif_number the low_content is set at $low_content the occurrence is more than $occurrence and the minimum length is $MOT_LENGTH\n";


#parsing ==========================
#parsing DNA is already done
my $motif2info = parse_motif_file($trawler_file, $K_cutoff, $avoid, $MOT_LENGTH,  $motif_number, $low_content);
my %reference = %$motif2info; #used to get the info

#==================================

#=========================================================================
#first go throught the motif and whatever is the same cluster already but keep the most abondant
#=========================================================================
my @array_motif_tmp;
foreach my $motif (keys %$motif2info) {
    push @array_motif_tmp, $motif;
    print "$motif\n" if $DEBUG;
}

my $c = @array_motif_tmp;
my %deletedmots;
for (my $i=0; $i<$c; $i++) {
    my $mot1 = _reg_exp($array_motif_tmp[$i]); # regular expression fwd motif
    my $motif1 = $array_motif_tmp[$i];
    next if (exists $deletedmots{$motif1});

    print "$i motif out of $c done\n" if $DEBUG;
    for (my $j=$i+1; $j<$c; $j++) {
        my $mot2 = _reg_exp($array_motif_tmp[$j]); # regular expression fwd motif
        my $motif2 = $array_motif_tmp[$j];
        next if (exists $deletedmots{$motif2});

        if ($mot1 =~ /$mot2/ && $mot2 =~ /$mot1/) {
            my $k1 = $reference{$motif1}->{"K"};
            my $k2 = $reference{$motif2}->{"K"};
            if ($k1 >= $k2) {
                delete $$motif2info{$motif2};
                $deletedmots{$motif2}++;
            }
            else {
                delete $$motif2info{$motif1};
                $deletedmots{$motif1}++;
            }
            print "$motif1 est $motif2 $motif2 is deleted\n" if $DEBUG;
        }
    }
}

#CHARLES: compute match location only once
my %motif2location;
foreach my $motif (keys %$motif2info) {
    my $rv_motif = _reverse_complement($motif);
    my $mot = _reg_exp($motif); #regular expression fwd motif
    my $mot_rv = _reg_exp($rv_motif);
    my ($motifs, $starts) =  _get_motif_loc($DNA, $mot);
    my ($motifsrv, $startsrv) = _get_motif_loc($DNA, $mot_rv);
    $motif2location{$motif}->{"loc"}->{"1"} = $starts;
    $motif2location{$motif}->{"mot"}->{"1"} = $motifs;
    $motif2location{$motif}->{"loc"}->{"-1"} = $startsrv;
    $motif2location{$motif}->{"mot"}->{"-1"} = $motifsrv;
}


my $count_nb_cluster = 0;
my $pwm_file;

my @array_motif;
my $count = 0;


#=========================================================================
#second got through all the motifs found and locate them in the sequences
#and create the graph with only the nodes
#=========================================================================
#

foreach my $motif (keys %$motif2info) {
    push @array_motif, $motif;
    $count ++;
}

#=========================================================================
#then pairewise compare the position
#=========================================================================
#
# print STDERR "remains $count motif\n";
my @result;
my %som;
my $gr;

for (my $i=0; $i<$count; $i++) {

    $$gr[$i]=$i;
    #testing
     #$$gr[$i+$count]=$i+$count;

    my $motif1 = $array_motif[$i];
    my $location1 = $motif2location{$motif1}->{"loc"}->{"1"};
    my $location1bis = $motif2location{$motif1}->{"loc"}->{"-1"};
    print "done $i motif out of $count\n" if $DEBUG;
    for (my $j=$i+1; $j<$count; $j++) {

        my $motif2 = $array_motif[$j];
        if ($motif1 && $motif2) {

            my $location2 = $motif2location{$motif2}->{"loc"}->{"1"};
            my $location2bis = $motif2location{$motif2}->{"loc"}->{"-1"};

            my $overlap_A = _calculate_overlap($location1, $location2, $motif1, $motif2);
            my $overlap_B = _calculate_overlap($location1bis, $location2, $motif1, $motif2);
            my $overlap_C = _calculate_overlap($location1, $location2bis, $motif1, $motif2);

            my $som_overlap = 0;
            if ($overlap_A > $overlap_B && $overlap_A > $overlap_C) {
                $som_overlap = $overlap_A;
                print "arrive at possibility 1\n" if $DEBUG;
            }
            elsif ($overlap_B > $overlap_C) {
                $som_overlap = $overlap_B;
                print "arrive at possibility 2\n" if $DEBUG;
            }
            else {
                $som_overlap =$overlap_C;
                print "arrive at possibility 3\n" if $DEBUG;
            }

            $som{$i}->{$j}=$som_overlap;
            $som{$j}->{$i}=$som_overlap;

            #testing this give an outgroup to the som algorithm !
            #$som{$i}->{$j+$count}=0;
            #$som{$i+$count}->{$j}=0;
            #$som{$i+$count}->{$j+$count}=100;
            #$som{$j+$count}->{$i+$count}=100;

        }
    }
}


#=========================================================================
#if the round is the first one, then look at all the subgraph and clasify
#the motifs according to the subgraph number
#=========================================================================


#=========================================================================
#go over the 2D array and cluster using SOM
#=========================================================================

my $nxgrid = @$gr;
my $nygrid = $nxgrid;

my $cluster2motif_cluster;

if ($nb_of_cluster eq "0"){
    print STDERR "performing SCC clustering\n";
    $cluster2motif_cluster = get_SCC($gr, \%som, \@array_motif, $nxgrid, $nygrid); #SCC = strongly connected componment from http://blog.urbanomic.com/robin2/
}

elsif ($nb_of_cluster eq "som") {
    print STDERR "performing SOM clustering\n";
    $cluster2motif_cluster = get_som_cluster($gr, \%som, \@array_motif, $nxgrid, $nygrid);
}
elsif ($nb_of_cluster > 0) {
    print STDERR "performing kmean clustering\n";
    $cluster2motif_cluster = get_k_mean($gr, \%som, \@array_motif, $nxgrid, $nygrid, $nb_of_cluster);
}
else {
    die "please the number of cluster parameter should be an integer above 0.\n";
}


foreach my $cluster (keys %$cluster2motif_cluster) {
    my $motif_cluster = $$cluster2motif_cluster{$cluster};

    #now align the motifs of one cluster together resepctive to the master motif
    my $master_motif = shift @$motif_cluster; #remove the first element of the array as it is the master motif
    #get the number of the connected components from the first round
    #my $cc_number =  $reference{$master_motif}->{"connected_component"};
    my $best_motif;
    my $highest_score = 0;

    if (@$motif_cluster) { #if the cluster is not only one motif then align
        $count_nb_cluster ++;
        my $clus = "family_".$count_nb_cluster;
        my $aligned_motifs = align_motifs_from_cluster($motif_cluster, $master_motif);

        $pwm_file = return_PWM_from_cluster($DNA, $aligned_motifs, $clus);
        print "dealing with family $clus\n" if $DEBUG;
        foreach my $mot (@$aligned_motifs) {
            my $tmp_mot = $mot;
            $tmp_mot =~ s/\.//g;
            my $k = $reference{$tmp_mot}->{"K"};
            #check which motif has the higest score
            if ($k > $highest_score) { $highest_score = $k; $best_motif = $mot; }
            my $reverse_tmp_mot = _reverse_complement($tmp_mot);
            my ($occ_real, $occ_bkg);
            if ($reference{$tmp_mot}->{"occurrence"}) {
                $occ_real = $reference{$tmp_mot}->{"occurrence"};
                $occ_bkg = $reference{$tmp_mot}->{"occurrence_bkg"};
                print CLUST "$mot\t$clus\t$k\t$occ_real\t$occ_bkg\t1\n";
            }
            elsif ($reference{$reverse_tmp_mot}->{"occurrence"}) {
                $occ_real = $reference{$reverse_tmp_mot}->{"occurrence"};
                $occ_bkg = $reference{$reverse_tmp_mot}->{"occurrence_bkg"};
                my $k = $reference{$reverse_tmp_mot}->{"K"};
                print CLUST "$mot\t$clus\t$k\t$occ_real\t$occ_bkg\t-1\n";
            }
        }
    }
}
close(CLUST) or croak "Can't close file '$cluster_file': $!";


# _compare file
my $compare_file = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{COMPARE_FILE_EXT});

my $command_TREG1 = "\"" . $tcst{TREG_COMP} . "\" -q \"$pwm_file\" -f raw -d \"" . $tcst{TREG_MATRIX} . "\" -ori both -c $treg_threshold -out \"$compare_file\"";
print "treg cmd: $command_TREG1\n" if $INFO;
system("perl $command_TREG1") == 0 or croak "FAILURE command $command_TREG1:\n $!\n with return code $?\n";

print "arrived at the end !!\n\n" if $DEBUG;


#==============================================================================
# Subroutines
#==============================================================================

sub run_trawler {

    my ($background, $sample, $outputfile, $occurrence, $WC) = @_;

    #my $command = "java -Xmx500m -cp pecan.jar bp.trawler.Trawler -E $sample -F $background -N $outputfile -J $occurrence -G ACGTMRWSYKN -I $WC";
    my $command = "java " . $tcst{JAVA_OPTS} . " -cp \"" . $tcst{LIB_PECAN} . "\" bp.trawler.Trawler -E \"$sample\" -F \"$background\" -N \"$outputfile\" -J $occurrence -G ACGTMRWSYKN -I $WC";
    print "trawler cmd: $command\n" if $INFO;
    print CLUST "$command\n";
    #run the commande now
    system($command) == 0 or croak "FAILURE command $command:\n $!\n with return code $?\n";
}

sub return_PWM_from_cluster {

    my ($DNA, $array_motif, $clus) = @_;
    my %tmp;my %tmp_rc; my $longer_motif = 0;
    #create all the files :
    #=======================================
    my $file_weblogo = File::Spec->catfile($html_imgs_dir, $clus . $tcst{MOTIF_PNG_EXT});
    my $tmp_out = File::Spec->catfile($tmp_result_dir, "pwm.fasta");
    my $pwm_file = File::Spec->catfile($tmp_result_dir, $tmp_dir_name . $tcst{PWM_FILE_EXT});
    open(OUT_TMP, ">$tmp_out") or croak "Cannot open file $tmp_out: $!";
    open(OUT_PWM, ">>$pwm_file") or croak "Cannot open file $pwm_file: $!";
    #=======================================

    foreach my $motif (@$array_motif) {
        #found the longer motif :
        my $motif_length = length($motif);
        if ($motif_length > $longer_motif) {
            $longer_motif = $motif_length;
        }

        #found the motif in the sequence
        my $reg_exp = _reg_exp($motif);
        my ($motifs, $starts) = _get_motif_loc($DNA, $reg_exp);
        foreach my $start (@$starts) {
            $tmp{$start}++;
        }
        #CHARLES: found the rc_motif in the sequence
        my $rv_motif = _reverse_complement($motif);

        my $reg_exp_rc = _reg_exp($rv_motif);
        my ($rc_motifs, $rc_starts) = _get_motif_loc($DNA, $reg_exp_rc);
        foreach my $rc_start (@$rc_starts) {
            $tmp_rc{$rc_start}++;
        }

    }
    #go throught all the location and get the sites.
    my $size_motif = $longer_motif + (2*$EXTRA_SIDE_LEN);
    my $As = create_array($size_motif);
    my $Ts = create_array($size_motif);
    my $Cs = create_array($size_motif);
    my $Gs = create_array($size_motif);

    foreach my $location (sort {$a<=>$b} keys %tmp) {
        my $motif = substr($DNA, $location-$EXTRA_SIDE_LEN, $size_motif);
        print "$motif\n" if $DEBUG;
        #it might happen that we get back something like atXgMOTIFctXc => X in the extra length, we need to check for this to be correct
        #note : it can t happen to a X in each side so we suppose that only one X could in there
        if($motif=~/^([^X]*)X([^X]*)$/) {
            my $prefix = $1;
            my $suffix = $2;
            #which one should be canceled?
            if(length($prefix)<=length($suffix)) {
                #the prefix contains the X
                $prefix = "X" x length($prefix);
            }
            else {
                $suffix = "X" x length($suffix);
            }
            $motif = $prefix."X".$suffix;
        }

        print OUT_TMP ">$location\n$motif\n";

        #==============create the pwm===================
        for(my $o=0; $o<$size_motif; $o++) {
            my $nt = substr($motif, $o,1);
            if ($nt eq "A") { $$As[$o]++; }
            elsif ($nt eq "T") { $$Ts[$o]++; }
            elsif ($nt eq "C") { $$Cs[$o]++; }
            elsif ($nt eq "G") { $$Gs[$o]++; }
        }
        #=======finish creating the pwm==================

    }
    #CHARLES: get footprints from Reverse strand as well
    foreach my $location (sort {$a<=>$b} keys %tmp_rc) {
        my $motif_rc = substr($DNA, $location-$EXTRA_SIDE_LEN, $size_motif);
        my $motif = _reverse_complement($motif_rc);
        print "$motif\n$motif_rc (rc)\n" if $DEBUG;
        #it might happen that we get back something like atXgMOTIFctXc => X in the extra length, we need to check for this to be correct
        #note : it can t happen to a X in each side so we suppose that only one X could in there
        if($motif=~/^([^X]*)X([^X]*)$/) {
            my $prefix = $1;
            my $suffix = $2;
            #which one should be canceled?
            if(length($prefix)<=length($suffix)) {
                #the prefix contains the X
                $prefix = "X" x length($prefix);
            }
            else {
                $suffix = "X" x length($suffix);
            }
            $motif = $prefix."X".$suffix;
        }

        print OUT_TMP ">$location\n$motif\n";

        #==============create the pwm===================
        for(my $o=0; $o<$size_motif; $o++) {
            my $nt = substr($motif, $o,1);
            if ($nt eq "A") { $$As[$o]++; }
            elsif ($nt eq "T") { $$Ts[$o]++; }
            elsif ($nt eq "C") { $$Cs[$o]++; }
            elsif ($nt eq "G") { $$Gs[$o]++; }
        }
        #=======finish creating the pwm==================

    }
    #go through the pwm and write the OUT_PWM:
    print OUT_PWM ">$clus\n";
    for(my $o=0; $o<$size_motif; $o++) {
        my $a = $$As[$o]; my $t = $$Ts[$o]; my $c = $$Cs[$o]; my $g = $$Gs[$o];

        print OUT_PWM "$a\t$c\t$g\t$t\n";
    }

    close(OUT_TMP) or croak "Can't close file '$tmp_out': $!";
    close(OUT_PWM) or croak "Can't close file '$pwm_file': $!";

    # Run seqlogo
    my $commande = "perl \"" . $tcst{SEQLOGO} . "\" -F PNG -f \"$tmp_out\" -c > \"$file_weblogo\"";
    print "seqlogo cmd: $commande\n" if $INFO;
    system($commande) == 0 or croak "FAILURE command $commande:\n $!\n with return code $?\n";
    unlink $tmp_out;

    return ($pwm_file);
}

# TODO[YH]: method not used ! should be removed ??
#sub return_PWM { # do not use
#    my ($DNA, $motif, $clus) = @_;
#
#    my $tmp_out = "pwm.fasta";
#    open(OUT_PWM, ">$tmp_out") or croak "Cannot open file $tmp_out: $!";
#    my $reg_exp = _reg_exp($motif);
#    print STDERR "exp :$reg_exp\n";
#    my ($motifs, $starts) =  get_motif_loc($DNA, $reg_exp);
#    my $size_motif = length($motif);
#    #my @A= create_array($size_motif);  my @T = create_array($size_motif);  my @C  create_array($size_motif);  my @G  create_array($size_motif);
#    my $count = 0;
#    foreach my $mot (@$motifs) {
#        $count++;
#        print OUT_PWM ">$count\n$mot\n";
#   #for(my $i=0; $i<$size_motif; $i++)
#   #{
#   #    my $nt =substr($mot, $i,1);
#   #    TROUVE :
#   #    {
#   #   if ($nt eq "A"){$A[$i]++; last TROUVE;}
#   #   if ($nt eq "T"){$T[$i]++; last TROUVE;}
#   #   if ($nt eq "C"){$C[$i]++; last TROUVE;}
#   #   if ($nt eq "G"){$G[$i]++; last TROUVE;}
#   #    }
#   #}
#    }
#    close OUT_PWM;
#    my $file_weblogo = $clus."_".$motif.".png";
#    my $commande = $tcst{SEQLOGO} . " -F PNG -f $tmp_out -c > $file_weblogo";
#    system($commande) == 0 or croak "FAILURE command $commande:\n $!\n with return code $?\n";
#    unlink $tmp_out;
#
#    #return (\@A, \@C, \@T, @\G);
#}

sub align_motifs_from_cluster {

    my ($array, $master_motif) = @_;

    #go through all the motif in array and compare to master motif
    my $count_motif=0; my @array_result; my @array_result_stop;
    foreach my $mot (@$array) {
        my $motrv = _reverse_complement($mot);
        my ($position_mot, $best_score) = slide_motifs($mot, $master_motif);
        my ($position_motrv, $best_scorerv) = slide_motifs($motrv, $master_motif);

        if ($best_scorerv>$best_score) { #the alignment with the reverse complement is better
            $position_mot = $position_motrv;
            $$array[$count_motif] = $motrv;
        }
        $array_result[$count_motif] = $position_mot;
        $array_result_stop[$count_motif] = $position_mot + length($mot) -1;

        $count_motif++;
    }
    #create the coordinate for all the motifs in the cluster
    my $lower_coordinate = 0;
    my $higher_coordinate = length($master_motif)-1;
    foreach my $pos (@array_result) {

        if ($pos<$lower_coordinate){$lower_coordinate = $pos;}
        #if ($pos>$higher_coordinate){$higher_coordinate = $pos;}
    }
    foreach my $pos (@array_result_stop) {
        print "$pos $higher_coordinate\n" if $DEBUG;
        if ($pos>$higher_coordinate){$higher_coordinate = $pos;}
    }
    $lower_coordinate = abs($lower_coordinate);#absolute value again
    #put the master_motif in the appropriate coordinates
    my @return_result;
    my $dots1 = "."x $lower_coordinate;

    #transform higher
    $higher_coordinate += $lower_coordinate;
    print "$higher_coordinate <---\n" if $DEBUG;
    my $new_master = $dots1.$master_motif;

    if(length($new_master)<=$higher_coordinate) {
        my $dots2 = "."x ($higher_coordinate - length($new_master) +1);
        $new_master = $new_master.$dots2;
    }
    push @return_result, $new_master;
    #put all the other motifs in the appropriate coordinates
    for (my $i = 0; $i<$count_motif; $i++) {
        my $relative_pos = $array_result[$i];
        my $motif = $$array[$i];
        my $dots1 = "."x($lower_coordinate+$relative_pos);

        my $new_motif = $dots1.$motif;
        if(length($new_motif)<=$higher_coordinate) {
            my $dots2 = "."x ($higher_coordinate - length($new_motif) +1);
            $new_motif = $new_motif.$dots2;
        }

        push @return_result, $new_motif;

    }
    return \@return_result;
}

sub slide_motifs {

    my ($mot, $master_motif) = @_;

    my @mot_tmp = split//, $mot;
    my $size_mot = @mot_tmp;
    my $count = 1;

    my @master_tmp = split//, $master_motif;
    my $size_master = @master_tmp;

    my $best_score = 0;
    my ($position_master, $position_mot);
    #first slide mot against master
    for (my $i=$size_mot-1; $i>=0; $i--) {
        my $sub_mot = substr($mot, $i, $count);
        my $score = compare_clustered_motif($sub_mot, $master_motif);
        if ($score >$best_score) {$best_score = $score; $position_mot = -$i;}
        $count ++;
    }
    #then slice master against mot
    $count = $size_master;
    for (my $j=1; $j<$size_master; $j++) {
        my $sub_master =  substr($master_motif, $j, $count);
        my $score = compare_clustered_motif($sub_master, $mot);
        if ($score >$best_score) {$best_score = $score; $position_mot = $j;}
        $count--;
    }
    return ($position_mot, $best_score);
}

sub compare_clustered_motif {

    my ($mot1, $mot2) = @_;

    #this methods calculate the score for this specific alignment
    my @tmp1 = split//, $mot1;
    my @tmp2 = split//, $mot2;
    my $size1 = @tmp1; my $size2 = @tmp2; my $score =0; my $size;
    if ($size1>=$size2) {
        $size = $size1;
    }
    else {
        $size =$size2;
    }
    for (my $i=0; $i<$size; $i++) {
        my $nt1 = _reg_exp($tmp1[$i]);
        my $nt2 = _reg_exp($tmp2[$i]);

        if (($nt1&&$nt2) && ($nt1 =~ /$nt2/ || $nt2 =~/$nt1/)) {
            $score++;
        }
    }

    return $score;
}

sub get_SCC {
    my ($array_subgraph, $som, $array_motif, $nyxgrid, $nygrid, $nb_of_cluster) = @_;
    my %graph; my @nodes; my %result;
    
    foreach my $vertex (@$array_subgraph) {
      
	push @nodes, $vertex;
    }
    my $size_node = @nodes;
    my @distances; my @motifs; my @mask;my @weight;
    
    #===============================                                                                                
    #now construct the 2D array.                                                                                    
    #===============================                                                                                
    for (my $i = 0; $i<$size_node; $i++) {
        my $node1 = $nodes[$i];
        push @motifs, $node1;

        for (my $j = 0; $j<$size_node; $j++) {
            my $node2 = $nodes[$j];
	    if ($node1 ne $node2)
	    {
		my $d = $$som{$node1}->{$node2};
		if ($d > $overlap_percentage)
		{
		    push @{$graph{$node1}}, $node2;
		}
	    }
	}
    }
    my $clusters = strong_components(\%graph);
    foreach my $clus (keys %$clusters)
    {
       
	my $cluster = "cluster_".$clus;
	my @array = @{$$clusters{$clus}};
	foreach my $i (@array)
	{
	    my $motif = $$array_motif[$i];
	    if ($motif) {
		push @{$result{$cluster}}, $motif;
	    }
	}
    }
    return (\%result)
}

sub strong_components
{
    # Graph represented as a hash mapping vertices
    # to arrays of directly reachable vertices
    my ($graph1) = @_;
    
    my %times = ();
    my %components = ();
    my $current_time = 0;
    my $current_component = 0;
    my %roots = ();
    my @stack = ();
    my %result;
    foreach my $v (keys %$graph1)
    {
        if (!defined($times{$v}))
        {
            visit($v);
        }
    }
    
    sub visit
    {
	my ($v) = @_;
       
	$times{$v} = ++$current_time;
	$roots{$v} = $v;
	push @stack, $v;
	
	foreach my $w (@{$graph1->{$v}})
	{
	    if (!defined($times{$w}))
	    {
		visit($w);
	    }
	    
	    if (!defined($components{$w}) &&
		$times{$roots{$w}} < $times{$roots{$v}})
	    {
		$roots{$v} = $roots{$w};
	    }
	}
	
	if ($roots{$v} eq $v)
	{
	
	    
	    my $w = undef;
	    while ($w ne $v)
	    {
		$w = pop @stack;
	      
		push @{$result{$current_component}},$w;
		$components{$w} = $current_component;
	    }
	    
	    $current_component++;
	}
     
    }
return \%result;
}



sub get_k_mean {

    my ($array_subgraph, $som, $array_motif, $nyxgrid, $nygrid, $nb_of_cluster) = @_;

    my @nodes;
    #==================================================
    #first get all the nodes in the subgraph analysed
    #==================================================

    foreach my $vertex (@$array_subgraph) {
        push @nodes, $vertex;
    }
    my $size_node = @nodes;
    my @distances; my @motifs; my @mask;my @weight;


    #===============================
    #now construct the 2D array.
    #===============================
    for (my $i = 0; $i<$size_node; $i++) {
        my $node1 = $nodes[$i];
        push @motifs, $node1;

        for (my $j = 0; $j<$size_node; $j++) {
            my $node2 = $nodes[$j];
            if ($i == $j) {
                $distances[$i][$j] = 100;
            }
            else {
                my $d = $$som{$node1}->{$node2};
                if ($d >$overlap_percentage) {
                    $distances[$i][$j]=$d;
                }
                else {
                    $distances[$i][$j]=0;
                }
            }
        }
    }

    my $size = @distances;
    my @mask; my @weight;
    for (my $t=0; $t<$size; $t++) {
        $mask[$t] = 1;
        $weight[$t] = 1.0;
    }

print "the number of clsuter is equal to : $nb_of_cluster\n";
#------------------
# Define the params we want to pass to kcluster
    my %params = (
      nclusters => $nb_of_cluster,
      transpose => 0,
      npass     => 100,
      method    => 'a',
      dist      => 'e', #the clsutering algo to use DEFAULT = e (eukledian distance) but there is also 'a' = absolute value of correlation, 'u' uncentred correlation 'x' absolute uncetred correlation 's' spearman's rank correlation 'k' Kendall's Tau 'b' city block
      data      => \@distances,
      mask      => \@mask,
      weight    => \@weight,
    );
    print "start kmean\n" if $DEBUG;
    my ($clusters, $error, $found) = Algorithm::Cluster::kcluster(%params);
    print "end kmean\n" if $DEBUG;


    my $i=0;
    my %result;
    foreach my $cluster (@{$clusters}) {
        $cluster = "clsuter_".$cluster;
        my $motif_nb = $motifs[$i];
        my $motif = $$array_motif[$motif_nb];
        print "==>$cluster $motif\n" if $DEBUG;
        $i++;
        if ($motif) {
            push @{$result{$cluster}}, $motif;
        }
    }
    return \%result;

}

sub get_som_cluster {
    my ($array_subgraph, $som, $array_motif, $nyxgrid, $nygrid) = @_;
    my @nodes;
    #==================================================
    #first get all the nodes in the subgraph analysed
    #==================================================

    foreach my $vertex (@$array_subgraph) {
        push @nodes, $vertex;
    }
    my $size_node = @nodes;
    my @distances; my @motifs; my @mask;my @weight;


    #===============================
    #now construct the 2D array.
    #===============================
    for (my $i = 0; $i<$size_node; $i++) {
        my $node1 = $nodes[$i];
        push @motifs, $node1;

        for (my $j = 0; $j<$size_node; $j++) {
            my $node2 = $nodes[$j];
            if ($i == $j) {
                $distances[$i][$j] = 100;
            }
            else {
                my $d = $$som{$node1}->{$node2};
                if ($d >$overlap_percentage) {
                    $distances[$i][$j]=$d;
                }
                else {
                    $distances[$i][$j]=0;
                }
            }
        }
    }

    #print STDERR "som :\n";
    #foreach my $e (@distances) { print STDERR "@$e\n"; }
    #print STDERR "\n\n";

    #===========
    #run som
    #===========
    my %params = (
                  transpose => 0,
                  dist      => 'e',  #the clsutering algo to use DEFAULT = e (eukledian distance) but there is also 'a' = absolute value of correlation, 'u' uncentred correlation 'x' absolute uncetred correlation 's' spearman's rank correlation 'k' Kendall's Tau 'b' city block
                  data      => \@distances,
                  mask      => '',
                  #weight    => '',
                  #inittau   => 0.02,
                  #nxgrid    =>  $nxgrid,
                  #nygrid    =>  $nygrid,
                  #niter     =>       100,
          );
    print "start SOM\n" if $DEBUG;
    my $clusterid = Algorithm::Cluster::somcluster(%params);
    print "end SOM\n" if $DEBUG;
    my $i = 0;
    my %result;
    foreach my $xy (@{$clusterid}) {
        my $x = $$xy[0];
        my $y = $$xy[1];
        my $cluster = $x.$y;
        print "som cluster name $cluster\n" if $DEBUG;
        my $motif_nb = $motifs[$i];
        my $motif = $$array_motif[$motif_nb];
        $i++;
        if ($motif) {
            push @{$result{$cluster}}, $motif;
        }
    }
    return \%result;
}

#=====================================
#
#routines that deals with graph theory
#
#=====================================

sub get_the_higest_connected_components_and_neighbor {

    my ($array_subgraph, $graph, $array_motif) = @_;

    my $highest_connection = -1;
    my $result;
    #get the vertex with the highest connection (if many, well take the first)
    foreach my $vertex (@$array_subgraph) {
        my @e = $graph->edges_at($vertex);
        my $size_connection = @e;
        if ($size_connection > $highest_connection) {
            $highest_connection = $size_connection;
            $result = $vertex;
        }
    }
    my $motif = $$array_motif[$result];
    my @e = $graph->edges_at($result);
    my @cluster_result;
    my @remaining;

    #add the highest connected motifs===========
    push @cluster_result, $motif;
    #===========================================

    foreach my $element (@e) {
        my $v1 = $$element[0]; my $v2 = $$element[1]; my $neighbor;
        if ($v1 ne $result) { $neighbor = $$array_motif[$v1];  push @cluster_result, $neighbor; }
        elsif ($v2 ne $result) { $neighbor = $$array_motif[$v2]; push @cluster_result, $neighbor; }
    }

    #now return the motif cluster (array of motif) with the highest connected motif being the first instance
    return (\@cluster_result);
}

sub create_array {

    my ($length) = @_;

    print "caL: $length\n" if $DEBUG;
    my @result;
    for (my $i=0; $i<$length; $i++) {
        $result[$i] = 0;
    }
    print "caR: @result\n" if $DEBUG;
    return \@result;
}

sub parse_motif_file {

    my ($trawler_file, $cutoff, $avoid, $MOT_LENGTH,  $motif_number, $low_content) = @_;

    my %result;
    open(FILE, $trawler_file) or croak "Cannot open file $trawler_file: $!";
    my $new_trawler_file = $trawler_file."_short";
    my %tmp;
    #first filter the motifs that are not good
    foreach my $line (<FILE>) {
        chomp $line;
        my @array = split/\s+/, $line;
        my $motif = $array[3];
        my $length_motif = length($motif);
        my $test_low_content = 0;
        if ($low_content) {
            $test_low_content = check_for_low_content($motif);
        }
        if ($array[0]>=1 && $length_motif>= $MOT_LENGTH && $array[2]>=$cutoff && $test_low_content ==0) { #get only motif of x pb
            if ($avoid) {
                my $avoid_rv = _reverse_complement($avoid);
                if ($motif !~/$avoid/ && $motif !~/$avoid_rv/) {
                    my $value = $array[2];
                    push @{$tmp{$value}}, $line;
                }
            }
            else {
                my $value = $array[2];
                push @{$tmp{$value}}, $line;
            }
        }
    }
    close(FILE) or croak "Can't close file '$trawler_file': $!";
    #then order the motifs according the the score (decsending)
    my @tmp;
    foreach my $k (sort {$b<=>$a} keys %tmp) {
        my @a = @{$tmp{$k}};
        foreach my $line (@a) {
            push @tmp, $line;
        }
    }
    #then get only the $number_motif with the best score only
    my $actual_nb_of_motif = @tmp;
    if ($actual_nb_of_motif < $motif_number) { $motif_number = $actual_nb_of_motif;}

    #if trawler has not find anything satisfying the motif criteria, then Trawler dies.
    if ($actual_nb_of_motif == 0) {
        die "================================\nNo motif satisfying your criteria has been found. Please re-run Trawler_standalone using different parameters.\n==================\n";
    }

    open(OUT_TRAWLER, ">$new_trawler_file") or croak "can't open the trawler output $new_trawler_file for making the sort file: $!";
    for (my $i=0; $i<$motif_number; $i++) {
        my $line = $tmp[$i];
        print OUT_TRAWLER "$line\n";
        my @array = split/\s+/, $line;
        my $motif = $array[3];
        my $K = $array[2];
        my $occurrence_bkg = $array[1];
        my $occurrence = $array[0];
        $result{$motif}->{"K"}=$K;
        $result{$motif}->{"occurrence"}=$occurrence;
        $result{$motif}->{"occurrence_bkg"}=$occurrence_bkg;
    }
    close(OUT_TRAWLER) or croak "Can't close trawler output '$new_trawler_file': $!";
    return \%result;
}

sub check_for_low_content {
    my ($motif) = @_;
    my $reg_exp_motif = _reg_exp($motif);
    $reg_exp_motif =~ s/\[//g;
    $reg_exp_motif =~ s/\]//g;
    my %hash;
    my @tmp = split //, $reg_exp_motif;
    foreach my $nt (@tmp) {
        $hash{$nt}++;
    }
    my $size = keys %hash;
    my $test;
    if ($size <=2 ) { $test = 1; }
    else { $test = 0; }
    return $test;
}

sub msg_pipeline {
  my $msg_pipeline = <<"MSG";
  ==========
  Running Trawler...
  ==========
MSG

  print $msg_pipeline . "\n";
}

1;
