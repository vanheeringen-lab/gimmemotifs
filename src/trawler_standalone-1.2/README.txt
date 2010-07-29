===============
Welcome to Trawler standalone 1.2
===============

The GNU General Public License (GPL)
     Authors: L. Ettwiller, B. Paten, M. Ramialison, Y. Haudry


Reference: Ettwiller L., Paten B.,Ramialison M., Birney E., Wittbrodt J.
     Trawler: de novo regulatory motif discovery pipeline for chromatin immunoprecipitation.
          Nat Methods. 2007 Jul;4(7):563-5. Epub 2007 Jun 24.
          PMID: 17589518


More information about Trawler can be found at:

     http://ani.embl.de/laurence/trawler/


*****************************************************************
***                                                           ***
***               Installation procedure                      ***
***                                                           ***
***                                                           ***
***                                                           ***
***       System Requirements                                 ***
***       Installation                                        ***
***       Dependencies                                        ***
***       Licensing                                           ***
***       Usage                                               ***
***                                                           ***
***    Documentation available through                        ***
***       http://ani.embl.de/trawler/doc                      ***
***                                                           ***
*****************************************************************

*********************************
***                           ***
***     System Requirements   ***
***                           ***
*********************************

The program has been tested successfully on Unix/Linux (Fedora 9),
Mac OS X (versions 10.4, 10.5) and Windows (XP).

The following software packages are required for running Trawler:
perl 5.6 or above
Java 5.0 or above

*********************************
***                           ***
***     Installation          ***
***                           ***
*********************************

please see the INSTALL.txt file for detailed installation instructions.


*********************************
***                           ***
***     Dependencies          ***
***                           ***
*********************************

The following software packages are required:
- Algorithm::Cluster (CPAN: http://search.cpan.org/~mdehoon/)

The following software packages are directly provided in the Trawler distribution:
- treg_comparator (http://treg.molgen.mpg.de)
- Jalview (http://www.jalview.org/)
- WebLogo (http://weblogo.berkeley.edu/)


*********************************
***                           ***
***     Licensing             ***
***                           ***
*********************************

Please see the file called LICENSE.txt


*********************************
***                           ***
***     USAGE                 ***
***                           ***
*********************************

USAGE
=====

trawler -sample [file containing the enriched sequences] -background [file containing the background sequences]

    -sample (FASTA format) better to be repeat-masked.
    -background (FASTA format)

OPTIONAL PARAMETERS
===================

    [MOTIF DISCOVERY]
    -occurrence (optional) minimum occurrence in the sample sequences. [DEFAULT = 10]
    -mlength (optional) minimum motif length. [DEFAULT = 6]
    -wildcard (optional) number of wild card in motifs. [DEFAULT = 2]
    -strand (optional) single or double. [DEFAULT = double]

    [CLUSTERING]
    -overlap (optional) in percentage. [DEFAULT = 70]
    -motif_number (optional) total number of motifs to be clustered. [DEFAULT = 200]
    -nb_of_cluster (optional) fixed number of cluster. if this option is used, 
the k-mean clustering algorithm with fixed k will be used instead of the strongly 
connected component (SCC). Alternatively this option can be set to 'som' and in 
this case, the self organizing map (SOM) will be use. [DEFAULT = NULL]


    [VISUALIZATION]
    -directory (optional) output directory. [DEFAULT = "$TRAWLER_HOME/myResults"]
    -dir_id (optional) gives an id to the results directory. [DEFAULT = NULL]
    -xtralen (optional) add bases around the motifs for the logo. [DEFAULT = 0]
    -alignments (optional) file containing the list of files containing the aligned sequences (see README file for more info) [DEFAULT = NULL]
    -ref_species (optional) name of the reference species [DEFAULT = NULL]
    -clustering (optional) if 1 the program clusters the instances, if 0 no clustering. [DEFAULT = 1]
    -web (optional) if 1 the output will be a web page with all the information [DEFAULT = 1]

USAGE EXAMPLES
==============

Trawler-standalone comes in two flavors: [1] with the conservation information or
[2] without. For the conservation you need to provide Trawler-standalone with a file
containing all the orthologous sequences (not aligned).

[1] with conservation

> trawler.pl -sample file_sample -background file_background -alignments file_alignments -ref_species reference_species_name

[2] without conservation

> trawler.pl -sample file_sample -background file_background


By default Trawler-standalone performs the full analysis nevertheless it is possible to

[1] get only the over-represented motifs:

> trawler.pl -sample file_sample -background file_bakground -clustering 0

[2] get only the over-represented motifs and the cluster (PWM):

> trawler.pl -sample file_sample -background file_bakground -web 0

OPTIONS DESCRIPTION
===================

    [MOTIF DISCOVERY]

    -occurrence : [DEFAULT = 10] The minimum number of time the motif occurs in the sample sequence.
            The lower limit defined by the suffix tree is 2 motifs. Nevertheless if your chromatin IP
            contains more than just a few sequences, you would expect to see the motif occurring quite often.
            For typical chromatin IP data, a good minimum number of motif occurrence is about 10 to 20.
            If you do not get anything interesting with this values, try a higher or lower number according
            to your sample size.
    -mlength : [DEFAULT = 6] Trawler-standalone does not assess motifs of fixed length.
            The maximum length has been set to 20 nucleotides but the minimum length is a user defined parrameter.
            A good minimum length is around 5-6 bp for typical transcription factor binding sites.
    -wildcard : [DEFAULT = 2] The number of wildcards is the maximum number of positions in the motif
            that can have a degenerate nucleotide. In IUPAC code, the possible wildcards used in Trawler-standalone are
            the following :
            M = [AC], R = [AG], Y = [CT], W = [AT], S = [CG], K = [GT], N = [ATCG]
            (the latest counts as 2 mismatches). A possible motif with -wildcard option set
            to 2 can be AGCMTWA for example.
    -strand : [DEFAULT = double] The analysis can be performed either on single-stranded or double-stranded sequences.

    [CLUSTERING]

    -overlap : [DEFAULT = 70] The minimum percentage overlap for two motifs to be considered in the same cluster.
            The value is by default 70 %; lowering down this value will cluster more distant motifs.
    -motif_number : [DEFAULT = 200] The maximum number of motifs to be considered for clustering.
            If the option is set to 200, the best 200 motifs returned by the discovery step will be used
            for the clustering step.
    -nb_of_cluster : [DEFAULT = NULL] By default, the number of cluster(s) is defined by the
            number of strongly connected components (SCC). By setting this option to a integer, 
	    the user set the number of expected cluster(s). In this case, the SCC is replaced 
	    by the k-mean clustering algorithm. if the option is set to 'som' the SOM algorithm is run
	    instead.

    [VISUALIZATION]

    -directory : [DEFAULT = "$TRAWLER_HOME/myResults"] name of the directory where all the outputs
            will be stored. If no directory is specified, then a default directory "myResults" will be created,
            Individual results will be stored in a separate directory with the following pattern
            "tmp_yyyy-mm-dd_hh:mm_ss" (ex : tmp_2009-05-20_12h13_20).
    -dir_id : [DEFAULT = NULL] gives a meaningful ID to the results directory. For instance, if this option
            is set to "myID", this will create a directory named "myID_yyyy-mm-dd_hh:mm_ss".
    -xtralen : [DEFAULT = 0] The number of additional nucleotides included in the PWM flanking the core motif.
    -alignments : [DEFAULT = NULL] Alignment file in the correct format. If alignments are provided,
            the motifs are visualized in the context of an alignment using Jalview.
            If the -alignment option is activated, the option -ref_species needs to be set.
    -ref_species : [DEFAULT = NULL] Name of the reference species. Only use this option if
            the -alignment option is activated.
    -clustering : [DEFAULT = 1]. By default, Trawler proceeds to the clustering step.
            To only obtain the over-represented motifs, set the -clustering option to 0.
            If the -clustering option is set to 0, no webpage is produced.
    -web : [DEFAULT = 1] By default, Trawler proceeds to create a web page summarizing the result.
            If this option is set to 0, Trawler proceeds until the clustering step only.
            Useful for batch  analysis.


*********************************
***                           ***
*** Additional features       ***
***                           ***
*********************************


Additional reference matrices can be added to the known matrices used by Treg comparator by adding
to the file matrix_set_all your PWM in the following format :

ID#name#text#text#0.3548828125 ...
Matrix entries (float) are separated by a whitespace. The matrix entries
start with the frequency of "A", followed by that of "C", "G", and "T" in the
first position (5'), then the same for the second position and so on.


*********************************
***                           ***
*** Trawler Team              ***
***                           ***
*********************************

Laurence Ettwiller  <ettwille@embl.de>
Benedict Paten      <benedict@soe.ucsc.edu>
Mirana Ramialison   <ramialis@embl.de>
Yannick Haudry      <haudry@embl.de>
