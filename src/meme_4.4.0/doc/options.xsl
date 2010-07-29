<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output 
    method="xml"
    omit-xml-declaration="yes"
    doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"
    doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd" 
  />

  <!-- 
    This templates in this stylesheet contain the text describing the command line
    options for all of the programs in the MEME Suite. This is done to standardize
    the descriptions across the programs, avoiding duplication and inconsitencies.

    Some options, particularly once represented by single letters, have different
    meanings in different applilcats ("type" for example). Use the program 
    attribute to distinguish them if needed (program='gendb' name='type' for  example).
  -->

  <xsl:template match="@*|node()">
    <xsl:copy>
      <xsl:apply-templates select="@*|node()"/>
    </xsl:copy>
  </xsl:template>

  <xsl:template match="command-option[@name='allow-weak-motifs']">
    <code>--allow-weak-motifs</code> 
    - In p-value score mode, weak motifs
    are defined as ones where the best possible hit has a p-value greater
    than the p-value threshold.  Such motifs cannot contribute to a match
    in p-value score mode. By default, the program rejects any search
    results containing weak motifs, unless the
    <code>--allow-weak-motifs</code> switch is given.  In that case, the
    search will proceed, but the weak motifs will never appear in any
    matches.  <b>Note</b>:This switch only applies to p-value score mode.
  </xsl:template>

  <xsl:template match="command-option[@name='bg']"> 
    <code>--bg <i>&lt;float&gt;</i></code>
    - The mutation rate for sites in the
    background model.  The default value is 1.
  </xsl:template>

  <xsl:template match="command-option[@name='bgfile']"> 
    <code>--bgfile <i>&lt;bfile&gt;</i></code> 
    - Read background frequencies from
    <code><i>&lt;bfile&gt;</i></code>.  
    The file should be in 
    <a href="bfile-format.html">MEME background file format</a>.
    The default is to use frequencies
    embedded in the application from the non-redundant database.  If the
    argument is the keyword <code>motif-file</code>, then the frequencies
    will be taken from the motif file.
  </xsl:template>

  <xsl:template match="command-option[@name='bgweight']"> 
    <code>--bgweight <i>&lt;weight&gt;</i></code> 
    - Add <i>&lt;weight&gt;</i> times the background frequency to
    the corresponding letter counts in each motif when
    converting them to postion specific scoring matrices.
    The default value is 4.0.
  </xsl:template>

  <xsl:template match="command-option[@name='blocksize']"> 
    <code>--blocksize <i>&lt;n&gt;</i></code> 
    - Read sequences in blocks of size <code><i>n</i></code>. 
    Default behavior is
    to read one full sequence at a time.  
  </xsl:template>

  <xsl:template match="command-option[@name='column-freqs']">
    <code>--column-freqs [simulated|empirical]</code>
    - The way to compute the frequencies of all possible columns.
    <ul>
      <li><code>simulated</code> - Use the evolutionary model to compute the frequency of each
      possible column of letters.</li>
      <li><code>empirical</code> - Count the numbers of each column in the input multiple
      alignments.  All alignments using the same (sub)set of species are 
      counted together.  Frequencies are computed by dividing each by the 
      total counts for that (sub)set of species.</li>
    </ul>
    The default is <code>simulated</code>.  These frequencies are used for
    computing <i>p</i>-values for scores.  If <code>simulated</code> is used, 
    the accuracy of the <i>p</i>-values depends strongly on the accuracy of the 
    evolutionary model.
  </xsl:template>

  <xsl:template match="command-option[@name='comp']">
    <code>--comp</code>
    - Adjust p-values and E-values for sequence composition.
  </xsl:template>

  <xsl:template match="command-option[@name='dna']">
    <code>-dna</code>
    - Sets alphabet to DNA (default protein).
  </xsl:template>

  <xsl:template match="command-option[@name='dummy']">
    <code>--dummy</code> 
    - Print a "dummy" sequence record before the generated sequences. 
    The "dummy" sequence record is a a FASTA header line listing the 
    <code>gendb</code> parameters but not followed by any sequence lines.
  </xsl:template>

  <xsl:template match="command-option[@name='eg-cost']">
    <code>--eg-cost <i>&lt;cost&gt;</i></code>
    - Scale the expected cost of a
    random gap to be <i><code>&lt;cost&gt;</code></i>
    times the expected score of a
    random hit. By default, gap costs are essentially zero. The larger
    you set <i><code>&lt;cost&gt;</code></i>, 
    the more gaps will be penalized. 
    This can only be used in conjunction with <code>--max-gap</code>. 
    This may not be used in conjunction with <code>--min-score</code>.
  </xsl:template>

  <xsl:template match="command-option[@name='e-thresh']">
    <code>--e-thresh <i>&lt;ev&gt;</i></code>
    - Only print results with E-values less than 
    <code><i>&lt;ev&gt;</i></code>.
  </xsl:template>

  <xsl:template match="command-option[@program='mcast' and @name='e-thresh']">
    <code>--e-thresh <i>&lt;ev&gt;</i></code>
    - Only print results with E-values less than 
    <code><i>&lt;ev&gt;</i></code>.
    The default threshold is 10.0.
  </xsl:template>

  <xsl:template match="command-option[@name='fg']">
    <code>--fg <i>&lt;float&gt;</i></code>
    - The mutation rate for sites in the foreground model(s).
    The default value is 1.
  </xsl:template>

  <xsl:template match="command-option[@name='fim']">
    <code>--fim</code>
    - Gaps between motifs are not penalized.
    Spacer states between motifs are represented as free-insertion modules (FIM).
    A FIM is an insert state with 1.0 probability of self-transition 
    and 1.0 probability of exit transition. 
    Thus, traversing such a state has zero transition cost.
    Specifying this option causes all spacers to be represented using FIMs.
  </xsl:template>

  <xsl:template match="command-option[@name='gap']">
    <code>--gap <i>&lt;method&gt;</i></code>
    - Specifies the gap handling strategy.
    Allowed values for method are:
    <ul>
      <li> <code>skip</code> Skip those sites where any position in the alignment
      window contains a gap. This is the default gap handling strategy.
      </li>
      <li> <code>fixed</code> Sites containing gaps are assigned a fixed
      score, specified by <code>--gap-cost</code>.
      </li>
      <li> <code>wildcard</code> The gap character matches any base, and the
      score is the product of the corresponding probabilities.
      </li>
      <li> <code>minimum</code> The gap character is assigned the score
      corresponding to the least likely letter at the given position.
      </li>
      <li> <code>model</code> Use model-specific gap handling. Currently, the only
      model that supports this is <code>f81_gap</code>.
      </li>
    </ul>
  </xsl:template>

  <xsl:template match="command-option[@name='gap-cost']">
    <code>--gap-cost <i>&lt;float&gt;</i></code>
    - Specifies the costs for gaps when
      using the <code>fixed</code> gap handling strategy. Default is 0.0.
  </xsl:template>

  <xsl:template match="command-option[@name='gap-extend']">
    <code>--gap-extend <i>&lt;cost&gt;</i></code>
    - This switch causes
    <b>all</b> spacer self-loop log-odds scores to be set to
    <code><i>&lt;cost&gt;</i></code>. 
    In addition, it causes all other transitions out of a
    spacer to be set to zero. Together with the <code>--gap-open</code>
    switch, this allows you to specify an affine gap penalty function,
    overriding the gap penalty implicit in the model (self-loop transition
    probabilities of gap states).
  </xsl:template>

  <xsl:template match="command-option[@name='gap-freq']">
    <code>--gap-freq &lt;float&gt;</code> 
    - Specifies the background frequency
    for gaps. Default is derived from the alignment.
  </xsl:template>

  <xsl:template match="command-option[@name='gap-open']">
    <code>--gap-open <i>&lt;cost&gt;</i></code>
    - This switch causes
    <b>all</b> transitions into a spacer state to be assigned a log-odds
    score equal to <code><i>&lt;cost&gt;</i></code>.
    Together with the <code>--gap-extend</code> switch,
    this allows you to specify an affine gap penalty function,
    overriding the gap penalty implicit in the model
    (transition probabilities into and out of gap states).
  </xsl:template>

  <xsl:template match="command-option[@name='global']">
    <code>--global</code>
    - Scores are computed for the match between the entire sequence 
    and the model (the default is to use the maximal local score).
  </xsl:template>

  <xsl:template match="command-option[@name='hb']">
    <code>--hb</code> 
    - Use the Halpern-Bruno modification to the evolutionary model.
  </xsl:template>

  <xsl:template match="command-option[@name='help']">
    <code>--help</code> 
    - Print a usage statement.
  </xsl:template>

  <xsl:template match="command-option[@name='indices']">
  <code>--indices</code>
  - Limits output to indices of motifs in the provided motif file.
  </xsl:template>

  <xsl:template match="command-option[@name='keep-unused']">
    <code>--keep-unused</code>
    - By default all inter-motif transitions that are
    not observed in the data are removed from the transition probability matrix.
    This option allows those transitions to be retained.
    This option is only relevant if the model has a
    completely connected topology.
  </xsl:template>

  <xsl:template match="command-option[@name='list']">
    <code>--list</code>
    - Treat the second required input as a list of
    alignments, rather than a single alignment.
  </xsl:template>

  <xsl:template match="command-option[@name='many']">
    <code>-many</code>
    - Read all sequences into memory at start up.
  </xsl:template>

  <xsl:template match="command-option[@name='max-gap']">
    <code>--max-gap <i>&lt;max-gap&gt;</i></code>
    - The value of <code><i>&lt;max-gap&gt;</i></code> specifies
    the longest distance allowed between two hits in a match.
    Hits separated by more than <code><i>&lt;max-gap&gt;</i></code> 
    will be placed in different matches. The default value is 50.
    <b>Note</b>:
    Large values of <code><i>&lt;max-gap&gt;</i></code> combined with
    large values of <i>pthresh</i> may prevent <code>MCAST</code> from
    computing <i>E</i>-values.
  </xsl:template>

  <xsl:template match="command-option[@name='max-seqs']">
    <code>--max-seqs <code><i>&lt;max&gt;</i></code></code>
    - Print results for no more than <code><i>&lt;max&gt;</i></code>
    sequences.  By default, all matches are reported, up to the specified
    E-value threshold (see <code>--e-thresh</code>).
  </xsl:template>

  <xsl:template match="command-option[@name='max-seq-length']">
    <code>--max-seq-length <code><i>&lt;max&gt;</i></code></code>
    - Set the maximum length allowed for input sequences.
    By default the maximum allowed length is 250000000.
  </xsl:template>

  <xsl:template match="command-option[@name='last']">
    <code>--last<code><i>&lt;max&gt;</i></code></code>
    -Use only scores of (up to) last &lt;n&gt; sequence
    positions to compute AMA. If the sequence is shorter
    than this value the entire sequence is scored. If the
    motif is longer than this value it will not be scored. 
  </xsl:template>

  <xsl:template match="command-option[@name='max-stored-scores']">
    <code>--max-stored-scores <code><i>&lt;max&gt;</i></code></code>
    - Set the maximum number of scores that will be stored.
    Precise calculation of q-values depends on having a complete list of scores.
    However, keeping a complete list of scores may exceed available memory.
    Once the number of stored scores reaches the maximum allowed, 
    the least significant 50% of scores will be dropped,
    and approximate q-values will be calculated.
    By default the maximum number of stored matches is 100,000.
  </xsl:template>

  <xsl:template match="command-option[@name='maxseq']">
    <code>--maxseq &lt;n&gt;</code>
    - Maximum sequence length. Default value: 2000.
  </xsl:template>

  <xsl:template match="command-option[@name='min-score']">
    <code>--min-score <i>&lt;minscore&gt;</i></code>
    - This switch allows you to
    specify the threshold for the repeated match algorithm used by
    <code>miao</code>.
    Matches must have a score of at least
    <code><i>&lt;minscore&gt;</i></code> to be detected.
    Matches containing internal regions with scores less than minus 
    'threshold' will be split and reported as two separate matches.
  </xsl:template>

  <xsl:template match="command-option[@name='minseq']">
    <code>--minseq &lt;n&gt;</code>
    - Minimum sequence length. Default value: 50.
  </xsl:template>

  <xsl:template match="command-option[@name='model']">
    <xsl:choose>
      <xsl:when test="@pgm='motiph'">
        <code>--model [single|average|jc|k2|f81|f84|hky|tn]</code> 
      </xsl:when>
      <xsl:otherwise>
        <code>--model [jc|k2|f81|f84|hky|tn]</code> 
      </xsl:otherwise>
    </xsl:choose>
    - The evolutionary model to use.
    The available models are:
    <ul>
    <xsl:choose>
      <xsl:when test="@pgm='motiph'">
        <li>single - score first sequence: compute standard log-odds score of
	  first sequence in the alignment; ignores tree but does NOT remove gaps.</li>
        <li>average - compute average of standard log-odds score of aligned sites.</li>
      </xsl:when>
    </xsl:choose>

    <li>jc - Jukes-Cantor: equilibrium base frequencies are all 1/4;
    the only free parameter is the mutation rate.</li>
    <li>k2 - Kimura 2-parameter: equilibrium base frequencies are all 1/4;
    the free parameters are the mutation rate and the transition/transversion rate ratio.</li>
    <li>f81 - Felsenstein 1981: equilibrium base frequencies are taken from the alignment;
    the only free parameter is the mutation rate.</li>
    <li>f84 - Felsenstein 1984: equilibrium base frequencies are taken
    from the alignment; the free parameters are the mutation rate and the
    transition/transversion rate ratio.  The ratio of purine-purine to
    pyrimidine->pyrimidine transitions is assumed to be 1.</li>
    <li>hky - Hasegawa-Kishino-Yano: equilibrium base frequencies are
    taken from the alignment; the free parameters are the mutation rate
    and the transition/transversion rate ratio.  The ratio of
    purine-purine to pyrimidine-pyrimidine transitions is assumed to be
    equal to the  ratio of purines to pyrimidines.</li>
    <li>tn - Tamura-Nei: equilibrium base frequencies are taken from the
    alignment;  the free parameters are the mutation rate, the
    transition/transversion rate ratio, and the ratio of purine-purine
    transitions to pyrimidine-pyrimidine transitions.</li>

    </ul>
    <p>
    The default model is f81.  A description of the f81 model is available
    in chapter 13 of <i>Statistical Methods in Bioinformatics</i> by Ewens
    and Grant. The other models are described in chapters 9 and 13 of
    <i>Inferring Phylogenies</i> by Felsenstein.
    </p>
  </xsl:template>

  <xsl:template match="command-option[@name='model-file']">
    <code>--model-file <i>&lt;model file&gt;</i></code>
    - Creation of the HMM will be skipped,
    and the HMM will be read from the file instead.
  </xsl:template>

  <xsl:template match="command-option[@name='motif']">
    <code>--motif <i>&lt;id&gt;</i></code>
    - Use only the motif identified by <code><i>&lt;id&gt;</i></code>. 
    This option may be repeated.
  </xsl:template>

  <xsl:template match="command-option[@name='motif-e-thresh']">
    <code>--motif-e-thresh <i>&lt;ev&gt;</i></code>
    - Only motifs with E-values less than 
    <code><i>&lt;ev&gt;</i></code> will be used to build the HMM.
  </xsl:template>

  <xsl:template match="command-option[@name='motif-p-thresh']">
    <code>--motif-p-thresh <i>&lt;pv&gt;</i></code>
    - Only motif occurences with p-values less than 
    <code><i>&lt;pv&gt;</i></code> will be used to build the HMM.
  </xsl:template>

  <xsl:template match="command-option[@name='motif-pseudo']">
    <code>--motif-pseudo <i>&lt;float&gt;</i></code>
    - A pseudocount to be added to each count in the motif matrix,
    after first multiplying by the corresponding background frequency (default=0.1).
  </xsl:template>

  <xsl:template match="command-option[@name='no-pvalue']">
    <code>--no-pvalue</code>
    - Skip the p-value calculation.  This 
    switch will be necessary when a large number <i>n</i> of species are
    in the tree, because the memory requirement is
    4<sup><i>n</i></sup>. This also disables computation of q-values.
  </xsl:template>

  <xsl:template match="command-option[@name='no-qvalue']">
    <code>--no-qvalue</code>
    - Do not compute a q-value for each p-value.  The q-value 
    calculation is that of Benjamini and Hochberg (1995).
    By default, q-values are computed.
  </xsl:template>

  <xsl:template match="command-option[@name='norc']">
    <code>--norc</code>
    - Do not score the reverse complement DNA strand.  Both
	strands are scored by default.
  </xsl:template>

  <xsl:template match="command-option[@name='no-search']">
    <code>--no-search</code>
    - This option turns off the search phase of <code>beadstring</code>.
    The HMM will be stored if the <code>--model</code> option is specified.
  </xsl:template>

  <xsl:template match="command-option[@name='nspacer']">
    <code>--nspacer <i>&lt;value&gt;</i></code>
    - By default each spacer is modeled using a single insert state.
    The distribution of spacer lengths produced by a single insert
    state is exponential in form.
    A more reasonable distribution would be a bell-shaped curve such as a Gaussian.
    Modeling the length distribution explicitly is computationally expensive;
    however, a Gaussian distribution can be
    approximated using multiple insert states to represent a single
    spacer region. The <code>--nspacer</code> option specifies the number of
    insert states used to represent each spacer.
  </xsl:template>

  <xsl:template match="command-option[@name='numseqs']">
    <code>--numseqs &lt;n&gt;</code>
    - Number of sequences to generate. Default value: 10.
  </xsl:template>

  <xsl:template match="command-option[@name='o']">
    <code>--o <i>&lt;dir name&gt;</i></code>
    - Specifies the output directory.
    If the directory already exists, the contents will not be overwritten.
  </xsl:template>

  <xsl:template match="command-option[@name='oc']">
    <code>--oc <i>&lt;dir name&gt;</i></code>
    - Specifies the output directory.
    If the directory already exists, the contents will be overwritten.
  </xsl:template>

  <xsl:template match="command-option[@name='order']">
    <code>--order <i>&lt;string&gt;</i></code>
    - The given string specifies the order and spacing of the motifs within the
    model, and has the format "l=n=l=n=...=l=n=l", where "l" is the
    length of a region between motifs, and "n" is a motif index. Thus,
    for example, the string "34=3=17=2=5" specifies a two-motif linear
    model, with motifs 3 and 2 separated by 17 letters and flanked by
    34 letters and 5 letters on the left and right. If the motif file
    contains motif occurrences on both strands, then the motif IDs in
    the order string should be preceded by "+" or "-" indicating the
    strandedness of the motif.
  </xsl:template>

  <xsl:template match="command-option[@name='output-pthresh']">
    <code>--output-pthresh <i>&lt;float&gt;</i></code>
    - The p-value threshold for displaying search results.
    If the p-value of a match is greater than this value,
    then the match will not be printed. 
    Using the <code>--output-pthresh</code> 
    option will set the q-value threshold to 1.0.
    The default p-value threshold is 1e-4.
  </xsl:template>

  <xsl:template match="command-option[@name='output-qthresh']">
    <code>--output-qthresh <i>&lt;float&gt;</i></code>
    - The q-value threshold for displaying search results.
    If the q-value of a match is greater than this value,
    then the match will not be printed.
    Using the <code>--output-qthresh</code> 
    option will set the p-value threshold to 1.0.
    The default q-value threshold is 1.0.
  </xsl:template>

  <xsl:template match="command-option[@name='p-score']">
    <code>--p-score <i>&lt;float&gt;</i></code>
    - The <code>--p-score</code> switch
    activates p-value score mode with the given threshold. (The default
    score mode is called "log-odds score mode".) In p-value score mode,
    motif match scores are converted to their p-values. They are then
    converted to bit scores as follows:
    <div style="margin-left: 2em">S = -log<sub>2</sub>(p/T)</div>
    where S is the bit score of the hit,
    p is the p-value of the log-odds score,
    and T is the p-value threshold.
    In this way,
    only hits more significant than the p-value threshold get positive scores.
    The p-value threshold, T, must be in the range 0&lt;T&lt;=1.
  </xsl:template>

  <xsl:template match="command-option[@name='pam']">
    <code>--pam <i>&lt;distance&gt;</i></code>
    - By default, target probabilities are derived
    from the <em>distance-250</em> PAM matrix for proteins, and from a
    <code><i>&lt;distance&gt;-1</i></code>
    transition/transversion matrix for DNA.  
    With the <code>-pam</code> switch, 
    you can specify a different <b>integer</b> distance
    from 1 to 500. (This can be overridden with 
    the <code>--score-file</code> switch below). 
    The <code><i>&lt;distance&gt;-1</i></code> transition/transversion joint
    probability matrix for DNA is given below:

    <pre>
           A    C    G    T    
      A  .990 .002 .006 .002
      C  .002 .990 .002 .006
      G  .006 .002 .990 .002
      T  .002 .006 .002 .990
    </pre>
  </xsl:template>

  <xsl:template match="command-option[@name='paths']">
    <code>--paths single|all</code>
    - This option determines how the program computes raw scores.
    With the <code>single</code> option,
    the program computes the Viterbi score,
    which is the log-odds score associated with the single
    most likely match between the sequence and the model. The
    <code>all</code> option yields the total log-odds score,
    which is the sum of the log-odds of all sequence-to-model matches.
    The default is Viterbi scoring.
  </xsl:template>

  <xsl:template match="command-option[@name='progress']">
    <code>--progress <i>&lt;value&gt;</i></code>
    - Print to standard error a progress
    message approximately every <code><i>&lt;value&gt;</i></code> seconds.
  </xsl:template>

  <xsl:template match="command-option[@name='pseudocount']">
    <code>--pseudocounts &lt;float&gt;</code> 
    - A pseudocount to be added to each count in the motif matrix, 
    weighted by the background frequencies for the nucleotides (Dirichlet prior), 
    before converting the motif to probabilities. 
    The default value is 0.1.
  </xsl:template>

  <xsl:template match="command-option[@program='mcast' and @name='p-thresh']">
    <code>--p-thresh <i>&lt;pv&gt;</i></code>
    - Only motif occurences with p-values less than 
    <code><i>&lt;pv&gt;</i></code> will be considered in computing
    the match score. 
    The default value is 5e-4.
  </xsl:template>

  <xsl:template match="command-option[@name='pur-pyr']">
    <code>--pur-pyr &lt;float&gt;</code>
    - The ratio of the purine transition rate to pyrimidine transition rate.
    This parameter is used by the Tamura-nei model.
    The default value is 1.0.
  </xsl:template>

  <xsl:template match="command-option[@name='score-file']">
    <code>--score-file <i>&lt;score file&gt;</i></code>
    - Cause a score file (in BLAST format) to be read and used instead of
    the built-in PAM (for proteins) or transition/transversion (for DNA)
    score file. Several score files are provided (including BLOSUM62) in the
    directory <code>doc</code>.
    Other, user-provided score files may be specified as well,
    as long as they are in the proper format.
  </xsl:template>

  <xsl:template match="command-option[@name='seed']">
    <code>--seed &lt;n&gt;</code>
    - Seed for random number generator.
  </xsl:template>

  <xsl:template match="command-option[@name='sep']">
    <code>--sep</code>
    - Score reverse complement DNA strand as a separate sequence.
  </xsl:template>

  <xsl:template match="command-option[@name='seqp']">
    <code>--seqp</code>
    - Use SEQUENCE p-values for motif thresholds (default: use POSITION p-values).
  </xsl:template>

  <xsl:template match="command-option[@name='spacer-pseudo']">
    <code>--spacer-pseudo <i>&lt;value&gt;</i></code>
    - Specify the value of the
    pseudocount used in converting transition counts to spacer
    self-loop probabilities.
  </xsl:template>

  <xsl:template match="command-option[@name='synth']">
    <code>--synth</code>
    - Create synthetic sequences for estimating
    <i>E</i>-values. This is useful with small input databases where
    not enough match scores are found to estimate <i>E</i>-values.
    The <code>--bgfile</code> option must also be set when
    using this option.
  </xsl:template>

  <xsl:template match="command-option[@name='text']">
    <code>--text</code>
    Limits output to plain text sent to standard out.  For FIMO, the
    text output is unsorted, and q-values are not reported.  This mode
    allows the program to search an arbitrarily large database,
    because results are not stored in memory.
  </xsl:template>

  <xsl:template match="command-option[@name='text' and @program='mcast']">
    <code>--text</code>
    Output is plain text rather then HTML.
  </xsl:template>

  <xsl:template match="command-option[@name='trans-pseudo']">
    <code>--trans-pseudo <i>&lt;value&gt;</i></code>
    - Specify the value of the
    pseudocount used in converting transition counts to transition
    probabilities.
  </xsl:template>

  <xsl:template match="command-option[@name='transfac']">
    <code>--transfac</code>
    - The input motif file is assumed to be in 
    <a href="transfac-format.html">TRANSFAC format</a>
    and is converted to 
    <a href="meme-format.html">MEME format</a>
    before being used.
  </xsl:template>
  <xsl:template match="command-option[@name='transition-transversion']">
    <code>--transition-transversion &lt;float&gt;</code>
    - The ratio of the transition rate to the transversion rate.
    This parameter is used by the Kimura 2-parameter, F84, HKY,
    and Tamura-nei models.
    The default value is 0.5.
  </xsl:template>

  <xsl:template match="command-option[@name='type']">
    <code>--type [complete|star]</code>
    - This option specifies the topology of the model.
    The <code>complete</code> topology includes
    transitions from the end of each motif to the beginning of every
    other motif in the model (with a spacer model along each
    transition). 
    This allows for motifs that are repeated, deleted or shuffled.
    In the <code>star</code> topology the transitions from each motif lead to
    the intra-motif state. 
    The default for <code>miao</code> is the <code>complete</code> topology.
  </xsl:template>

  <xsl:template match="command-option[@name='verbosity']">
    <code>--verbosity 1|2|3|4</code>
    - Set the verbosity of status reports to standard error.
    The default level is 2.
  </xsl:template>

  <xsl:template match="command-option[@name='version']">
    <code>--version</code>
    -  Causes the program to report its version number and exit.
  </xsl:template>

  <xsl:template match="command-option[@name='zselo']">
    <code>--zselo</code>
    Spacer emission log-odds scores to be set to zero.
    This prevents regions of unusual base/residue composition 
    matching spacers well when
    the spacer emission frequencies are different than the background frequencies.
    It is particularly useful with DNA models.
  </xsl:template>

  <xsl:template match="command-option[@program='fasta-get-markov' and @name='m']">
    <code>-m &lt;n&gt;</code> 
    - Order of Markov model to use. Default value is 0.
  </xsl:template>

  <xsl:template match="command-option[@program='fasta-get-markov' and @name='norc']">
    <code>-norc</code> 
    - Do not combine forward and reverse complement frequencies.
  </xsl:template>

  <xsl:template match="command-option[@program='fasta-get-markov' and @name='p']">
    <code>-p</code> 
    - Use protein alphabet. Default is to use DNA alphabet.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='d']">
    <code>-d &lt;smin&gt; &lt;smax&gt;</code> 
    - Print length vs p-value for
    <code>&lt;smin&gt; &lt;= score &lt;= &lt;smax&gt;</code>.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='h']">
    <code>-h &lt;H&gt;</code> 
    - Initial value  of H. Default value: 1.0.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='ilh']">
    <code>-ilh &lt;ilh&gt;</code> 
    - Maximum iterations in outer L-H loop.
    Default value: 10.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='inr']">
    <code>-inr &lt;inr&gt;</code> 
    - Maximum iterations in inner N-R loop.
    Default value: 20.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='lr']">
    <code>-lr &lt;lr&gt;</code> 
    - Minimum ration of shortest to longest in  group.
    Default value: 1.50.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='p']">
    <code>-p &lt;n&gt;</code> 
    - Print p-values. Default value: off.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='r']">
    <code>-r &lt;n&gt;</code>
    - use <code>&lt;n&gt;</code> randomly selected scores.
    Default value: all.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='s']">
    <code>-s &lt;s&gt;</code>
    - Size of length groups. Default value: 10000.
  </xsl:template>

  <xsl:template match="command-option[@program='fitevd' and @name='t']">
    <code>-t</code>
    - Turn on debug tracing. Default value: off.
  </xsl:template>

  <xsl:template match="command-option[@program='gendb' and @name='bfile']"> 
    <code>--bfile <i>&lt;file&gt;</i></code>
    - File containg a background model in 
    <a href="bfile-format.html">background model format</a>.
  </xsl:template>

  <xsl:template match="command-option[@program='gendb' and @name='order']"> 
    <code>--order <i>&lt;n&gt;</i></code>
    - Use Markov model of order <code>&lt;n&gt;</code>.
    Default value: the order of the background model if one is provided,
    otherwise 0.
  </xsl:template>

  <xsl:template match="command-option[@program='gendb' and @name='type']"> 
    <code>--type [0|1|2|3|4]</code>
    - The type of sequence to generate.
    Allowed types are:
    <ul>
      <li>0 = protein (default)</li>
      <li>1 = DNA</li>
      <li>2 = codons</li>
      <li>3 = DNA without ambiguous characters</li>
      <li>4 = protein without ambigouous characters</li>
    </ul>
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='codons']">
    <code>-codons</code>
    - Print frame 0 codon usage.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='dna']">
    <code>-dna</code>
    - Print DNA frequencies in <a href="bfile-format.html">background file format</a>.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='f']">
    <code>-f</code>
    - Print letter frequences.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='ft']">
    <code>-ft</code>
    - Print letter frequences in a LaTex table.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='l']">
    <code>-l</code>
    - Just print the length of each sequence.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='nd']">
    <code>-f</code>
    - Do not print warnings about duplicate sequences.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='prot']">
    <code>-prot</code>
    - Print protein frequencies in <a href="bfile-format.html">background file format</a>.
  </xsl:template>

  <xsl:template match="command-option[@program='getsize' and @name='x']">
    <code>-x</code>
    - Translate DNA in 6 frames.
  </xsl:template>

  <xsl:template match="command-option[@name='nostatus']">
    <code>--nostatus</code>
    -  Suppresses the process information.
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='nobalance']">
    <code>--nobalance </code>
    - Disables the averaging of complementary nucleotide 
       background frequencies (AT and CG).
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='o-format']">
    <code>--o-format gff|cisml</code>
    - Output file format (default cisml).
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='bg']"> 
    <code>--bg 1|2|3</code>
    - Source used to determine background frequencies
      (1 = sequence file (default), 2 = meme motif file,
       3 = specified background file)
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='bgfile']"> 
    <code>--bgfile <i>&lt;bfile&gt;</i></code> 
    - Read background frequencies from
    <code><i>&lt;bfile&gt;</i></code>. The default is to obtain frequencies
    from the specified <i>&lt;sequence file&gt;</i>.
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='print-bg']">
    <code>--print-bg <i>&lt;destination&gt;</i> </code>
    - Prints the background frequencies to the destination.
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='scoring']">
    <code>--scoring avg-odds|max-odds </code>
    - Indicates whether the average or the maximum likelihood ratio (odds)
        score should be calculated (default avg-odds).
        If max-odds is chosen, no <i>p</i>-value will be printed.
  </xsl:template>
  
  <xsl:template match="command-option[@program='ama' and @name='pvalues']">
  <code>--pvalues</code>
    - Print the <i>p</i>-value of the average odds score in the output 
        file.  The <i>p</i>-score for a score is normally computed
        (but see <code>--gcbins</code>) assuming the
	sequences were each generated by the 0-order Markov model specified
	by the background file frequencies.
        By default, no <i>p</i>-value will be printed.
	This option is ignored if max-odds scoring is used.
  </xsl:template>

  <xsl:template match="command-option[@program='ama' and @name='rma']">
  <code>--rma</code>
    - Scale the motif affinity score by the maximum achievable score for
        each motif. This is termed the <b>R</b>elative <b>M</b>otif 
	<b>A</b>ffinity score. This allows for direct comparison between
	different motifs. By default, affinity scores are not scaled.
  </xsl:template>

  <xsl:template match="command-option[@program='ama' and @name='gcbins']">
  <code>--gcbins <i>&lt;bins&gt;</i></code>
    - Compensate <i>p</i>-values for the GC content of each
	sequence independently.  This is done by computing the
	score distributions for a range of GC values.  Using
	41 bins (recommended) computes distributions at intervals
	of 2.5% GC content.  The computation assumes that the
	ratios of G to C and A to T are both equal to 1.  This
	assumption will fail if a sequence contains far more of
	a letter than its complement.
	This option sets the <code>--pvalues</code> option.
        By default, uncompensated <i>p</i>-values are printed.
	This option is ignored if max-odds scoring is used.
  </xsl:template>

  <xsl:template match="command-option[@program='gomo' and @name='dag']">
    <code>--dag <i>&lt;go dag file&gt;</i></code>
    - Path to the optional <a href="godag-format.html" >Gene Ontology DAG</a> to be used for 
    identifying the most specific terms in the gomo xml output
    so they can be highlighted in the html output. 
  </xsl:template>
  <xsl:template match="command-option[@program='gomo' and @name='text']">
    <code>--text</code>
      - Output in tab separated values format to standard out. 
        Will not create an output directory or files.
  </xsl:template>
  <xsl:template match="command-option[@program='gomo' and @name='gs']">
  <code>--gs</code> 
    - Indicates that gene scores contained in the cisml file should be 
        used for the calculations. The default is to use the gene <i>p</i>-values.        
  </xsl:template>

  <xsl:template match="command-option[@program='gomo' and @name='shuffle_scores']">
  <code>--shuffle_scores <i>&lt;n&gt;</i></code>
  - Number of times to shuffle the sequence = score assignment and use the shuffled
      scores to generate empirical p-values.
  </xsl:template>
  
  <xsl:template match="command-option[@program='gomo' and @name='score_E_thresh']">
  <code>--score_E_thresh <i>&lt;n&gt;</i></code>
    - Threshold used on the gene score <i>E</i>-values above which all 
        <i>E</i>-values become maximal in order to reduce the impact of noise. 
        Subsequently, this results in all genes having <i>E</i>-values above 
        the threshold to obtain the same rank in the ranksum statistics.
        The threshold will be ignored when gene scores are used (--gs).         
  </xsl:template>

  <xsl:template match="command-option[@program='gomo' and @name='t']">
  <code>--t <i>&lt;n&gt;</i></code>
    - Threshold used on the q-values above which results are not considered significant and
    subsequently will not be reported. Defaults to 0.05 . To show all results use a value of 1.0 .
  </xsl:template>

  <xsl:template match="command-option[@program='gomo' and @name='min_gene_count']">
  <code>--min_gene_count <i>&lt;n&gt;</i></code>
    - Filter out GO-terms which are annotated with less genes. Defaults to 1 which shows all results.
  </xsl:template>

</xsl:stylesheet>
