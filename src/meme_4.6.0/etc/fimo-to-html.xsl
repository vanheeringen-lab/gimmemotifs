<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" 
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:cis="http://zlab.bu.edu/schema/cisml"
  xmlns:fimo="http://noble.gs.washington.edu/schema/cisml"
  xmlns:mem="http://noble.gs.washington.edu/meme"
>

  <xsl:output method="html" indent="yes" 
    doctype-public="-//W3C//DTD HTML 4.01 Transitional//EN"
    doctype-system="http://www.w3.org/TR/html4/loose.dtd"
  />

  <!-- Stylesheet processing starts here -->
  <xsl:template match="/fimo">
    <html>
      <xsl:call-template name="html_head"/>
      <body bgcolor='#D5F0FF'>  
        <a name="top_buttons"></a>
        <hr />
        <table summary="buttons" align="left" cellspacing="0">
          <tr>
            <td bgcolor="#00FFFF">
              <a href='#database_and_motifs'><b>Database and Motifs</b></a>
            </td>
            <td bgcolor="#DDFFDD">
              <a href="#sec_i"><b>High-scoring Motif Occurences</b></a>
            </td>
            <td bgcolor="#DDDDFF">
              <a href="#debugging_information"><b>Debugging Information</b></a>
            </td>
          </tr>
        </table>
        <br />
        <br />
        <!-- Create the various sub-sections of the document -->
        <xsl:call-template name="version"/>
        <xsl:call-template name="database_and_motifs"/>
        <xsl:call-template name="high_scoring_motif_occurrences"/>
        <!-- Diagrams and annotations break for long sequences 
        (H. sap. chr21 for example -->
        <!--
        <xsl:call-template name="motif_diagrams"/>
        <xsl:call-template name="annotations"/>
        -->
        <xsl:call-template name="debugging_information"/>
        <!-- Buttons are pointless since there is only one section
        <xsl:call-template name="button_help"/>
        -->
      </body>
    </html>
  </xsl:template>

  <xsl:template name="html_head">
    <!-- This template prints the HTML head block, including the document level CSS. -->
    <head>
      <meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"/>
      <title>FIMO Results</title>
      <style type="text/css">
        td.left {text-align: left;}
        td.right {text-align: right; padding-right: 1cm;}
      </style>
    </head>
  </xsl:template>

  <xsl:template name="version">
    <!-- 
      This template prints the HTML describing the version of FIMO 
      that generated this document. 
    -->
    <hr />
    <center>
      <big><b>FIMO - Motif search tool</b></big>
    </center>
    <hr />
    <p>
      FIMO version <xsl:value-of select="@version"/>, 
      (Release date: <xsl:value-of select="@release"/>)
    </p>
    <p>
      For further information on how to interpret these results
      or to get a copy of the FIMO software please access
      <a href="http://meme.nbcr.net">http://meme.nbcr.net</a>
    </p>
  </xsl:template>

  <xsl:template name="reference">
    <!-- This template prints the instructions for citing FIMO. -->
    <hr/>
    <center>
      <a name="reference"/>
      <big><b>REFERENCE</b></big>
    </center>
    <hr/>
    <p>
    If you use this program in your research, please cite:
    </p>
    <p>
    </p>
  </xsl:template>

  <xsl:template name="database_and_motifs">
    <!-- 
      This template prints the HTML descibing the input sequences
      and motifs that were used to generate this document.
     -->
    <hr/>
    <center>
      <big>
        <b><a name="database_and_motifs">DATABASE AND MOTIFS</a></b>
      </big>
    </center>
    <hr/>
    <div style="padding-left: 0.75in; line-height: 1em; font-family: monospace;">
      <p>
        DATABASE 
        <xsl:value-of select="/fimo/settings/setting[@name='sequence file name']"/>
        <br />
        Database contains 
        <xsl:value-of select="/fimo/sequence-data/@num-sequences"/>
        sequences,
        <xsl:value-of select="/fimo/sequence-data/@num-residues"/>
        residues
      </p>
      <p>
        MOTIFS 
        <xsl:value-of select="/fimo/settings/setting[@name='MEME file name']"/>
        <xsl:text> </xsl:text>
        (<xsl:value-of select="/fimo/alphabet"/>)
        <table>
          <thead>
          <tr>
            <th style="border-bottom: 1px dashed;">MOTIF</th>
            <th style="border-bottom: 1px dashed; padding-left: 1em;">WIDTH</th>
            <th style="border-bottom: 1px dashed; padding-left: 1em;" align="left">
              BEST POSSIBLE MATCH
            </th>
          </tr>
          </thead>
          <tbody>
          <xsl:for-each select="/fimo/motif">
            <tr>
              <td align="right"><xsl:value-of select="./@name"/></td>
              <td align="right" style="padding-left: 1em;"><xsl:value-of select="./@width"/></td>
              <td align="left" style="padding-left: 1em;">
                <xsl:value-of select="./@best-possible-match"/>
              </td>
            </tr>
          </xsl:for-each>
          </tbody>
        </table>
      </p>
      <p>
        Random model letter frequencies 
        (from <xsl:value-of select="/fimo/background/@source"/>):
        <br />
        <xsl:call-template name="bg_freq"/>
      </p>
    </div>
  </xsl:template>

  <xsl:template name="bg_freq">
    <!-- This template prints the HTML listing the background frequencies for each symbol in the alphabet. -->
    <xsl:for-each select="/fimo/background/value">
      <xsl:value-of select="./@letter"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="format-number(., '0.000')"/>
      <xsl:text> </xsl:text>
      <!-- For amino acids insert a line break very 9 characters. -->
      <xsl:variable name="letter_index" select="position()" />
      <xsl:if test="$letter_index mod 9 = 0">
        <xsl:text>&#10;</xsl:text>
        <br/>
        <xsl:text>&#10;</xsl:text>
      </xsl:if>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="high_scoring_motif_occurrences">
    <xsl:variable name="alphabet">
      <xsl:value-of select="/fimo/alphabet"/>
    </xsl:variable>
    <xsl:variable name="num_qvalues">
      <xsl:value-of
      select="count(document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence/cis:matched-element/mem:qvalue)"/>
    </xsl:variable>
    <!-- This template prints the matches to the pattern. -->
    <hr/>
    <center>
      <big><b><a name="sec_i">SECTION I: HIGH-SCORING MOTIF OCCURENCES</a></b></big>
    </center>
    <hr/>
    <ul>
    <li> Each of the following 
    <xsl:value-of select="count(document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence/cis:matched-element)"/>
    motif occurrences has 
    <xsl:choose>
      <xsl:when test="/fimo/settings/setting[@name='output q-threshold set'] = 'true'">
        q-value less than
        <xsl:value-of select="/fimo/settings/setting[@name='output q-value threshold' = 'true']"/>
      </xsl:when>
      <xsl:otherwise>
        p-value less than 
        <xsl:value-of select="/fimo/settings/setting[@name='output p-value threshold']"/>
      </xsl:otherwise>
    </xsl:choose>
    </li>
    <li> The p-value of a motif occurrence is defined as the
    probability of a random sequence of the same length as the motif
    matching that position of the sequence with as good or better a score.
    </li>
    <li> The score for the match of a position in a sequence to a motif
    is computed by summing the appropriate entries from each column of
    the position-dependent scoring matrix that represents the motif.
    </li>
    <xsl:if test="$num_qvalues > 0">
    <li> The q-value of a motif occurrence is defined as the
    false discovery rate if the occurrence is accepted as significant.
    </li>
    </xsl:if>
    <li>The table is sorted by increasing p-value.</li>
    <li>If the start position is larger than the end position,
    the motif occurrence is on the reverse strand.
    </li>
    </ul>
    <table border="1">
      <thead>
        <tr>
          <th>Motif</th>
          <th>Sequence Name</th>
          <xsl:if test="$alphabet='nucleotide'">
            <th>Strand</th>
          </xsl:if>
          <th>Start</th>
          <th>End</th>
          <th>p-value</th>
          <xsl:if test="$num_qvalues > 0">
            <th>q-value</th>
          </xsl:if>
          <th>Matched Sequence</th>
        </tr>
      </thead>
      <tbody>
        <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence/cis:matched-element">
          <xsl:sort select="@pvalue" data-type="number" order="ascending"/>
          <xsl:variable name="seq_name" select="../@name"/>
          <tr>
            <td align="left">
             <a name="occurrences_{$seq_name}"><xsl:value-of select="../../@name"/></a>
            </td>
            <td align="left">
             <xsl:value-of select="../@name"/>
            </td>
            <!-- start  and stop positions -->
            <!-- CISML uses start > stop to indicate - strand -->
            <!-- For HTML we want start < stop, indicate strand in strand column -->
            <xsl:if test="$alphabet='nucleotide'">
              <xsl:choose>
                <xsl:when test="@start &gt; @stop">
                  <!-- start > stop so - strand -->
                 <td align="center"><xsl:text>&#8722;</xsl:text></td>
                 </xsl:when>
                 <xsl:otherwise>
                   <td align="center"><xsl:text>+</xsl:text></td>
                 </xsl:otherwise>
              </xsl:choose>
            </xsl:if>
             <xsl:choose>
               <xsl:when test="@start &gt; @stop">
                 <!-- start > stop so - strand -->
                 <td align="left">
                   <xsl:value-of select="@stop"/>
                 </td>
                 <td align="left">
                   <xsl:value-of select="@start"/>
                 </td>
               </xsl:when>
               <xsl:otherwise>
                 <td align="left">
                   <xsl:value-of select="@start"/>
                 </td>
                 <td align="left">
                   <xsl:value-of select="@stop"/>
                 </td>
               </xsl:otherwise>
            </xsl:choose>
            <td align="left">
             <xsl:value-of select="@pvalue"/>
            </td>
            <xsl:if test="$num_qvalues > 0">
              <td align="left">
               <xsl:value-of select="./mem:qvalue"/>
              </td>
            </xsl:if>
            <td align="left">
             <code style="font-size: x-large;"><xsl:value-of select="./cis:sequence"/></code>
            </td>
          </tr>
        </xsl:for-each>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="motif_diagrams">
    <!--
      This template draws the bead-on-string diagram for each
      of the scanned sequence.
    -->
    <hr/>
    <center>
      <big><b><a name="sec_ii">SECTION II: MOTIF DIAGRAMS</a></b></big>
    </center>
    <hr/>
    <table border="1">
      <thead>
      <tr>
        <th align="left"><br/>Name</th>
        <th align="left"><br/>Motifs</th>
        <th></th>
      </tr>
      </thead>
      <tbody>
        <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[not (@name = preceding::cis:scanned-sequence/@name)]">
          <xsl:sort select="@name" order="ascending"/>
          <xsl:variable name="sequence_width" select="@length"/>
          <xsl:variable name="sequence_name" select="@name"/>
          <tr>
            <td><a name="diagram_{$sequence_name}"><xsl:value-of select="$sequence_name"/></a></td>
            <td>
            <table>
            <tr>
              <xsl:call-template name="draw_diagram_for_sequence">
                <xsl:with-param name="sequence_name" select="$sequence_name"/>
                <xsl:with-param name="sequence_width" select="@length"/>
              </xsl:call-template>
            </tr>
            </table>
            </td>
          </tr>
        </xsl:for-each>
        <tr>
          <td rowspan="2" style="color: blue;">SCALE</td>
          <td>
            <table>
              <tr>
                <td>
                  <xsl:call-template name="ticks">
                    <xsl:with-param name="count" select="1"/>
                    <xsl:with-param name="limit">
                      <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence/@length">
                        <xsl:sort data-type="number" order="descending"/>
                        <xsl:if test="position()=1"><xsl:value-of select="."/></xsl:if>
                      </xsl:for-each>
                    </xsl:with-param>
                  </xsl:call-template>
                </td>
              </tr>
              <tr>
                <td>
                  <xsl:call-template name="scale">
                    <xsl:with-param name="count" select="1"/>
                    <xsl:with-param name="limit">
                      <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence/@length">
                        <xsl:sort data-type="number" order="descending"/>
                        <xsl:if test="position()=1"><xsl:value-of select="."/></xsl:if>
                      </xsl:for-each>
                    </xsl:with-param>
                  </xsl:call-template>
                </td>
              </tr>
            </table>
          </td>
        </tr>
      </tbody>
    </table>
  </xsl:template>

  <xsl:template name="draw_diagram_for_sequence">
    <!--
      This template draws the bead-on-string diagram for
      a single sequence.
    -->
    <xsl:param name="sequence_name"/>
    <xsl:param name="sequence_width"/>
    <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name = $sequence_name]/cis:matched-element">
      <xsl:sort select="number(@start)" data-type="number" order="ascending"/>
      <xsl:variable name="starting">
        <xsl:choose>
          <xsl:when test="number(@start) &lt; number(@stop)">
            <xsl:value-of select="@start"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="@stop"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="ending">
        <xsl:choose>
          <xsl:when test="number(@start) &lt; number(@stop)">
            <xsl:value-of select="@stop"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="@start"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="current_element" select="position()"/>
      <xsl:variable name="previous_end">
        <xsl:choose>
          <xsl:when test="position() = 1">0</xsl:when>
          <xsl:otherwise>
            <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name = $sequence_name]/cis:matched-element">
              <xsl:sort select="number(@start)" data-type="number" order="ascending"/>
              <xsl:if test="position() = ($current_element - 1)">
                <xsl:choose>
                  <xsl:when test="number(@start) &lt; number(@stop)">
                    <xsl:value-of select="@stop"/>
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="@start"/>
                  </xsl:otherwise>
                </xsl:choose>
              </xsl:if>
            </xsl:for-each>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="motif_name">
        <xsl:value-of select="../../@name"/>
      </xsl:variable>
      <xsl:variable name="motif_index">
        <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern">
          <xsl:if test="@name = $motif_name"><xsl:value-of select="position()"/></xsl:if>
        </xsl:for-each>
        <xsl:value-of select="@name"/>
      </xsl:variable>
      <xsl:variable name="motif_color">
        <xsl:call-template name="pick_color">
          <xsl:with-param name="index" select="$motif_index"/>
        </xsl:call-template>
      </xsl:variable>
      <xsl:variable name="font_color">
        <xsl:call-template name="pick_font_color">
          <xsl:with-param name="index" select="$motif_index"/>
        </xsl:call-template>
      </xsl:variable>
      <xsl:variable name="strand">
        <xsl:choose>
          <xsl:when test="number(@start) &lt; number(@stop)">
            <xsl:text>+</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>-</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="alphabet" select="document('fimo.xml')/fimo/alphabet"/>
      <td style="width: {2 * ($starting - $previous_end) }px;"><hr size="4" noshade="noshade"/></td>
       <td style="width: {2 * ($ending - $starting)}px; background: {$motif_color}; color: {$font_color}; text-align: center;">
         <xsl:if test="$alphabet = 'nucleotide'"><xsl:value-of select="$strand"/></xsl:if>
         <xsl:value-of select="../../@name"/>
      </td>
      <xsl:if test="position() = last()">
        <td style="width: {2 * ($sequence_width - $ending) }px;"><hr size="4" noshade="noshade"/></td>
      </xsl:if>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="ticks">
      <!--
        This template is used to print the tick marks in the 
        scale at the bottom of the bead-on-string daigrams.
      -->
      <xsl:param name="count" select="1"/>
      <xsl:param name="limit" select="1"/>
      <xsl:if test="$count &lt; $limit">
        <td style="width: 50px;" align="left"><font color="blue">|</font></td>
        <xsl:call-template name="ticks">
          <xsl:with-param name="count" select="$count + 24"/>
          <xsl:with-param name="limit" select="$limit"/>
        </xsl:call-template>
      </xsl:if>
  </xsl:template>

  <xsl:template name="scale">
      <!--
        This template is used to print the text in the 
        scale at the bottom of the bead-on-string daigrams.
      -->
      <xsl:param name="count" select="1"/>
      <xsl:param name="limit" select="1"/>
      <xsl:if test="$count &lt; $limit">
        <td style="width: 50px;" align="left"><font color="blue"><xsl:value-of select="$count"/></font></td>
        <xsl:if test="$count = 1">
          <xsl:call-template name="scale">
            <xsl:with-param name="limit" select="$limit"/>
            <xsl:with-param name="count" select="$count + 24"/>
          </xsl:call-template>
        </xsl:if>
        <xsl:if test="$count > 1">
          <xsl:call-template name="scale">
            <xsl:with-param name="limit" select="$limit"/>
            <xsl:with-param name="count" select="$count + 25"/>
          </xsl:call-template>
        </xsl:if>
      </xsl:if>
  </xsl:template>

  <xsl:template name="pick_color">
    <!-- 
      This template is used to set the color of motif markers
      in the bead-on-string diagrams.
    -->
    <xsl:param name="index" select="1"/>
    <xsl:choose>
      <xsl:when test="($index mod 16) = 0">white</xsl:when>
      <xsl:when test="($index mod 16) = 1">aqua</xsl:when>
      <xsl:when test="($index mod 16) = 2">blue</xsl:when>
      <xsl:when test="($index mod 16) = 3">red</xsl:when>
      <xsl:when test="($index mod 16) = 4">fuchsia</xsl:when>
      <xsl:when test="($index mod 16) = 5">yellow</xsl:when>
      <xsl:when test="($index mod 16) = 6">lime</xsl:when>
      <xsl:when test="($index mod 16) = 7">teal</xsl:when>
      <xsl:when test="($index mod 16) = 8">#444444</xsl:when>
      <xsl:when test="($index mod 16) = 9">green</xsl:when>
      <xsl:when test="($index mod 16) = 10">silver</xsl:when>
      <xsl:when test="($index mod 16) = 11">purple</xsl:when>
      <xsl:when test="($index mod 16) = 12">olive</xsl:when>
      <xsl:when test="($index mod 16) = 13">navy</xsl:when>
      <xsl:when test="($index mod 16) = 14">maroon</xsl:when>
      <xsl:when test="($index mod 16) = 15">black</xsl:when>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="pick_font_color">
    <!-- 
      This template is used to set the font color in motif markers
      in the bead-on-string diagrams.
    -->
    <xsl:param name="index" select="1"/>
    <xsl:choose>
      <xsl:when test="($index mod 16) = 0">black</xsl:when>
      <xsl:when test="($index mod 16) = 1">black</xsl:when>
      <xsl:when test="($index mod 16) = 2">white</xsl:when>
      <xsl:when test="($index mod 16) = 3">black</xsl:when>
      <xsl:when test="($index mod 16) = 4">black</xsl:when>
      <xsl:when test="($index mod 16) = 5">black</xsl:when>
      <xsl:when test="($index mod 16) = 6">black</xsl:when>
      <xsl:when test="($index mod 16) = 7">white</xsl:when>
      <xsl:when test="($index mod 16) = 8">white</xsl:when>
      <xsl:when test="($index mod 16) = 9">black</xsl:when>
      <xsl:when test="($index mod 16) = 10">black</xsl:when>
      <xsl:when test="($index mod 16) = 11">white</xsl:when>
      <xsl:when test="($index mod 16) = 12">white</xsl:when>
      <xsl:when test="($index mod 16) = 13">white</xsl:when>
      <xsl:when test="($index mod 16) = 14">white</xsl:when>
      <xsl:when test="($index mod 16) = 15">white</xsl:when>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="annotations">
    <!--
      This template draws the motif/sequence annotations for each
      of the scanned sequence less than 10kb.
    -->
    <hr/>
    <center>
      <big><b><a name="sec_iii">SECTION III: SEQUENCE/MOTIF ANNOTATIONS</a></b></big>
    </center>
    <hr/>
    <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[not (@name = preceding::cis:scanned-sequence/@name)]">
      <xsl:sort select="@name" order="ascending"/>
      <xsl:variable name="sequence_length" select="@length"/>
      <xsl:variable name="sequence_name" select="@name"/>
      <p>
        <a name="annotation_{$sequence_name}"><xsl:value-of select="$sequence_name"/></a>
        <br />
        <span style="background-color: ##DDDD88;"><a href="#occurrences_{$sequence_name}">S</a></span>
        <span style="background-color: #DDFFDD;"><a href="#diagram_{$sequence_name}">D</a></span>
        <span style="background-color: #FFFFFF;"><a href="#button_help">?</a></span>
      </p>
      <p>
      Length = <xsl:value-of select="$sequence_length"/>
      </p>
      <p>
      String Length = <xsl:value-of
      select="string-length(normalize-space(document('matched_sequences.xml')/matched-sequences/sequence[@name=$sequence_name]))"/>
      </p>
      <pre>
      <table style="font-family: monospace;">
      <tbody>
      <xsl:call-template name="wrap_sequence">
        <xsl:with-param name="sequence" 
           select="normalize-space(document('matched_sequences.xml')/matched-sequences/sequence[@name=$sequence_name])"/>
        <xsl:with-param name="sequence_index" select="1"/>
        <xsl:with-param name="sequence_name" select="$sequence_name"/>
      </xsl:call-template>
      </tbody>
      </table>
      </pre>
      <hr/>
    </xsl:for-each>
       
  </xsl:template>

  <xsl:template name="debugging_information">
    <!-- This template print the HTML describing the FIMO input parameters. -->
    <hr/>
    <center>
      <big><b><a name="debugging_information">DEBUGGING INFORMATION</a></b></big>
    </center>
    <hr/>
    <p>
    Command line:
    </p>
    <pre>
    <xsl:value-of select="/fimo/command-line"/>
    </pre>
    <p>
    Settings:
    </p>
    <pre>
    <table>
    <xsl:apply-templates select="/fimo/settings/setting[position() mod 3 = 1]"/>
    </table>
    </pre>
    <p>
      This information can be useful in the event you wish to report a
      problem with the FIMO software.
    </p>
    <hr />
    <span style="background-color: #DDDDFF"><a href="#top_buttons"><b>Go to top</b></a></span>

  </xsl:template>

  <xsl:template match="setting">
    <!-- This template prints the program settings in 2 columns -->
    <tr>
      <td style="padding-right: 2em"><xsl:value-of select="@name"/> = <xsl:value-of select="."/></td>
      <xsl:choose>
      <xsl:when test="count(following-sibling::*[name()=name(current())])">
        <td style="padding-left: 5em; padding-right: 2em">
          <xsl:value-of select="following-sibling::setting[1]/@name"/> = <xsl:value-of select="following-sibling::setting[1]"/>
        </td>
        <td style="padding-left: 5em; padding-right: 2em">
          <xsl:value-of select="following-sibling::setting[2]/@name"/> = <xsl:value-of select="following-sibling::setting[2]"/>
        </td>
      </xsl:when>
      <xsl:otherwise>
        <td style="padding-left: 5em; padding-right: 2em"></td>
        <td align="right"></td>
      </xsl:otherwise>
      </xsl:choose>
    </tr>
  </xsl:template>

  <xsl:template name="button_help">
    <hr />
    <center>
    <h3><a name="button_help">Button Help</a></h3>
    </center>
    <hr />
    <span style="background-color: #DDDDFF">S</span>
    Links to Entrez database at <a href="http://www.ncbi.nlm.nih.gov">NCBI</a>
    <span style="background-color: #DDDD88">S</span>
    Links to motif occurences table(<a href="#sec_i">SECTION I</a>)
    <br/>
    <span style="background-color: #DDFFDD">D</span>
    Links to motif occurences diagram (<a href="#sec_ii">SECTION II</a>)
    <br />
    <span style="background-color: #FFDDDD">A</span>
    Links to sequence/motif annotated alignments (<a href="#sec_iii">SECTION III</a>) 
    <br />
    <span style="background-color: #FFFFFF">?</span>
    This information 
    <br />

    <hr />
    <span style="background-color: #DDDDFF"><a href="#top_buttons"><b>Go to top</b></a></span>
  </xsl:template>

  <xsl:template name="wrap_sequence">
    <xsl:param name="sequence"/>
    <xsl:param name="sequence_index"/>
    <xsl:param name="sequence_name"/>
    <xsl:variable name="current_sequence">
      <xsl:value-of select="substring($sequence, 1, 75)"/>
    </xsl:variable>
    <tr>
      <td></td>
      <td>
        <xsl:call-template name="annotate_motif_names">
          <xsl:with-param name="sequence_index" select="$sequence_index"/>
          <xsl:with-param name="sequence_name" select="$sequence_name"/>
          <xsl:with-param name="current_sequence" select="$current_sequence" />
        </xsl:call-template>
      </td>
    </tr>
    <tr>
    <td style="width: 6em;"><xsl:value-of select ="$sequence_index"/></td>
    <td><xsl:value-of select="$current_sequence"/></td>
    </tr>
    <xsl:variable name="new_sequence" select="substring($sequence, 76)"/>
    <xsl:if test="string-length($new_sequence)">
      <xsl:call-template name="wrap_sequence">
        <xsl:with-param name="sequence" select="$new_sequence"/>
        <xsl:with-param name="sequence_index" select="number($sequence_index) + 75"/>
        <xsl:with-param name="sequence_name" select="$sequence_name"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="annotate_motif_names">
    <xsl:param name="sequence_index"/>
    <xsl:param name="sequence_name"/>
    <xsl:param name="input_sequence"/>
    <xsl:param name="output_sequence"/>
    <xsl:choose>
      <xsl:when test="count(document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name=$sequence_name]) = 1">
      </xsl:when>
    </xsl:choose>
      <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name=$sequence_name]">
        <xsl:for-each select="./cis:matched-element[(@start &gt;= $sequence_index) and (@start &lt; ($sequence_index + 75))]">
          <xsl:sort select="@start" data-type="number" />
          <xsl:text>     </xsl:text><xsl:value-of select="../../@name"/>
        </xsl:for-each>
      </xsl:for-each>
  </xsl:template>

  <xsl:template name="annotate_pvalues">
    <xsl:param name="sequence_index"/>
    <xsl:param name="sequence_name"/>
    <pre>
      <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name=$sequence_name]">
        <xsl:for-each select="./cis:matched-element[(@start &gt;= $sequence_index) and (@start &lt; ($sequence_index + 75))]">
          <xsl:sort select="@start" data-type="number" />
          <xsl:text>    </xsl:text><xsl:value-of select="@pvalue"/>
        </xsl:for-each>
      </xsl:for-each>
    </pre>
  </xsl:template>

  <xsl:template name="annotate_best_match">
    <xsl:param name="sequence_index"/>
    <xsl:param name="sequence_name"/>
    <pre>
      <xsl:for-each select="document('cisml.xml')/cis:cis-element-search/cis:pattern/cis:scanned-sequence[@name=$sequence_name]">
        <xsl:for-each select="./cis:matched-element[(@start &gt;= $sequence_index) and (@start &lt; ($sequence_index + 75))]">
          <xsl:sort select="@start" data-type="number" />
          <xsl:text>     </xsl:text><xsl:value-of select="document('fimo.xml')/fimo/motif[@name=./@name]/@best-possible-match"/>
        </xsl:for-each>
      </xsl:for-each>
    </pre>
  </xsl:template>
</xsl:stylesheet>

