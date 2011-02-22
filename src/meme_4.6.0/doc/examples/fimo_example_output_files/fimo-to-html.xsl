<?xml version="1.0"?>

<xsl:stylesheet version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
 xmlns:mem="http://noble.gs.washington.edu/meme"
>

  <xsl:output method="html" />

  <!-- A stylesheet to generate an HTML report of p-values and q-value for motifs
  scanned over a group of sequences using fimo -->

  <xsl:template match="cis:cis-element-search">
    <html>
      <head>
        <title>Sequence Analysis with <xsl:value-of select="cis:program-name"/></title>
        <link href="cisml.css" rel="stylesheet" type="text/css"></link>
      </head> 
      <body>
        <h2 align="center">Sequence Analysis with <xsl:value-of select="cis:program-name"/></h2>
        <p>The score for the match of a position in a sequence to a motif
        is computed by by summing the appropriate entries from each column of
        the position-dependent scoring matrix that represents the motif.
        </p>
        <p>The p-value of a motif occurrence is the
        probability of a random sequence of the same length as the motif
        matching that position of the sequence with a score at least as good.
        </p>
        <p> The q-value of a motif occurrence is the estimated
        false discovery rate if the occurrence is accepted as significant.
        See Storey JD, Tibshirani R. Statistical significance for genome-wide studies. 
        <i>Proc. Natl Acad. Sci. USA (2003) 100:9440–9445</i>
        </p>
        <p>The table is sorted by increasing p-value.</p>
        <p>If the start position is larger then the end position,
        the motif occurrence is on the reverse strand.
        </p>
        <table>
          <tr><th>Pattern Name</th><th>Sequence Name</th><th>Start</th><th>Stop</th>
          <th>Score</th><th>p-value</th><th>q-value</th><th>Matched Sequence</th></tr>
          <xsl:apply-templates select="cis:pattern/cis:scanned-sequence/cis:matched-element">
            <xsl:sort order="ascending" data-type="number" select="../../@name"/>
            <xsl:sort order="ascending" data-type="number" select="@pvalue"/>
          </xsl:apply-templates>
        </table>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="cis:matched-element">
    <tr>
      <td><xsl:value-of select="../../@name" /></td>
      <td><xsl:value-of select="../@name" /></td>
      <td><xsl:value-of select="@start" /></td>
      <td><xsl:value-of select="@stop" /></td>
      <td><xsl:value-of select="@score" /></td>
      <td><xsl:value-of select="@pvalue" /></td>
      <td><xsl:value-of select="mem:qvalue" /></td>
      <td><xsl:value-of select="cis:sequence" /></td>
    </tr>
  </xsl:template>

</xsl:stylesheet>
