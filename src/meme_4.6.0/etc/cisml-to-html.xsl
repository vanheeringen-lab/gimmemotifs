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
