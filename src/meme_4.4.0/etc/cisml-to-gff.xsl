<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
 xmlns:mem="http://noble.gs.washington.edu/meme"
>
  <xsl:output method="text" />

  <!-- A stylesheet to create a GFF formated file from CisML.-->

  <xsl:template match="/*">
    <xsl:apply-templates select="cis:pattern/cis:scanned-sequence/cis:matched-element">
      <xsl:sort order="ascending" data-type="number" select="../../@name"/>
      <xsl:sort order="ascending" data-type="number" select="@pvalue"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="cis:matched-element">
   <!-- sequence name -->
   <xsl:value-of select="../@name"/>
   <xsl:text>&#9;</xsl:text>
   <!-- source -->
   <xsl:value-of select="/cis:cis-element-search/cis:program-name" />
   <xsl:text>&#9;</xsl:text>
   <!-- feature-->
   <xsl:text>motif </xsl:text><xsl:value-of select="../../@name"/><xsl:text>&#9;</xsl:text>
   <!-- start -->
   <xsl:value-of select="@start"/>
   <xsl:text>&#9;</xsl:text>
   <!-- end -->
   <xsl:value-of select="@stop"/>
   <xsl:text>&#9;</xsl:text>
   <!-- score -->
   <xsl:value-of select="@score"/>
   <!-- strand -->
   <xsl:text>&#9;.</xsl:text>
   <!-- frame -->
   <xsl:text>&#9;.</xsl:text>
   <xsl:text>&#9;</xsl:text>
   <!-- attribute -->
   <!-- attribute pvalue -->
   <xsl:text>pvalue </xsl:text>
   <xsl:value-of select="@pvalue"/>
   <xsl:text>; </xsl:text>
   <!-- attribute qvalue -->
   <xsl:text>qvalue </xsl:text>
   <xsl:value-of select="./mem:qvalue"/>
   <xsl:text>; </xsl:text>
   <!-- attribute sequence -->
   <xsl:text>sequence </xsl:text>
   <xsl:value-of select="./cis:sequence"/>
   <xsl:text>;</xsl:text>
   <xsl:text>&#10;</xsl:text>
  </xsl:template>

</xsl:stylesheet>
