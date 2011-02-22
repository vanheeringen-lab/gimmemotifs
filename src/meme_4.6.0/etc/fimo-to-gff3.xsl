<?xml version="1.0"?>
<!-- A stylesheet to create a GFF3 formated file from FIMO CisML.-->
<xsl:stylesheet
 version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
 xmlns:mem="http://noble.gs.washington.edu/meme"
>
  <xsl:include href="cisml-to-gff3.xsl"/>
  <xsl:output method="text"/>

  <xsl:template match="/fimo">
    <xsl:apply-templates select="cisml-file"/>
  </xsl:template>

  <xsl:template match="cisml-file">
    <xsl:apply-templates select="document(text())">
      <xsl:with-param name="alphabet" select="../alphabet/text()"/>
    </xsl:apply-templates>
  </xsl:template>

</xsl:stylesheet>
