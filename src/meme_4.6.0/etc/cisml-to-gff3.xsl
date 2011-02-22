<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
 xmlns:mem="http://noble.gs.washington.edu/meme"
>
  <xsl:output method="text" />

  <!-- A stylesheet to create a GFF formated file from CisML.-->
  <xsl:param name="alphabet" select="nucleotide"/>

  <xsl:template match="/*">
    <xsl:text>##gff-version 3&#10;</xsl:text>
    <xsl:apply-templates select="cis:pattern/cis:scanned-sequence/cis:matched-element">
      <xsl:with-param name="alphabet" select="$alphabet"/>
      <xsl:sort order="ascending" data-type="text" select="../@name"/>
      <xsl:sort order="ascending" data-type="text" select="../../@name"/>
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
   <xsl:choose>
     <xsl:when test="$alphabet='nucleotide'">
        <xsl:text>nucleotide_motif</xsl:text>
     </xsl:when>
     <xsl:otherwise>
        <xsl:text>polypeptide_motif</xsl:text>
     </xsl:otherwise>
   </xsl:choose>
   <xsl:text>&#9;</xsl:text>
   <!-- start  and stop positions -->
   <!-- CISML uses start > stop to indicate - strand -->
   <!-- GFF always wants start < stop, indicate strand in strand column -->
   <xsl:choose>
     <xsl:when test="@start &lt; @stop">
       <!-- start < stop so + strand -->
       <xsl:value-of select="@start"/>
       <xsl:text>&#9;</xsl:text>
       <!-- end -->
       <xsl:value-of select="@stop"/>
       <xsl:text>&#9;</xsl:text>
       <!-- score -->
       <xsl:value-of select="@score"/>
       <!-- strand -->
       <xsl:text>&#9;+</xsl:text>
      </xsl:when>
      <xsl:otherwise>
       <!-- start > stop so - strand -->
       <!-- end -->
       <xsl:value-of select="@stop"/>
       <xsl:text>&#9;</xsl:text>
       <!-- start -->
       <xsl:value-of select="@start"/>
       <xsl:text>&#9;</xsl:text>
       <!-- score -->
       <xsl:value-of select="@score"/>
       <!-- strand -->
       <xsl:text>&#9;-</xsl:text>
      </xsl:otherwise>
   </xsl:choose>
   <!-- frame -->
   <xsl:text>&#9;.</xsl:text>
   <xsl:text>&#9;</xsl:text>
   <!-- attributes -->
   <!-- attribute motif_name -->
   <xsl:text>motif_name=</xsl:text>
   <xsl:value-of select="../../@name"/>
   <xsl:text>;</xsl:text>
   <!-- attribute pvalue -->
   <xsl:text>pvalue=</xsl:text>
   <xsl:value-of select="@pvalue"/>
   <xsl:text>;</xsl:text>
   <!-- attribute qvalue -->
   <xsl:text>qvalue=</xsl:text>
   <xsl:value-of select="./mem:qvalue"/>
   <xsl:text>;</xsl:text>
   <!-- attribute sequence -->
   <xsl:text>sequence=</xsl:text>
   <xsl:value-of select="./cis:sequence"/>
   <xsl:text>&#10;</xsl:text>
  </xsl:template>

</xsl:stylesheet>
