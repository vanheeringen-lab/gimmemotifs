<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
>
  <xsl:output method="text" />

  <!-- A stylesheet to create a WIG formated file from CisML file.-->

  <xsl:template match="/*">
    <xsl:apply-templates select="cis:pattern">
      <xsl:sort order="ascending" data-type="text" select="@name"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="cis:pattern">
    <!-- Each motif goes into a separate track -->
    <xsl:text>track type=wiggle_0 name=&#34;motif </xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>&#34; description=&#34;</xsl:text>
    <xsl:value-of select="../cis:program-name"/>
    <xsl:text> scan of motif </xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>&#34;</xsl:text>
    <xsl:text>&#10;</xsl:text>
    <xsl:apply-templates select="*/cis:matched-element">
      <xsl:sort order="descending" data-type="number" select="@score"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="cis:matched-element">
   <!-- 
     Each WIG format item is split over two lines. 
     The first line contains the chromosome name,
     and the width of the item. The second line contain
     the starting position and the score.
   -->
   <xsl:text>variablestep </xsl:text>
   <xsl:text>chrom=</xsl:text>
   <xsl:value-of select="../@name"/>
   <xsl:text> span=</xsl:text>
   <xsl:choose>
     <!-- 
       In CisML having @stop < @start indicates the motif occurence
       is on the '-' strand. For WIG format the strand is ignored,
       so the starting position is always the smaller of @start and @stop.                            
     -->
     <xsl:when test="@start &lt; @stop">
       <!-- Calculate the width of the item -->
       <xsl:value-of select="@stop - @start"/>
       <xsl:text>&#10;</xsl:text>
       <xsl:value-of select="@start"/>
       <xsl:text>&#9;</xsl:text>
     </xsl:when>
     <xsl:otherwise>
       <!-- Calculate the width of the item -->
       <xsl:value-of select="@start - @stop"/>
       <xsl:text>&#10;</xsl:text>
       <xsl:value-of select="@stop"/>
       <xsl:text>&#9;</xsl:text>
     </xsl:otherwise>
   </xsl:choose>
   <xsl:value-of select="@score"/>
   <xsl:text>&#10;</xsl:text>
  </xsl:template>

</xsl:stylesheet>
