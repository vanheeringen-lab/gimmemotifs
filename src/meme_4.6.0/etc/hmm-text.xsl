<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text"/>

<!-- This template matches the top level of a MEME-HMM document -->
<xsl:template match="/MEME_HMM">
<xsl:value-of select="/MEME_HMM/@version"/>
<xsl:text>
</xsl:text>
<xsl:text>type: </xsl:text>
<xsl:value-of select="/MEME_HMM/@type"/>
<xsl:text>
</xsl:text>
<xsl:text>num_states: </xsl:text>
<xsl:value-of select="/MEME_HMM/states/@num_states"/>
<xsl:text>
</xsl:text>
<xsl:text>num_spacers: </xsl:text>
<xsl:value-of select="/MEME_HMM/states/@num_spacers"/>
<xsl:text>
</xsl:text>
<xsl:text>spacer_length: </xsl:text><xsl:value-of select="/MEME_HMM/states/@spacer_length"/>
<xsl:text>
</xsl:text>
<xsl:text>alph_size: </xsl:text><xsl:value-of select="/MEME_HMM/alphabet/@length"/>
<xsl:text>
</xsl:text>
<xsl:text>alphabet: </xsl:text><xsl:for-each select="/MEME_HMM/alphabet/letter"><xsl:value-of select="@id"/></xsl:for-each>
<xsl:text>
</xsl:text>
<xsl:text>background: </xsl:text><xsl:call-template name="background_frequencies"/>
<xsl:text>
</xsl:text>
<xsl:text>motif_file: </xsl:text><xsl:value-of select="/MEME_HMM/@motif_file"/>
<xsl:for-each select="/MEME_HMM/states/state">
  <xsl:call-template name="state"><xsl:with-param name="state_id" select="./@id"/></xsl:call-template>
</xsl:for-each>
<xsl:call-template name="transition_matrix"/>
<xsl:text>End of MHMM
</xsl:text>
</xsl:template>

<!-- This template prints out the array of background frequencies -->
<xsl:template name="background_frequencies">
  <xsl:for-each select="/MEME_HMM/alphabet/letter">
    <xsl:variable name="alphabet_id"><xsl:value-of select="./@id"/></xsl:variable>
    <xsl:value-of select="/MEME_HMM/background_frequencies/alphabet_array/value[@letter_id=$alphabet_id]"/><xsl:text> </xsl:text>
  </xsl:for-each>
</xsl:template>

<!-- This template prints out the description of a state -->
<xsl:template name="state">
<xsl:param name="state_id"/>
<xsl:text>
</xsl:text>
<xsl:text>State: </xsl:text>
    <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@index"/>
    <xsl:text>
    </xsl:text>
    <xsl:text>type: </xsl:text>
    <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@type"/>
    <xsl:text>
    </xsl:text>
    <xsl:text>i_motif: </xsl:text>
    <xsl:choose>
      <xsl:when test="/MEME_HMM/states/state[@id=$state_id]/@motif_index">
        <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@motif_index"/>
      </xsl:when>
      <xsl:otherwise>-1</xsl:otherwise>
    </xsl:choose>
    <xsl:text>
    </xsl:text>
    <xsl:text>motif_id: </xsl:text>
    <xsl:choose>
      <xsl:when test="/MEME_HMM/states/state[@id=$state_id]/@motif_id">
        <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@motif_id"/>
      </xsl:when>
      <xsl:otherwise>---</xsl:otherwise>
    </xsl:choose>
    <xsl:text>
    </xsl:text>
    <xsl:text>i_position: </xsl:text>
    <xsl:choose>
      <xsl:when test="./@position">
        <xsl:value-of select="./@position"/>
      </xsl:when>
      <xsl:otherwise>-1</xsl:otherwise>
    </xsl:choose>
    <xsl:text>
    </xsl:text>
    <xsl:text>id_char: </xsl:text>
    <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@id_char"/>
    <xsl:text>
    </xsl:text>
    <xsl:text>num_sites: </xsl:text>
    <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/@num_sites"/>
    <xsl:if test="/MEME_HMM/states/state[@id=$state_id]/emission_probabilities">
      <xsl:call-template name="emission_probabilities">
        <xsl:with-param name="state_id" select="./@id"/>
      </xsl:call-template>
    </xsl:if>
</xsl:template>

<!-- This template prints out the array of emission probabilities -->
<xsl:template name="emission_probabilities">
  <xsl:param name="state_id"/>
    <xsl:text>
    </xsl:text>
    <xsl:text>emit: </xsl:text>
    <xsl:for-each select="/MEME_HMM/alphabet/letter">
      <xsl:variable name="alphabet_id"><xsl:value-of select="./@id"/></xsl:variable>
      <xsl:value-of select="/MEME_HMM/states/state[@id=$state_id]/emission_probabilities/alphabet_array/value[@letter_id=$alphabet_id]"/>
      <xsl:text> </xsl:text>
    </xsl:for-each>
</xsl:template>

<!-- This template prints out the matrix of transition probabilities -->
<xsl:template name="transition_matrix">
<xsl:text>
</xsl:text>
<xsl:text>Transition probability matrix (</xsl:text>
<xsl:value-of select="/MEME_HMM/states/@num_states"/> 
<xsl:text> x </xsl:text>
<xsl:value-of select="/MEME_HMM/states/@num_states"/>
<xsl:text>):</xsl:text>
<xsl:text>
</xsl:text>
<xsl:for-each select="/MEME_HMM/states/state">
  <xsl:variable name="state_index"><xsl:value-of select="./@index"/></xsl:variable>
  <xsl:text>    </xsl:text>
  <xsl:for-each select="/MEME_HMM/states/state">
    <xsl:variable name="to_state_index"><xsl:value-of select="./@index"/></xsl:variable>
    <xsl:choose>
      <xsl:when test="/MEME_HMM/states/state[$state_index=@index]/transition[$to_state_index=@to_state_index]">
        <xsl:value-of select="/MEME_HMM/states/state[$state_index=@index]/transition[$to_state_index=@to_state_index]"/>
        <xsl:text> </xsl:text>
      </xsl:when>
      <xsl:otherwise>0 </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
  <xsl:text>
</xsl:text>
</xsl:for-each>
</xsl:template>

</xsl:stylesheet>
