<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE xsl:stylesheet [
<!ENTITY nbsp "&#160;">
<!ENTITY space " ">
<!ENTITY newline "&#10;">
<!ENTITY tab "&#9;">
<!ENTITY more "&#8615;">
]><!-- define nbsp as it is not defined in xml, only lt, gt and amp are defined by default -->
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="html" indent="yes" 
    doctype-public="-//W3C//DTD HTML 4.01 Transitional//EN"
    doctype-system="http://www.w3.org/TR/html4/loose.dtd"
    />
  <!-- contains variable site_url -->
  <xsl:include href="constants.xsl"/>
  <!-- css includes -->
  <xsl:include href="meme.css.xsl"/>
  <xsl:include href="spamo-to-html.css.xsl"/>
  <!-- js includes -->
  <xsl:include href="delay_draw.js.xsl"/>
  <xsl:include href="motif_logo.js.xsl"/>
  <xsl:include href="spamo_graph.js.xsl"/>
  <xsl:include href="spamo-to-html.js.xsl"/>

  <!-- Stylesheet processing starts here -->
  <xsl:template match="/">
    <!-- Possible way to output a html 5 doctype -->
    <!--
    <xsl:text disable-output-escaping='yes'>&lt;!DOCTYPE HTML&gt;&newline;</xsl:text>
    -->
    <html>
      <xsl:call-template name="html-head"/>
      <xsl:call-template name="html-body"/>
    </html>
  </xsl:template>
  
  <xsl:template name="html-head">
    <head>
      <title>SpaMo - Spaced Motif Analysis</title>
      <xsl:call-template name="html-css"/>
      <xsl:call-template name="html-script"/>
    </head>
  </xsl:template>

  <xsl:template name="html-css">
    <style type="text/css">
    <xsl:call-template name="meme.css"/>
    <xsl:call-template name="spamo-to-html.css"/>
    </style>
  </xsl:template>

  <xsl:template name="html-script">
    <script type="text/javascript">
      <xsl:call-template name="delay_draw.js"/>
      <xsl:call-template name="motif_logo.js"/>
      <xsl:call-template name="spamo_graph.js"/>
      <xsl:call-template name="spamo-to-html.js"/>
      <xsl:for-each select="/spamo/primary_motif">
        <xsl:call-template name="store-motif">
          <xsl:with-param name="id" select="position()"/>
          <xsl:with-param name="tabbing" select="'      '"/>
        </xsl:call-template>
      </xsl:for-each>
    </script>
  </xsl:template>

  <xsl:template name="html-body">
    <body onscroll="delayed_process_draw_tasks()" onresize="delayed_process_draw_tasks()">
      <!--<xsl:call-template name="help-popups"/>-->
      <xsl:call-template name="top"/>
      <xsl:call-template name="primary-summary"/>
      <xsl:call-template name="sequence-summary"/>
      <xsl:call-template name="secondary-summary"/>
      <xsl:call-template name="primary-motifs"/>
      <xsl:call-template name="program"/>
      <!--<xsl:call-template name="documentation"/>-->
      <xsl:call-template name="footer"/>
    </body>
  </xsl:template>


  <xsl:template name="top">
    <a name="top"/>
    <div class="pad1"><!-- FIXME url currently local for testing -->
      <h1><img src="{$site_url}/doc/images/spamo_logo.png" alt="Spaced Motif analysis tool (SpaMo)" /></h1>
      <p class="spaced">
        For further information on how to interpret these results or to get a 
        copy of the MEME software please access 
        <a href="http://meme.nbcr.net/">http://meme.nbcr.net</a>. 
      </p>
    </div>
    <!-- navigation -->
    <div class="pad2">
      <a class="jump" href="#primary_summary">Primary Motif</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#sequence_summary">Sequence Database</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#secondary_summary">Secondary Databases</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#spacing_analysis">Spacing Analysis</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#program">Program information</a><!--
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#doc">Explanation</a>-->
    </div>
  </xsl:template>

  <xsl:template name="sequence-summary">
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Sequence Database'"/>
      <xsl:with-param name="self" select="'sequence_summary'"/>
      <xsl:with-param name="next" select="'secondary_summary'"/>
      <xsl:with-param name="prev" select="'primary_summary'"/>
    </xsl:call-template>
    <xsl:variable name="seq" select="/spamo/files/sequence_db"/>
    <div class="box" >
      <table class="preview">
        <tr>
          <th>Name</th>
          <th>Last Modified</th>
          <th>Loaded</th>
          <th>Too Short</th>
          <th>No Primary</th>
          <th>Too Similar</th>
          <th>Used</th>
        </tr>
        <tr>
          <td><xsl:value-of select="$seq/@name"/></td>
          <td><xsl:value-of select="$seq/@last_modified"/></td>
          <td><xsl:value-of select="$seq/@loaded"/></td>
          <td><xsl:value-of select="$seq/@excluded_too_short"/></td>
          <td><xsl:value-of select="$seq/@excluded_no_match"/></td>
          <td><xsl:value-of select="$seq/@excluded_similar"/></td>
          <td><xsl:value-of select="$seq/@loaded - $seq/@excluded_too_short - $seq/@excluded_no_match - $seq/@excluded_similar"/></td>
        </tr>
      </table>
    </div>
  </xsl:template>

  <xsl:template name="primary-summary">
    <xsl:variable name="top-list" select="20"/>
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Primary Motif'"/>
      <xsl:with-param name="self" select="'primary_summary'"/>
      <xsl:with-param name="next" select="'sequence_summary'"/>
    </xsl:call-template>
    <div class="box" >
      <table class="preview">
        <tr>
          <th>Name<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_primary_motif_name')"/>--></th>
          <th>Preview<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_primary_motif_preview')"/>--></th>
          <th>Significant Secondaries<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_primary_motif_sig_motifs')"/>--></th>
          <th>List</th>
        </tr>
        <xsl:for-each select="spamo/primary_motif">
          <xsl:variable name="primary_pos" select="position()"/>
          <tr>
            <td><xsl:call-template name="motif-name"/></td>
            <td>
              <xsl:call-template name="motif-logo">
                <xsl:with-param name="motif_id" select="position()"/>
                <xsl:with-param name="replace_id" select="concat('preview_primary_', position())"/>
                <xsl:with-param name="title_text" select="concat('Primary: ', motif/@name)"/>
                <xsl:with-param name="already_stored" select="'true'"/>
              </xsl:call-template>
            </td>
            <td><xsl:value-of select="count(secondary_motif)"/></td>
            <td>
              <xsl:for-each select="secondary_motif[position() &lt;= $top-list]">
                <xsl:variable name="link" select="concat('s_', $primary_pos, '_', position())"/>
                <a href="#{$link}">
                  <xsl:call-template name="motif-name"/>
                </a>
                <xsl:if test="position() != last() and position() &lt; $top-list">
                  <xsl:text>, </xsl:text>
                </xsl:if>
              </xsl:for-each>
            </td>
          </tr>
        </xsl:for-each>
      </table>
    </div>
  </xsl:template>

  <xsl:template name="secondary-summary">
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Secondary Databases'"/>
      <xsl:with-param name="self" select="'secondary_summary'"/>
      <xsl:with-param name="prev" select="'sequence_summary'"/>
      <xsl:with-param name="next" select="'spacing_analysis'"/>
    </xsl:call-template>
    <div class="box">
      <table class="preview">
        <tr>
          <th>Name<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary_databases_name')"/>--></th>
          <th>Last Modified<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary_databases_last_modified')"/>--></th>
          <th>Number of Motifs<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary_databases_num_motifs')"/>--></th>
          <th>Motifs Significant<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary_databases_sig_motifs')"/>--></th>
          <th>Motifs Redundant<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary_databases_redundant_motifs')"/>--></th>
        </tr>
        <xsl:for-each select="/spamo/files/motif_db[@id != 'primary_file']">
          <xsl:variable name="db_id" select="@id"/>
          <tr>
            <td><xsl:value-of select="@name"/></td>
            <td><xsl:value-of select="@last_modified"/></td>
            <td><xsl:value-of select="@loaded - @excluded"/></td>
            <td><xsl:value-of select="count(/spamo/primary_motif/secondary_motif/motif[@db = $db_id])"/></td>
            <td><xsl:value-of select="count(/spamo/primary_motif/secondary_motif/redundant/secondary_motif/motif[@db = $db_id])"/></td>
          </tr>
        </xsl:for-each>
      </table>
    </div>
  </xsl:template>

  <xsl:template name="primary-motifs">
      <a name="spacing_analysis"/>
      <xsl:for-each select="/spamo/primary_motif">
        <xsl:variable name="primary_pos" select="position()"/>
        <xsl:variable name="primary_last" select="position() = last()"/>
          <xsl:for-each select="secondary_motif">
            <xsl:choose>
              <xsl:when test="$primary_last and position() = last()">
               <a name="last_spacing_analysis"/>
             </xsl:when>
             <xsl:when test="position() = last()">
               <a name="l_{$primary_pos}"/>
             </xsl:when>
            </xsl:choose>
            <xsl:call-template name="secondary-motif">
              <xsl:with-param name="primary_index" select="$primary_pos"/>
              <xsl:with-param name="id" select="concat('s_', $primary_pos, '_', position())"/>
              <xsl:with-param name="prev">
                <xsl:choose>
                  <xsl:when test="$primary_pos = 1 and position() = 1">secondary_summary</xsl:when>
                  <xsl:when test="position() = 1">
                    <xsl:value-of select="concat('l_', $primary_pos - 1)"/>
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="concat('s_', $primary_pos, '_', position() - 1)"/>
                  </xsl:otherwise>
                </xsl:choose>
              </xsl:with-param>
              <xsl:with-param name="next">
                <xsl:choose>
                  <xsl:when test="$primary_last and position() = last()">
                    <xsl:text>program</xsl:text>
                  </xsl:when>
                  <xsl:when test="position() = last()">
                    <xsl:value-of select="concat('s_', $primary_pos + 1, '_1')"/>
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="concat('s_', $primary_pos, '_', position() + 1)"/>
                  </xsl:otherwise>
                </xsl:choose>
              </xsl:with-param>
            </xsl:call-template>
          </xsl:for-each>
      </xsl:for-each>
  </xsl:template>

  <xsl:template name="secondary-motif">
    <xsl:param name="primary_index"/>
    <xsl:param name="id"/>
    <xsl:param name="prev"/>
    <xsl:param name="next"/>


    <xsl:variable name="title">
      <xsl:text>Spacings of &quot;</xsl:text>
      <xsl:call-template name="motif-name"/>
      <xsl:text>&quot; relative to &quot;</xsl:text>
      <xsl:call-template name="motif-name"><xsl:with-param name="motif" select="../motif"/></xsl:call-template>
      <xsl:text>&quot;</xsl:text>
    </xsl:variable>

    <xsl:call-template name="header">
      <xsl:with-param name="title" select="$title"/>
      <xsl:with-param name="self" select="$id"/>
      <xsl:with-param name="prev" select="$prev"/>
      <xsl:with-param name="next" select="$next"/>
    </xsl:call-template>
    <div class="box" >
      <!-- table showing primary and secondary motifs -->
      <table>
        <tr>
          <th>Primary: <xsl:call-template name="linked-motif-name">
              <xsl:with-param name="motif" select="../motif"/>
            </xsl:call-template><!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_primary')"/>-->
          </th>
          <th>Secondary: <xsl:call-template name="linked-motif-name"/><!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary')"/>--></th>
        </tr>
        <tr>
          <td>
            <xsl:call-template name="motif-logo">
              <xsl:with-param name="motif_id" select="$primary_index"/>
              <xsl:with-param name="replace_id" select="concat($id, '_pm')"/>
              <xsl:with-param name="motif" select="../motif"/>
              <xsl:with-param name="title_text" select="concat('Primary: ', ../motif/@name)"/>
              <xsl:with-param name="already_stored" select="'true'"/>
            </xsl:call-template>
          </td>
          <td>
            <xsl:call-template name="motif-logo">
              <xsl:with-param name="motif_id" select="$id"/>
              <xsl:with-param name="replace_id" select="concat($id, '_sm')"/>
              <xsl:with-param name="title_text" select="concat('Secondary: ', motif/@name)"/>
            </xsl:call-template>
          </td>
        </tr>
      </table>
      <table>
        <tr>
          <th colspan="2">Motif Spacing Histogram<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_histogram')"/>--></th>
          <th>&nbsp;</th><th colspan="2">Significant Motif Spacings<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_spacings')"/>--></th><th>&nbsp;</th>
        </tr>
        <tr>
          <th style="width:175px">Upstream</th><th style="width:175px">Downstream</th>
          <th>&nbsp;</th>
          <th style="width:130px">Upstream</th><th style="width:130px">Downstream</th>
          <th>Other Details</th>
        </tr>
        <tr>
          <td colspan="2">
            <xsl:call-template name="secondary-motif-histogram">
              <xsl:with-param name="graph_id" select="concat($id, '_hist')"/>
            </xsl:call-template>
          </td>
          <th>
            <xsl:call-template name="strand-titles"/>
          </th>
          <td colspan="2">
            <xsl:call-template name="secondary-motif-spacings"/>
          </td>
          <td class="details" style="vertical-align:top">
            <h4>Total sequences with primary and secondary motif<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_seqtotal')"/>--></h4>
            <xsl:value-of select="histogram/@total"/><br/>
            <h4>Motif Database<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_motifdb')"/>--></h4><!-- TODO add popup -->
            <xsl:value-of select="id(motif/@db)/@name"/><br/>
          <xsl:if test="redundant">
          <fieldset class="sm_group">
            <h4>Secondary motifs with similar spacings<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_redundant')"/>--></h4>
            <xsl:variable name="ck_scroll">
              <xsl:choose>
                <xsl:when test="count(redundant/secondary_motif) &gt; 6">ck_scroll</xsl:when>
                <xsl:otherwise>ck_fixed</xsl:otherwise>
              </xsl:choose>
            </xsl:variable>
            <div class="{$ck_scroll}">
            <xsl:for-each select="redundant/secondary_motif">
            <label><input type="checkbox" autocomplete="off" onclick="toggle_group(&quot;{concat($id, '_', position())}&quot;)"/>
              <xsl:call-template name="motif-name"/>
            </label>
            </xsl:for-each>
            </div>
          </fieldset>
          </xsl:if>
          </td>
        </tr>
      </table>
      <xsl:for-each select="redundant/secondary_motif">
      <xsl:variable name="sid" select="concat($id, '_', position())"/>
      <div id="{$sid}">
        <table>
          <tr>
            <th colspan="6"><xsl:text>Similar Secondary: </xsl:text><xsl:call-template name="linked-motif-name"/></th>
          </tr>
          <tr>
            <td colspan="2" style="width:350px;">
              <xsl:call-template name="secondary-motif-histogram">
                <xsl:with-param name="graph_id" select="concat($sid, '_hist')"/>
                <xsl:with-param name="img" select="concat('hist_', ../../../motif/@name, '_', motif/@db, '_', motif/@name, '.png')"/>
                <xsl:with-param name="group_id" select="$sid"/>
              </xsl:call-template>
            </td>
            <th>
              <xsl:call-template name="strand-titles"/>
            </th>
            <td colspan="2">
              <xsl:call-template name="secondary-motif-spacings"/>
            </td>
            <td class="details" style="vertical-align:top">
              <h4>Total sequences with primary and secondary motif<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_seqtotal')"/>--></h4>
            <xsl:value-of select="histogram/@total"/>
            <xsl:choose>
              <xsl:when test="/spamo/model/bin_size = 1">
                <h4>Alignment by most significant spacings<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_alignment')"/>--></h4>
                  <xsl:call-template name="secondary-motif-alignment">
                    <xsl:with-param name="s1_id" select="$id"/>
                    <xsl:with-param name="s2_id" select="$sid"/>
                  </xsl:call-template>
              </xsl:when>
              <xsl:otherwise>
                <h4>Secondary Motif<!--&nbsp;<div class="help" onclick="help_popup(this, 'pop_secondary')"/>--></h4>
                <xsl:call-template name="motif-logo">
                  <xsl:with-param name="motif_id" select="$sid"/>
                  <xsl:with-param name="replace_id" select="concat($sid, '_ssm')"/>
                  <xsl:with-param name="group_id" select="$sid"/>
                  <xsl:with-param name="title_text" select="concat('Similar Secondary: ', motif/@name)"/>
                </xsl:call-template>
              </xsl:otherwise>
            </xsl:choose>
            </td>
          </tr>
        </table>
      </div>
      <script>
        document.getElementById(&quot;<xsl:value-of select="$sid"/>&quot;).style.display = "none";
      </script>
      </xsl:for-each>
    </div>
  </xsl:template>

  <xsl:template name="motif-name">
    <xsl:param name="motif" select="motif"/>
    <xsl:value-of select="translate($motif/@name, '_', ' ')"/>
    <xsl:if test="$motif/@alt">
      <xsl:value-of select="concat(' (', translate($motif/@alt, '_', ' '), ')')"/>
    </xsl:if>
  </xsl:template>

  <xsl:template name="linked-motif-name">
    <xsl:param name="motif" select="motif"/>
    <xsl:choose>
      <xsl:when test="$motif/@url">
        <a href="{$motif/@url}"><xsl:call-template name="motif-name"><xsl:with-param name="motif" select="$motif"/></xsl:call-template></a>
      </xsl:when>
      <xsl:otherwise>
        <xsl:call-template name="motif-name"><xsl:with-param name="motif" select="$motif"/></xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="strand-titles">
    <div class="st_box" >
      <div class="st_same deg90rotate">Same Strand</div>
      <div class="st_oppo deg90rotate">Opposite Strand</div>
    </div>
  </xsl:template>

<!--
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Creates a span containing the consensus sequence of the secondary logo which
is replaced by a canvas based logo when possible.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-->
  <xsl:template name="motif-logo">
    <xsl:param name="motif_id"/>
    <xsl:param name="replace_id" select="$motif_id"/>
    <xsl:param name="motif" select="motif"/>
    <xsl:param name="group_id" select="''"/>
    <xsl:param name="title_text" select="''"/>
    <xsl:param name="already_stored" select="'false'"/>

    <xsl:call-template name="motif-consensus">
      <xsl:with-param name="motif" select="$motif"/>
      <xsl:with-param name="id" select="$replace_id"/>
    </xsl:call-template>
    <script>
      <xsl:if test="$already_stored = 'false'">
        <xsl:call-template name="store-motif">
          <xsl:with-param name="id" select="$motif_id"/>
        </xsl:call-template>
      </xsl:if>
      <xsl:call-template name="ready-motif">
        <xsl:with-param name="replace_id" select="$replace_id"/> 
        <xsl:with-param name="motif_id" select="$motif_id"/>
        <xsl:with-param name="title_text" select="$title_text"/>
        <xsl:with-param name="group_id" select="$group_id"/>
      </xsl:call-template>
    </script>
  </xsl:template>

<!--
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Creates a table with the main secondary motif and similar secondary motif 
aligned. If canvas or javascript is not avaliable then it uses an aligned
consensus.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-->
  <xsl:template name="secondary-motif-alignment">
    <xsl:param name="s1_id"/>
    <xsl:param name="s2_id"/>

    <!-- calculate indents -->
    <xsl:variable name="sp1" select="../../spacing[1]"/>
    <xsl:variable name="mot1" select="../../motif"/>
    <xsl:variable name="sp2" select="spacing[1]"/>
    <xsl:variable name="mot2" select="motif"/>

    <!--
    -->
    <xsl:variable name="s1_strand">
      <xsl:choose>
        <xsl:when test="$sp2/@strand = $sp1/@strand">
          <xsl:text>false</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>true</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <!-- 
    calculate the indent given that the motifs are oriented as per the best peak's orientation
    then later we can use this to orient them so the bottom motif is always displayed normal orientation 
    note this assumes a bin size of 1
    -->
    <xsl:variable name="natural_indent">
      <!-- calculate the trimmed length -->
      <xsl:variable name="mot1_tlen" select="$mot1/@length - $mot1/@ltrim - $mot1/@rtrim"/>
      <xsl:variable name="mot2_tlen" select="$mot2/@length - $mot2/@ltrim - $mot2/@rtrim"/>
      <xsl:variable name="mot1_trim">
        <xsl:choose>
          <xsl:when test="$sp1/@strand = 'same'">
            <xsl:value-of select="$mot1/@ltrim"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$mot1/@rtrim"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:variable name="mot2_trim">
        <xsl:choose>
          <xsl:when test="$sp2/@strand = 'same'">
            <xsl:value-of select="$mot2/@ltrim"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$mot2/@rtrim"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:variable>
      <xsl:choose>
        <xsl:when test="$sp1/@side = 'left'"><!-- both are upstream -->
          <xsl:variable name="sp1_offset" select="$sp1/@bin + $mot1_tlen + $mot1_trim"/>
          <xsl:variable name="sp2_offset" select="$sp2/@bin + $mot2_tlen + $mot2_trim"/>
          <xsl:value-of select="$sp1_offset - $sp2_offset"/>
        </xsl:when>
        <xsl:otherwise><!-- both are downstream -->
          <xsl:variable name="sp1_offset" select="$sp1/@bin - $mot1_trim"/>
          <xsl:variable name="sp2_offset" select="$sp2/@bin - $mot2_trim"/>
          <xsl:value-of select="$sp2_offset - $sp1_offset"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <!-- 
    calculate the indent required when keeping the redundant motif at the correct orientation
    -->
    <xsl:variable name="corrected_indent">
      <xsl:choose>
        <xsl:when test="$sp2/@strand = 'same'"><!-- alignment can be left as natural -->
          <xsl:value-of select="$natural_indent"/>
        </xsl:when>
        <xsl:otherwise><!-- alignment must be flipped -->
          <xsl:value-of select="$mot1/@length - $mot2/@length - $natural_indent"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <xsl:variable name="s1_indent">
      <xsl:choose>
        <xsl:when test="$corrected_indent &lt; 0"><xsl:value-of select="0 - $corrected_indent"/></xsl:when>
        <xsl:otherwise>0</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="s2_indent">
      <xsl:choose>
        <xsl:when test="$corrected_indent &gt; 0"><xsl:value-of select="$corrected_indent"/></xsl:when>
        <xsl:otherwise>0</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

    <table>
      <tr>
        <td>Best Similar<br/>Secondary</td>
        <td>
          <xsl:call-template name="motif-consensus">
            <xsl:with-param name="motif" select="../../motif"/>
            <xsl:with-param name="id" select="concat($s2_id, '_sm')"/>
            <xsl:with-param name="col" select="$s1_indent"/>
            <xsl:with-param name="rc" select="$s1_strand"/>
            <xsl:with-param name="box" select="'false'"/>
          </xsl:call-template>
        </td>
      </tr>
      <tr>
        <td>This Similar<br/>Secondary</td>
        <td>
          <xsl:call-template name="motif-consensus">
            <xsl:with-param name="motif" select="motif"/>
            <xsl:with-param name="id" select="concat($s2_id, '_ssm')"/>
            <xsl:with-param name="col" select="$s2_indent"/>
            <xsl:with-param name="box" select="'false'"/>
          </xsl:call-template>
        </td>
      </tr>
    </table>
    <script>
      <xsl:call-template name="store-motif">
        <xsl:with-param name="id" select="$s2_id"/>
      </xsl:call-template>
      
      <xsl:call-template name="ready-motif">
        <xsl:with-param name="replace_id" select="concat($s2_id, '_sm')"/> 
        <xsl:with-param name="motif_id" select="$s1_id"/>
        <xsl:with-param name="title_text" select="concat('Secondary: ', $mot1/@name)"/>
        <xsl:with-param name="group_id" select="$s2_id"/>
        <xsl:with-param name="rc" select="$s1_strand"/>
        <xsl:with-param name="indent" select="$s1_indent"/>
      </xsl:call-template>

      <xsl:call-template name="ready-motif">
        <xsl:with-param name="replace_id" select="concat($s2_id, '_ssm')"/> 
        <xsl:with-param name="motif_id" select="$s2_id"/>
        <xsl:with-param name="title_text" select="concat('Similar Secondary: ', $mot2/@name)"/>
        <xsl:with-param name="group_id" select="$s2_id"/>
        <xsl:with-param name="indent" select="$s2_indent"/>
      </xsl:call-template>

    </script>
  </xsl:template>

  <xsl:template name="secondary-motif-histogram">
    <xsl:param name="graph_id"/>
    <xsl:param name="replace_id" select="$graph_id"/>
    <xsl:param name="img" select="concat('hist_', ../motif/@name, '_', motif/@db, '_', motif/@name, '.png')"/>
    <xsl:param name="group_id" select="''"/>

    <noscript><img src="{$img}"/></noscript>
    <span id="{$replace_id}"></span>
    <script type="text/javascript">
      store_graph(
        &quot;<xsl:value-of select="$graph_id"/>&quot;,
        <xsl:call-template name="spamo-graph">
          <xsl:with-param name="smotif" select="."/>
          <xsl:with-param name="indent" select="'        '"/>
        </xsl:call-template>
      );
      ready_graph(
        &quot;<xsl:value-of select="$replace_id"/>&quot;,
        &quot;<xsl:value-of select="$graph_id"/>&quot;,
        false,
        &quot;<xsl:value-of select="$img"/>&quot;,
        &quot;&quot;
        <xsl:if test="$group_id != ''">
          ,&quot;<xsl:value-of select="$group_id"/>&quot;
        </xsl:if>
      );
    </script>
  </xsl:template>

  <xsl:template name="secondary-motif-spacings">
    <div class="sp_box">
      <div class="sp_same_up"> 
          <xsl:call-template name="secondary-motif-quad-spacings">
            <xsl:with-param name="strand" select="'same'"/>
            <xsl:with-param name="side" select="'left'"/>
          </xsl:call-template>
      </div>
      <div class="sp_same_down">
          <xsl:call-template name="secondary-motif-quad-spacings">
            <xsl:with-param name="strand" select="'same'"/>
            <xsl:with-param name="side" select="'right'"/>
          </xsl:call-template>
      </div>
      <div class="sp_oppo_up">
          <xsl:call-template name="secondary-motif-quad-spacings">
            <xsl:with-param name="strand" select="'opposite'"/>
            <xsl:with-param name="side" select="'left'"/>
          </xsl:call-template>
      </div>
      <div class="sp_oppo_down">
          <xsl:call-template name="secondary-motif-quad-spacings">
            <xsl:with-param name="strand" select="'opposite'"/>
            <xsl:with-param name="side" select="'right'"/>
          </xsl:call-template>
      </div>
    </div>
  </xsl:template>


  <xsl:template name="secondary-motif-quad-spacings">
    <xsl:param name="strand"/>
    <xsl:param name="side"/>
    <xsl:variable name="spacing_set" select="spacing[@strand = $strand and @side = $side]"/>
    <xsl:variable name="spacing_count" select="count($spacing_set)"/>
    <xsl:variable name="sp_scroll">
      <xsl:choose>
        <xsl:when test="$spacing_count &gt; 5">sp_scroll</xsl:when>
        <xsl:otherwise>sp_fixed</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:if test="$spacing_count &gt; 0"> 
    <div class="{$sp_scroll}">
      <table >
        <thead>
          <tr>
            <th>P-value</th>
            <th>Gap</th>
            <th>#</th>
            <th style="width:15px;">&nbsp;</th>
          </tr>
        </thead>
        <tbody>
          <xsl:for-each select="$spacing_set">
            <xsl:call-template name="secondary-motif-spacing"/>
          </xsl:for-each>
        </tbody>
      </table>
    </div>
    </xsl:if>
  </xsl:template>

  <xsl:template name="secondary-motif-spacing">
    <tr>
      <td class="sp_pvalue" ><xsl:value-of select="@pvalue"/></td>
      <td class="sp_gap">
        <xsl:choose>
          <xsl:when test="/spamo/model/bin_size = 1">
            <xsl:value-of select="@bin - 1"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:variable name="high" select="@bin * /spamo/model/bin_size - 1"/>
            <xsl:variable name="low" select="$high - /spamo/model/bin_size + 1"/>
            <xsl:value-of select="$low"/>
            <xsl:text>:</xsl:text>
            <xsl:value-of select="$high"/>
          </xsl:otherwise>
        </xsl:choose>
      </td>
      <td class="sp_count"><xsl:value-of select="@num"/></td>
      <td>&nbsp;</td>
    </tr>
  </xsl:template>

  <xsl:template name="program">
    <a name="program"/>
    <div class="bar">
      <div style="text-align:right;"><a href="#last_spacing_analysis">Previous</a><!--<xsl:text> </xsl:text>-->
        <!--<a href="#doc">Next</a>--><xsl:text> </xsl:text><a href="#top">Top</a></div>
      <div class="subsection">
        <a name="version"/>
        <h5>SpaMo version</h5>
        <xsl:value-of select="/spamo/@version"/> (Release date: <xsl:value-of select="/spamo/@release"/>)
      </div>
      <div class="subsection">
        <a name="reference"/>
        <h5>Reference</h5>
        "Inferring transcription factor complexes from ChIP-seq data", Tom Whitington, Martin C. Frith, James Johnson, and Timothy L. Bailey,
        submitted for publication.        
      </div>
      <div class="subsection">
        <a name="command" />
        <h5>Command line summary</h5>
        <textarea rows="1" style="width:100%;" readonly="readonly">
          <xsl:value-of select="/spamo/model/command_line"/>
        </textarea>
        <br />
        <xsl:text>Result calculation took </xsl:text>
        <xsl:variable name="total_seconds" select="/spamo/run_time/@real"/>
        <xsl:variable name="total_minutes" select="floor($total_seconds div 60)"/>
        <xsl:variable name="hours" select="floor($total_minutes div 60)"/>
        <xsl:variable name="minutes" select="$total_minutes - 60 * $hours"/>
        <!-- don't use mod here as seconds value has fractions of a second too -->
        <xsl:variable name="seconds" select="$total_seconds - 60 * $total_minutes"/>
        <xsl:if test="$hours &gt; 0">
          <xsl:value-of select="$hours"/> hour<xsl:if test="$hours != 1">s</xsl:if>
          <xsl:text> </xsl:text>
        </xsl:if>
        <xsl:if test="$hours &gt; 0 or $minutes &gt; 0">
          <xsl:value-of select="$minutes"/> minute<xsl:if test="$minutes != 1">s</xsl:if>
          <xsl:text> </xsl:text>
        </xsl:if>
        <xsl:value-of select="$seconds"/> second<xsl:if test="$seconds != 1">s</xsl:if>
        <br />
        <xsl:text>Note that the random number generator was initilized with a seed of </xsl:text>
        <xsl:value-of select="/spamo/model/seed"/> so you need
        &quot;-numgen <xsl:value-of select="/spamo/model/seed"/>&quot; in the list of arguments
        to replicate the experiment.
      </div>      
      <a href="javascript:show_hidden('model')" id="model_activator">show model parameters...</a>
      <div class="subsection" id="model_data" style="display:none;">
        <h5>Model parameters</h5>
        <xsl:text>&newline;</xsl:text>
        <textarea style="width:100%;" rows="{count(/spamo/model/*) - 1}" readonly="readonly">
          <xsl:variable name="spaces" select="'                    '"/>
          <xsl:text>&newline;</xsl:text>
          <xsl:for-each select="/spamo/model/*[name(.) != 'command_line']">
            <xsl:variable name="pad" select="substring($spaces,string-length(name(.)))"/>
            <xsl:value-of select="name(.)"/>
            <xsl:value-of select="$pad"/>
            <xsl:text> = </xsl:text>
            <xsl:choose>
              <xsl:when test="count(@*) &gt; 0">
                <xsl:for-each select="@*">
                  <xsl:value-of select="name(.)"/>
                  <xsl:text>: "</xsl:text>
                  <xsl:value-of select="."/>
                  <xsl:text>"</xsl:text>
                  <xsl:if test="position() != last()">
                    <xsl:text>, </xsl:text>
                  </xsl:if>
                </xsl:for-each>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="normalize-space(.)"/>
              </xsl:otherwise>
            </xsl:choose>
            <xsl:text>&newline;</xsl:text>
          </xsl:for-each>
        </textarea>
      </div>
      <a href="javascript:hide_shown('model')" style="display:none;" id="model_deactivator">hide model parameters...</a>
    </div>
  </xsl:template>

  <xsl:template name="help-popups">
    <div class="pop_content" id="pop_primary_motif_name">
      <p>The name of the primary motif.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_primary_motif_preview">
      <p>The logo of the primary motif.</p>
      <p>Sections of the motif with a gray background have been trimmed and were not used for scanning.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_primary_motif_sig_motifs">
      <p>The number of secondary motifs found that had significant spacings in the tested region.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary_databases_name">
      <p>The name of the motif database derived from the file name.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary_databases_last_modified">
      <p>The date that the motif database was last modified.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary_databases_num_motifs">
      <p>The number of motifs loaded from the motif database. Some motifs may have been excluded.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary_databases_sig_motifs">
      <p>The number of motifs with significant spacings that were not considered too similar to another motif.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary_databases_redundant_motifs">
      <p>The count of motifs that while having significant spacings were less signficant than another motif that matched
      most of the same sites.</p>
      <div style="float:right; bottom:0px;">[ <a href="#">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_primary">
      <p>The primary motif is used as the reference point for all spacing calculation.</p>
      <p>Sections of the motif with a gray background have been trimmed and were not used for scanning.</p>
      <div style="float:right; bottom:0px;">[ <a href="#doc_primary">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_secondary">
      <p>The secondary motif occurs at the spacings relative to the primary shown in the histogram below.</p>
      <p>Sections of the motif with a gray background have been trimmed and were not used for scanning.</p>
      <div style="float:right; bottom:0px;">[ <a href="#doc_secondary">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_histogram">
      <p>The histogram below shows the frequency of spacings from the primary motif to the secondary motif. Red bars
      indicate that a spacing has occured a statistically significant number of times.</p>
      <dl>
      <dt><b>Upstream</b></dt>
      <dd>These are sequences where the secondary motif occurs before the primary motif.</dd>
      <dt><b>Downstream</b></dt>
      <dd>These are sequences where the secondary motif occurs after the primary motif.</dd>
      <dt><b>Same Strand</b></dt>
      <dd>These are sequences where the secondary motif is on the same strand as the primary motif.</dd>
      <dt><b>Opposite Strand</b></dt>
      <dd>These are sequences where the secondary motif is on the opposite strand to the primary motif.</dd>
      </dl>
      <div style="float:right; bottom:0px;">[ <a href="#doc_histogram">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_spacings">
      <p>The details of the significant spacings are shown in the four quadrants as they relate to the quadrants of the graph.</p>
      <dl>
      <dt><b>P-value</b></dt> 
      <dd>is the probability that the spacing happened by chance.</dd> 
      <dt><b>Gap</b></dt> 
      <dd>is the space between the primary and secondary motifs where a value of zero means there is no space between them.</dd> 
      <dt><b>#</b></dt> 
      <dd>is the count of sequences where that spacing was observed.</dd>
      </dl>
      <div style="float:right; bottom:0px;">[ <a href="#doc_spacings">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_seqtotal">
      <p>The total count of sequences that have a match for both the primary motif and this secondary motif.</p>
      <div style="float:right; bottom:0px;">[ <a href="#doc_seqtotal">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_redundant">
      <p>The list of secondary motifs which have been identified as possibly similar to this motif because their most significant spacings 
      overlaped with the most significant spacing of this motif. Check the boxes to show the possibly similar secondary motifs.</p>
      <div style="float:right; bottom:0px;">[ <a href="#doc_redundant">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>
    <div class="pop_content" id="pop_alignment">
      <p>This shows the first secondary motif aligned with the current secondary motif by the most significant spacing.</p>
      <p>Sections of the motif with a gray background have been trimmed and were not used for scanning.</p>
      <div style="float:right; bottom:0px;">[ <a href="#doc_alignment">more</a> | <a href="javascript:help_popup()">close</a> ]</div>
    </div>

  </xsl:template>

  <xsl:template name="documentation">
    <span class="explain">
      <xsl:call-template name="header">
        <xsl:with-param name="title" select="'Explanation of SpaMo Results'"/>
        <xsl:with-param name="self" select="'doc'"/>
        <xsl:with-param name="prev" select="'program'"/>
      </xsl:call-template>
    
      <div class="box">
        <h4>The SpaMo results consist of</h4>
        <ul>
          <li>The <a href="#doc_spacing_analysis">motif spacing analysis</a>. [<a href="#spacing_analysis">View</a>]</li>
          <li>The <b>program</b> details including:
            <ol>
              <li>The <b>version</b> of SpaMo and the date it was released. [<a href="#version">View</a>]</li>
              <li>The <b>reference</b> to cite if you use SpaMo in your research. [<a href="#reference">View</a>]</li>
              <li>The <b>command line summary</b> detailing the parameters with which you ran SpaMo. [<a href="#command">View</a>]</li>
            </ol>
          </li>
          <li>This <b>explanation</b> of how to interpret SpaMo results.</li>
        </ul>
        <a name="doc_spacing_analysis"/>
        <h4>Spacing Analysis</h4>
        <p>All secondary motifs which were found to have a significant spacing relative to the primary motif are listed in order of increasing p-value.</p>
        <a name="doc_primary"/>
        <h5>Primary Motif</h5>
        <div class="doc">
          <p>The primary motif is the motif for which possible interacting motifs are discovered. 
          All spacings use as their reference point the best match to the primary motif in the 
          sequence which leaves a margin on either side which the interacting motifs are discovered.</p>
        </div>
        <a name="doc_secondary"/>
        <h5>Secondary Motif</h5>
        <div class="doc">
          <p>The secondary motif is a candidate for interaction with the primary motif as to be 
          listed it must have a significant spacing relative to the primary motif.</p>
        </div>
        <a name="doc_histogram"/>
        <h5>Motif Spacing Histogram</h5>
        <div class="doc">
          <p>The histogram shows the frequency of spacings of the primary to the secondary motif. 
            Any spacing with a statistically significant count is highlighted in red.</p>
          <p>The histogram is divided into four quadrants which give the position and orientation
            of the secondary motif relative to the primary motif.</p>
        </div>
        <a name="doc_spacings"/>
        <h5>Significant Motif Spacings</h5>
        <div class="doc">
          <p>The significant spacings in the histogram are shown divided up by the quadrant of the graph. Within each quadrant the spacings are sorted 
          in order of significance.</p>
          <dl>
            <dt>P-value</dt>
            <dd>The probability that this peak occurred due to chance.</dd>
            <dt>Gap</dt>
            <dd>The number of bases separating the primary and secondary motif. A gap of zero means the two motifs are up against each other.</dd>
            <dt>#</dt>
            <dd>The number of sequences that had the best primary motif match separated from the best 
              secondary motif match in the margin by exactly Gap bases.</dd>
          </dl>
        </div>
        <a name="doc_seqtotal"/>
        <h5>Total sequences with primary and secondary motif</h5>
        <div class="doc">
          <p>The number of sequences shown on the graph. The graph shows only the sequences which had 
            binding sites for both the primary and secondary motifs.</p>
        </div>
        <h4>Secondary motifs with similar spacings</h4>
        <p>These are motifs which have overlapping best primary peaks.</p>
<!--
        <a name="anchor_goes_here"/>
        <h5>Subtitle</h5>
        <div class="doc">
          <p>Description goes here</p>
          <dl>
            <a name="item_anchor"/>
            <dt>Item Name</dt>
            <dd>The item description.</dd>
          </dl>
        </div>
-->
<!--
        <h4>Title</h4>
        <p>Description goes here.</p>
        <a name="anchor_goes_here"/>
        <h5>Subtitle</h5>
        <div class="doc">
          <p>Description goes here</p>
          <dl>
            <a name="item_anchor"/>
            <dt>Item Name</dt>
            <dd>The item description.</dd>
          </dl>
        </div>
-->
      </div>
    </span>
  </xsl:template>

  <xsl:template name="footer">
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
    <br/>
  </xsl:template>


  <xsl:template name="store-motif">
    <xsl:param name="id"/>
    <xsl:param name="motif" select="motif"/>
    <xsl:param name="tabbing" select="''"/>

    <xsl:variable name="dtabbing" select="concat($tabbing, '  ')"/>

    <xsl:text>&newline;</xsl:text><xsl:value-of select="$tabbing"/>
    <xsl:text>store_motif(&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:text>&quot;</xsl:text><xsl:value-of select="$id"/><xsl:text>&quot;,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:text>&quot;</xsl:text><xsl:value-of select="motif/@name"/><xsl:text>&quot;,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:value-of select="motif/@ltrim"/><xsl:text>,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:value-of select="motif/@rtrim"/><xsl:text>,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:call-template name="motif-pspm">
      <xsl:with-param name="motif" select="motif"/>
      <xsl:with-param name="indent" select="$dtabbing"/>
    </xsl:call-template>
    <xsl:text>&newline;</xsl:text>
    <xsl:value-of select="$tabbing"/>
    <xsl:text>);&newline;</xsl:text>
  </xsl:template>

  <xsl:template name="ready-motif">
    <xsl:param name="replace_id"/>
    <xsl:param name="motif_id"/>
    <xsl:param name="cache" select="'false'"/>
    <xsl:param name="indent" select="0"/>
    <xsl:param name="rc" select="'false'"/>
    <xsl:param name="title_text" select="''"/>
    <xsl:param name="group_id" select="''"/>
    <xsl:param name="tabbing" select="''"/>

    <xsl:variable name="dtabbing" select="concat($tabbing, '  ')"/>

    <xsl:text>&newline;</xsl:text><xsl:value-of select="$tabbing"/>
    <xsl:text>ready_motif(&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:text>&quot;</xsl:text><xsl:value-of select="$replace_id"/><xsl:text>&quot;,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:text>&quot;</xsl:text><xsl:value-of select="$motif_id"/><xsl:text>&quot;,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:value-of select="$cache"/><xsl:text>,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:value-of select="$indent"/><xsl:text>,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:value-of select="$rc"/><xsl:text>,&newline;</xsl:text>
    <xsl:value-of select="$dtabbing"/>
    <xsl:text>&quot;</xsl:text><xsl:value-of select="$title_text"/><xsl:text>&quot;</xsl:text>
    <xsl:if test="$group_id != ''">
      <xsl:text>,&newline;</xsl:text>
      <xsl:value-of select="$dtabbing"/>
      <xsl:text>&quot;</xsl:text>
      <xsl:value-of select="$group_id"/>
      <xsl:text>&quot;&newline;</xsl:text>
    </xsl:if>
    <xsl:text>&newline;</xsl:text>
    <xsl:value-of select="$tabbing"/>
    <xsl:text>);&newline;</xsl:text>
  </xsl:template>


  <xsl:template name="motif-pspm">
    <xsl:param name="motif"/>
    <xsl:param name="indent" select="''"/>

    <xsl:text>&quot;letter-probability matrix: alength= 4 w= </xsl:text><xsl:value-of select="$motif/@length"/>
    <xsl:text> nsites= </xsl:text><xsl:value-of select="$motif/@nsites"/><xsl:text> E= </xsl:text>
    <xsl:value-of select="$motif/@evalue"/><xsl:text>\n&quot;</xsl:text>
    <xsl:for-each select="$motif/pos">
      <xsl:variable name="psum" select="@A + @C + @G + @T"/>
      <xsl:variable name="delta" select="0.01"/>
      <xsl:if test="($psum &gt; 1 and $psum - 1 &gt; $delta) or ($psum &lt; 1 and 1 - $psum &gt; $delta)">
        <xsl:message>
          Warning: Motif <xsl:value-of select="@id"/> probabilities at pos <xsl:value-of select="@i"/> don't sum to 1 
          (delta <xsl:value-of select="$delta"/>).
        </xsl:message>
      </xsl:if>
      <xsl:text> + &newline;</xsl:text>
      <xsl:value-of select="$indent"/><xsl:text>&quot;</xsl:text>
      <xsl:value-of select="@A"/><xsl:text> </xsl:text><xsl:value-of select="@C"/><xsl:text> </xsl:text>
      <xsl:value-of select="@G"/><xsl:text> </xsl:text><xsl:value-of select="@T"/>
      <xsl:text>\n&quot;</xsl:text>
    </xsl:for-each>
  </xsl:template>


  <xsl:template name="spamo-graph-bins">
    <xsl:param name="bins" />
    <xsl:param name="indent" select="''"/>

    <xsl:value-of select="$indent"/><xsl:text>&quot;</xsl:text>
    <xsl:for-each select="$bins">
      <xsl:choose>
        <xsl:when test="@p">
          <xsl:value-of select="concat(@n, ' ', @p)"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="concat(@n, ' 1')"/>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:choose>
        <xsl:when test="position() = last()">
          <xsl:text>\n\n&quot;</xsl:text>
        </xsl:when>
        <xsl:when test="position() mod 40 = 0">
          <xsl:text>\n&quot; + &newline;</xsl:text>
          <xsl:value-of select="$indent"/>
          <xsl:text>&quot;</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text> </xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="spamo-graph">
    <xsl:param name="smotif"/>
    <xsl:param name="indent" select="''"/>
    
    <xsl:text>&quot;% Input Data\n&quot; + &newline;</xsl:text>
    <xsl:value-of select="$indent"/>
    <xsl:text>&quot;</xsl:text>
    <xsl:value-of select="/spamo/model/margin"/><xsl:text> margin </xsl:text>
    <xsl:value-of select="$smotif/motif/@length - $smotif/motif/@ltrim - $smotif/motif/@rtrim"/><xsl:text> motif-length </xsl:text>
    <xsl:value-of select="/spamo/model/bin_size"/><xsl:text> bin-size </xsl:text>
    <xsl:value-of select="/spamo/model/bin_max" /><xsl:text> bin-max </xsl:text>
    <xsl:value-of select="histogram/@total" /><xsl:text> sequences </xsl:text>
    <xsl:value-of select="/spamo/model/bin_pvalue_cutoff" /><xsl:text> threshold\n\n&quot; + &newline;</xsl:text>
    <xsl:value-of select="$indent"/>
    <xsl:text>&quot;% Same Left\n&quot; + &newline;</xsl:text>
    <xsl:call-template name="spamo-graph-bins">
      <xsl:with-param name="bins" select="histogram/same_strand/left_side/bin"/>
      <xsl:with-param name="indent" select="$indent"/>
    </xsl:call-template>
    <xsl:text> + &newline;</xsl:text>
    <xsl:value-of select="$indent"/>
    <xsl:text>&quot;% Same Right\n&quot; + &newline;</xsl:text>
    <xsl:call-template name="spamo-graph-bins">
      <xsl:with-param name="bins" select="histogram/same_strand/right_side/bin"/>
      <xsl:with-param name="indent" select="$indent"/>
    </xsl:call-template>
    <xsl:text> + &newline;</xsl:text>
    <xsl:value-of select="$indent"/>
    <xsl:text>&quot;% Opposite Left\n&quot; + &newline;</xsl:text>
    <xsl:call-template name="spamo-graph-bins">
      <xsl:with-param name="bins" select="histogram/opposite_strand/left_side/bin"/>
      <xsl:with-param name="indent" select="$indent"/>
    </xsl:call-template>
    <xsl:text> + &newline;</xsl:text>
    <xsl:value-of select="$indent"/>
    <xsl:text>&quot;% Opposite Right\n&quot; + &newline;</xsl:text>
    <xsl:call-template name="spamo-graph-bins">
      <xsl:with-param name="bins" select="histogram/opposite_strand/right_side/bin"/>
      <xsl:with-param name="indent" select="$indent"/>
    </xsl:call-template>

  </xsl:template>

  <xsl:template name="motif-consensus">
    <xsl:param name="motif"/>
    <xsl:param name="id" select="''"/>
    <xsl:param name="col" select="0"/>
    <xsl:param name="rc" select="'false'"/>
    <xsl:param name="box" select="'true'"/>
    <xsl:variable name="seqclass">
      <xsl:text>seq</xsl:text>
      <xsl:if test="$box = 'true'">
        <xsl:text> seqbox</xsl:text>
      </xsl:if>
    </xsl:variable>
    <div class="{$seqclass}" id="{$id}">
      <xsl:call-template name="motif-consensus-pad"><xsl:with-param name="pad" select="$col"/></xsl:call-template>
      <xsl:choose>
        <xsl:when test="$rc = 'true'">
          <xsl:for-each select="$motif/pos">
            <xsl:sort select="@i" data-type="number" order="descending"/>
            <xsl:call-template name="pos-consensus-complement"/>
          </xsl:for-each>
        </xsl:when>
        <xsl:otherwise>
          <xsl:for-each select="$motif/pos"><xsl:call-template name="pos-consensus"/></xsl:for-each>
        </xsl:otherwise>
      </xsl:choose>
    </div>
  </xsl:template>

  <xsl:template name="motif-consensus-pad">
    <xsl:param name="pad" select="0"/>
    <xsl:if test="$pad &gt; 0">
      <xsl:text>&nbsp;</xsl:text>
      <xsl:call-template name="motif-consensus-pad"><xsl:with-param name="pad" select="$pad - 1"/></xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="pos-consensus">
    <xsl:variable name="trim">
      <xsl:if test="../@ltrim &gt;= @i or @i &gt; ../@length - ../@rtrim ">
        <xsl:text>trim</xsl:text>
      </xsl:if>
    </xsl:variable>
    <xsl:choose>
      <xsl:when test="@A &gt;= @C and @A &gt;= @G and @A &gt;= @T">
        <span class="A {$trim}">A</span>
      </xsl:when>
      <xsl:when test="@C &gt;= @A and @C &gt;= @G and @C &gt;= @T">
        <span class="C {$trim}">C</span>
      </xsl:when>
      <xsl:when test="@G &gt;= @A and @G &gt;= @C and @G &gt;= @T">
        <span class="G {$trim}">G</span>
      </xsl:when>
      <xsl:otherwise>
        <span class="T {$trim}">T</span>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="pos-consensus-complement">
    <xsl:variable name="trim">
      <xsl:if test="../@ltrim &gt;= @i or @i &gt; ../@length - ../@rtrim ">
        <xsl:text>trim</xsl:text>
      </xsl:if>
    </xsl:variable>
    <xsl:choose>
      <xsl:when test="@A &gt;= @C and @A &gt;= @G and @A &gt;= @T">
        <span class="T {$trim}">T</span>
      </xsl:when>
      <xsl:when test="@C &gt;= @A and @C &gt;= @G and @C &gt;= @T">
        <span class="G {$trim}">G</span>
      </xsl:when>
      <xsl:when test="@G &gt;= @A and @G &gt;= @C and @G &gt;= @T">
        <span class="C {$trim}">C</span>
      </xsl:when>
      <xsl:otherwise>
        <span class="A {$trim}">A</span>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="header">
    <xsl:param name="title" />
    <xsl:param name="self" select="$title" />
    <xsl:param name="prev" select="''" />
    <xsl:param name="next" select="''" />

    <a name="{$self}"/>
    <table width="100%" border="0" cellspacing="1" cellpadding="4" bgcolor="#FFFFFF">
      <tr>
        <td>
          <h2 class="mainh"><xsl:value-of select="$title"/></h2>
        </td>
        <td align="right" valign="bottom">
          <xsl:if test="$prev != ''"><a href="#{$prev}">Previous</a>&nbsp;</xsl:if>
          <xsl:if test="$next != ''"><a href="#{$next}">Next</a>&nbsp;</xsl:if>
          <a href="#top">Top</a>
        </td>
      </tr>
    </table>
  </xsl:template>

</xsl:stylesheet>
