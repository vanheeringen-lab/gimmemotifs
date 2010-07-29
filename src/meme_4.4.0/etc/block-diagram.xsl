<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="pos-integer-power">
    <xsl:param name="exponent"/>
    <xsl:param name="base"/>
    <xsl:choose>
      <xsl:when test="$exponent &gt; 1">
        <xsl:choose>
          <xsl:when test="($exponent mod 2) = 0">
            <!-- divide the exponent by 2 and recurse then square result -->
            <xsl:variable name="sqrt-result">
              <xsl:call-template name="pos-integer-power">
                <xsl:with-param name="exponent" select="$exponent div 2"/>
                <xsl:with-param name="base" select="$base"/>
              </xsl:call-template>
            </xsl:variable>
            <xsl:value-of select="$sqrt-result * $sqrt-result"/>
          </xsl:when>
          <xsl:otherwise>
            <!-- divide the exponent by 2 and recurse then square result and multiply by base -->
            <xsl:variable name="sqrt-result">
              <xsl:call-template name="pos-integer-power">
                <xsl:with-param name="exponent" select="($exponent - 1) div 2"/>
                <xsl:with-param name="base" select="$base"/>
              </xsl:call-template>
            </xsl:variable>
            <xsl:value-of select="$sqrt-result * $sqrt-result * $base"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$exponent = 1">
        <xsl:value-of select="$base" />
      </xsl:when>
      <xsl:otherwise><!-- assume exponent is zero because we can't handle any other value -->
        <xsl:value-of select="1" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="ln-approx-half">
    <xsl:param name="steps" />
    <xsl:param name="sum" />
    <xsl:param name="factor1" />
    <xsl:param name="factor2" />
    <xsl:param name="factor2multiplier" />
    <xsl:choose>
      <xsl:when test="$steps &gt; 0">
        <xsl:call-template name="ln-approx-half">
          <xsl:with-param name="steps" select="$steps - 1" />
          <xsl:with-param name="sum" select="$sum + ((1 div $factor1) * $factor2)" />
          <xsl:with-param name="factor1" select="$factor1 + 2" />
          <xsl:with-param name="factor2" select="$factor2 * $factor2multiplier" />
          <xsl:with-param name="factor2multiplier" select="$factor2multiplier" />
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$sum" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <!-- this is unfortunately too slow for very small numbers -->
  <xsl:template name="ln-approx2">
    <xsl:param name="base"/>
    <xsl:param name="iterations" select="1000"/>
    <xsl:variable name="factor2" select="($base - 1) div ($base + 1)"/>
    <xsl:variable name="factor2multiplier" select="$factor2 * $factor2"/>
    <xsl:variable name="halfln">
      <xsl:call-template name="ln-approx-half">
        <xsl:with-param name="steps" select="$iterations"/>
        <xsl:with-param name="sum" select="0"/>
        <xsl:with-param name="factor1" select="1"/>
        <xsl:with-param name="factor2" select="$factor2"/>
        <xsl:with-param name="factor2multiplier" select="$factor2multiplier"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:value-of select="2 * $halfln"/>
  </xsl:template>

  <!-- find the significand and exponent -->
  <xsl:template name="ln-approx-sig-exp">
    <xsl:param name="num" />
    <xsl:param name="exp" select="0" />
    <xsl:choose>
      <xsl:when test="$num = 0">
        <!-- shouldn't happen, log(0) is negative infinity -->
        <xsl:value-of select="0"/>
      </xsl:when>
      <xsl:when test="$num &lt; 1">
        <!-- multiply by 2, subtract 1 from exponent -->
        <xsl:call-template name="ln-approx-sig-exp">
          <xsl:with-param name="num" select="$num * 2"/>
          <xsl:with-param name="exp" select="$exp - 1"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:when test="$num &gt;= 2">
        <!-- divide by 2, add 1 to exponent -->
        <xsl:call-template name="ln-approx-sig-exp">
          <xsl:with-param name="num" select="$num div 2"/>
          <xsl:with-param name="exp" select="$exp + 1"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
        <!-- return significand and exponent -->
        <xsl:value-of select="$num"/>
        <xsl:text> </xsl:text>
        <xsl:value-of select="$exp"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="ln-approx">
    <xsl:param name="base"/>
    <xsl:param name="iterations" select="100"/>
    <xsl:variable name="parts">
      <xsl:call-template name="ln-approx-sig-exp">
        <xsl:with-param name="num" select="$base"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="sig" select="substring-before($parts, ' ')"/>
    <xsl:variable name="exp" select="substring-after($parts, ' ')"/>
    <xsl:variable name="ln_2" select="0.6931471805599453094172321214581765680755001343602552541206"/>
    <xsl:variable name="ln_sig">
      <xsl:call-template name="ln-approx2">
        <xsl:with-param name="base" select="$sig"/>
        <xsl:with-param name="iterations" select="$iterations"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:value-of select="$ln_sig + $exp * $ln_2"/>
  </xsl:template>

  <xsl:template name="site">
    <xsl:param name="max_seq_len"/>
    <xsl:param name="max_log_pvalue"/>
    <xsl:param name="position"/>
    <xsl:param name="width"/>
    <xsl:param name="index"/>
    <xsl:param name="strand" select="''"/>
    <xsl:param name="pvalue" select="''"/>
    <xsl:param name="frame" select="''"/>
    <xsl:param name="name"/>
    <xsl:variable name="log_pvalue">
      <xsl:choose>
        <xsl:when test="$pvalue != ''">
          <xsl:call-template name="ln-approx">
            <xsl:with-param name="base" select="$pvalue"/>
            <xsl:with-param name="iterations" select="5"/><!-- a very rough approximation of the log based mostly on the exponent -->
          </xsl:call-template>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="$max_log_pvalue" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="max_block_height" select="12" />
    <xsl:variable name="divider_height" select="1" />
    <xsl:variable name="motif_fract" select="round(($width div $max_seq_len) * 100000) div 1000" />
    <xsl:variable name="offset_fract" select="round((($position - 1) div $max_seq_len) * 100000) div 1000" />
    <xsl:variable name="motif_impact" >
      <!-- assume both logs are negative (as they should be for pvalues)  -->
      <xsl:choose>
        <xsl:when test="$log_pvalue &gt; $max_log_pvalue">
          <xsl:value-of select="round($max_block_height * ($log_pvalue div $max_log_pvalue) * 10) div 10" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="$max_block_height" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="motif_color">
      <xsl:call-template name="pick_color">
        <xsl:with-param name="index" select="$index"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="motif_border">
      <xsl:call-template name="pick_border">
        <xsl:with-param name="index" select="$index"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="offset_top">
      <xsl:choose>
        <xsl:when test="$strand='plus'">
          <xsl:value-of select="$max_block_height - $motif_impact - 1"/>
        </xsl:when>
        <xsl:when test="$strand='minus'">
          <xsl:value-of select="$max_block_height"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="$max_block_height - $motif_impact - 1"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <xsl:variable name="title">
      <xsl:text>Motif </xsl:text><xsl:value-of select="$name"/><xsl:text>    </xsl:text>
      <xsl:if test="$pvalue != ''">
        <xsl:text>p-value: </xsl:text><xsl:value-of select="$pvalue"/><xsl:text>    </xsl:text>
      </xsl:if>
      <xsl:if test="$frame != ''">
        <xsl:text>frame: </xsl:text><xsl:value-of select="$frame"/><xsl:text>    </xsl:text>
      </xsl:if>
      <xsl:text>starts: </xsl:text><xsl:value-of select="$position"/><xsl:text>    </xsl:text>
      <xsl:text>ends: </xsl:text><xsl:value-of select="$position + $width - 1"/><xsl:text>    </xsl:text>
    </xsl:variable>
    <div class="block_motif" 
      style="left:{$offset_fract}%; top:{$offset_top}px; width:{$motif_fract}%; height:{$motif_impact}px; background-color:{$motif_color}; border: {$motif_border};"
      title="{$title}">
    </div>
  </xsl:template>

  <xsl:template name="ruler">
    <xsl:param name="max" />              <!-- the maximum value on the ruler -->
    <xsl:call-template name="ruler-tic">
      <xsl:with-param name="max" select="$max"/>
      <xsl:with-param name="step">
        <xsl:choose>
          <xsl:when test="$max &lt; 50">
            <xsl:value-of select="'1'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 100">
            <xsl:value-of select="'2'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 200">
            <xsl:value-of select="'4'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 500">
            <xsl:value-of select="'10'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 1000">
            <xsl:value-of select="'20'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 2000">
            <xsl:value-of select="'40'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 5000">
            <xsl:value-of select="'100'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 10000">
            <xsl:value-of select="'200'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 20000">
            <xsl:value-of select="'400'"/>
          </xsl:when>
          <xsl:otherwise> 
            <xsl:value-of select="floor($max div 20000) * 400" />
          </xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
      <xsl:with-param name="period">
        <xsl:choose>
          <xsl:when test="$max &lt; 10">
            <xsl:value-of select="'1'"/>
          </xsl:when>
          <xsl:when test="$max &lt; 20">
            <xsl:value-of select="'2'"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="'5'"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:with-param>
    </xsl:call-template>
  </xsl:template>

  <xsl:template name="ruler-tic">
    <xsl:param name="pos" select="0" />   <!-- the position to put the next tic -->
    <xsl:param name="step" select="20" /> <!-- the step to increment between tics -->
    <xsl:param name="max" />              <!-- the maximum value on the ruler -->
    <xsl:param name="period" select="5" /><!-- the number of minor tics between major tics -->
    <xsl:param name="cycle" select="0" /> <!-- the current position in the cycle of major/minor tics this is always a mod of peroid -->

    <xsl:variable name="fract" select="round(($pos div $max) * 100000) div 1000" />
    <xsl:variable name="label-len" select="string-length($pos)" />
    <xsl:variable name="label-offset" select="$label-len div 2" />

    <xsl:choose>
      <xsl:when test="$cycle = 0">
        <div class="tic_major" style="left:{$fract}%">
          <div class="tic_label" style="left:-{$label-offset}em; width:{$label-len}em;">
            <xsl:value-of select="$pos"/>
          </div>
        </div>
      </xsl:when>
      <xsl:otherwise>
        <div class="tic_minor" style="left:{$fract}%;"/>
      </xsl:otherwise>
    </xsl:choose>

    <xsl:if test="($pos + $step) &lt;= $max">
      <xsl:call-template name="ruler-tic">
        <xsl:with-param name="pos" select="$pos + $step" />   
        <xsl:with-param name="step" select="$step" /> 
        <xsl:with-param name="period" select="$period" />
        <xsl:with-param name="cycle" select="($cycle + 1) mod $period" /> 
        <xsl:with-param name="max" select="$max" />
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="legend_motif">
    <xsl:param name="name" />
    <xsl:param name="index" />
    <xsl:variable name="motif_colour">
      <xsl:call-template name="pick_color">
        <xsl:with-param name="index" select="$index"/>
      </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="motif_border">
      <xsl:call-template name="pick_border">
        <xsl:with-param name="index" select="$index"/>
      </xsl:call-template>
    </xsl:variable>
    <table style="display:inline-block; padding:5px;">
    <tr>
      <td><div style="width: 20px; height: 20px; background-color:{$motif_colour}; border: {$motif_border}"></div></td>
      <td>Motif <xsl:value-of select="$name" /></td>
    </tr>
    </table>
  </xsl:template>


  <xsl:template name="pick_color">
    <xsl:param name="index" select="1"/>
    <xsl:variable name="opt" select="($index - 1) mod 15"/>
    <xsl:choose>
      <xsl:when test="$opt = 0">aqua</xsl:when>
      <xsl:when test="$opt = 1">blue</xsl:when>
      <xsl:when test="$opt = 2">red</xsl:when>
      <xsl:when test="$opt = 3">fuchsia</xsl:when>
      <xsl:when test="$opt = 4">yellow</xsl:when>
      <xsl:when test="$opt = 5">lime</xsl:when>
      <xsl:when test="$opt = 6">teal</xsl:when>
      <xsl:when test="$opt = 7">#444444</xsl:when>
      <xsl:when test="$opt = 8">green</xsl:when>
      <xsl:when test="$opt = 9">silver</xsl:when>
      <xsl:when test="$opt = 10">purple</xsl:when>
      <xsl:when test="$opt = 11">olive</xsl:when>
      <xsl:when test="$opt = 12">navy</xsl:when>
      <xsl:when test="$opt = 13">maroon</xsl:when>
      <xsl:when test="$opt = 14">white</xsl:when>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="pick_border">
    <xsl:param name="index" select="1"/>
    <xsl:variable name="opt" select="floor(($index - 1) div 15) mod 3"/>
    <xsl:choose>
      <xsl:when test="$opt = 0">1px solid black</xsl:when>
      <xsl:when test="$opt = 1">1px dashed grey</xsl:when>
      <xsl:when test="$opt = 2">1px dotted grey</xsl:when>
    </xsl:choose>
  </xsl:template>
  
</xsl:stylesheet>

