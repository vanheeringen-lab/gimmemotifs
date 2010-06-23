<?python
title = "Gimme motifs report"
?>
<html xmlns:py="http://purl.org/kid/ns#">
  <head>
      <title py:content="title">This is replaced.</title>
	    </head>
		  <body>
<table>
<tr>
<th>
Motif
</th>
<th colspan="2" py:if="random">
Sig vs random
</th>
<th colspan="2" py:if="genomic">
Sig vs genomic
</th>
<th>
Location
</th>
<th>
Closest match in JASPAR
</th>


</tr>
<tr py:for="motif in motifs">

<td>
<table>
<tr>
<td>
	<img height="80" py:attrs="motif.img"/>
</td>
</tr>
<tr>
<td align="center" py:content="motif.consensus">
rrrCwwGyyy
</td>
</tr>
</table>
</td>


<td py:if="random">
<table>
<tr><td>enrichment</td><td py:content="motif.random_e">100.8</td></tr>
<tr><td>p-value</td><td py:content="motif.random_p">1e-10</td></tr>
<tr><td>ROC_AUC</td><td py:content="motif.random_auc">0.91</td></tr>
</table>
</td>
<td py:if="random">
	<a href="roc_big.png"><img height="100" py:attrs="motif.random_roc_img"/></a>
</td>

<td py:if="genomic">
<table>
<tr><td>enrichment</td><td py:content="motif.genomic_e">100.8</td></tr>
<tr><td>p-value</td><td py:content="motif.genomic_p">1e-10</td></tr>
<tr><td>ROC_AUC</td><td py:content="motif.genomic_auc">0.91</td></tr>
</table>
</td>
<td py:if="genomic">
	<a href="roc_big.png"><img height="100" py:attrs="motif.genomic_roc_img"/></a>
</td>

<td>
<a py:attrs="motif.histogram_link"><img py:attrs="motif.histogram_img" height="100" /></a>
</td>

<td>
<table>
<tr>
<td>
	<img height="40" py:attrs="motif.match_img"/>
</td>
</tr>
<tr>
<td align="center" py:content="motif.match_id">
</td>
</tr>
<tr>
<td align="center" py:content="motif.match_pval">
</td>
</tr>
</table>
</td>

</tr>
</table>
</body>
</html>
