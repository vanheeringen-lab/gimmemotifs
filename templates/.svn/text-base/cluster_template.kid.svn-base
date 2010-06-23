<?python
title = "Gimme motifs cluster report"
?>
<html xmlns:py="http://purl.org/kid/ns#">
  <head>
      <title py:content="title">This is replaced.</title>
	    </head>
		  <body>
<table border="1">
<tr py:for="motif in motifs">
<td py:content="motif[0]"/>
<td>
<img height="80" py:attrs="motif[1]"/>
</td>
<td>
<table>
<tr py:for="m in motif[2]">
<td>
<img height="25" py:attrs="m"/>
</td>
<td py:content='m["alt"]'/>
</tr>
</table>
</td>
</tr>
</table>
</body>
</html>
