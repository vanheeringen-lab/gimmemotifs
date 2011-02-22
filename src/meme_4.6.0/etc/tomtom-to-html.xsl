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
  <xsl:include href="meme.css.xsl"/>
  <xsl:include href="motif_logo.js.xsl"/>

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
      <title>TOMTOM</title>
      <xsl:call-template name="html-css"/>
      <xsl:call-template name="html-script"/>
    </head>
  </xsl:template>

  <xsl:template name="html-css">
    <style type="text/css">
      <xsl:call-template name="meme.css"/>
      span.A {font-size:20px; color:red}
      span.C {font-size:20px; color:blue}
      span.G {font-size:20px; color:orange}
      span.T {font-size:20px; color:green}
      .pt {padding: 0 10px}
      .ac {text-align: center;}
      .downloadTd {padding-left:20px;}
      div.logo_container {position:relative; width:99%; height:285px; padding:0px; margin:0px;}
      img.logo {position:absolute; z-index:2; max-width:100%;}
      tr.tspace th, tr.tspace td {padding-top: 20px;}
    </style>
  </xsl:template>

  <xsl:template name="html-script">
    <script type="text/javascript">
      var show_opts_link;
      var task_queue = new Array();
      var task_delay = 100;
      var my_alphabet;
      var query_pspm;
      <xsl:call-template name="motif_logo.js"/>
      function LoadQueryTask(target_id, query_name) {
        this.target_id = target_id;
        this.query_name = query_name;
        this.run = LoadQueryTask_run;
      }
      function LoadQueryTask_run() {
        var alpha = new Alphabet("ACGT");
        var query_pspm = new Pspm(document.getElementById('q_' + this.query_name).value, this.query_name);
        replace_logo(logo_1(alpha, "Tomtom", query_pspm), this.target_id, 0.5, "Preview of " + this.query_name, "block");
      }
      function FixLogoTask(image_id, query_name, target_db_id, target_name, target_rc, target_offset) {
        this.image_id = image_id;
        this.query_name = query_name;
        this.target_db_id = target_db_id;
        this.target_name = target_name;
        this.target_rc = target_rc;
        this.target_offset = target_offset;
        this.run = FixLogoTask_run;
      }
      function FixLogoTask_run() {
        var image = document.getElementById(this.image_id);
        var canvas = create_canvas(image.width, image.height, image.id, image.title, image.style.display);
        if (canvas == null) return;
        var alpha = new Alphabet("ACGT");
        var query_pspm = new Pspm(document.getElementById('q_' + this.query_name).value, this.query_name);
        var target_pspm = new Pspm(document.getElementById('t_' + this.target_db_id + '_' + this.target_name).value, this.target_name);
        if (this.target_rc) target_pspm.reverse_complement(alpha);
        var logo = logo_2(alpha, "Tomtom", target_pspm, query_pspm, this.target_offset);
        draw_logo_on_canvas(logo, canvas, true, 1);
        image.parentNode.replaceChild(canvas, image);
      }
      function push_task(task) {
        task_queue.push(task);
        if (task_queue.length == 1) {
          window.setTimeout("process_tasks()", task_delay);
        }
      }
      function process_tasks() {
        if (task_queue.length == 0) return; //no more tasks
        //get next task
        var task = task_queue.shift();
        task.run();
        //allow UI updates between tasks
        window.setTimeout("process_tasks()", task_delay);
      }
      function showHidden(prefix) {
        document.getElementById(prefix + '_activator').style.display = 'none';
        document.getElementById(prefix + '_deactivator').style.display = 'block';
        document.getElementById(prefix + '_data').style.display = 'block';
      }
      function hideShown(prefix) {
        document.getElementById(prefix + '_activator').style.display = 'block';
        document.getElementById(prefix + '_deactivator').style.display = 'none';
        document.getElementById(prefix + '_data').style.display = 'none';
      }
      function show_more(query_index, match_index, query_id, target_id, offset, rc) {
        var link = document.getElementById("show_" + query_index + "_" + match_index);
        var match_table = document.getElementById("table_" + query_index);
        var match_row = document.getElementById("tr_" + query_index + "_" + match_index);
        var query_pspm = document.getElementById(query_id).value;
        var query_length = (new Pspm(query_pspm)).get_motif_length();
        var target_pspm = document.getElementById(target_id).value;
        var target_length = (new Pspm(target_pspm)).get_motif_length();
        var total_length = (offset &lt; 0 ? Math.max(query_length, target_length - offset) : Math.max(query_length + offset, target_length));
        var download_row = document.getElementById("tr_download");
        if (download_row) {
          //remove existing
          match_table.deleteRow(download_row.rowIndex);
          show_opts_link.style.visibility = "visible";
        }
        show_opts_link = link;
        show_opts_link.style.visibility = "hidden";
        download_row = match_table.insertRow(match_row.rowIndex + 1);
        download_row.id = "tr_download";
        var cell = download_row.insertCell(0);
        cell.colSpan = 2;
        cell.style.verticalAlign = "top";
        //create the download options
        var box = document.createElement("div");
        box.appendChild(create_title("Create custom LOGO", "#download_logo_doc"));
        var form = create_form("<xsl:value-of select="concat($site_url, '/cgi-bin/tomtom_request.cgi')"/>", "post");
        add_hfield(form, "mode", "logo");
        add_hfield(form, "version", "<xsl:value-of select="/tomtom/@version"/>");
        add_hfield(form, "background", "<xsl:apply-templates select="/tomtom/model/background"/>");
        add_hfield(form, "target_id", target_id);
        add_hfield(form, "target_length", target_length);
        add_hfield(form, "target_pspm", target_pspm);
        add_hfield(form, "target_rc", rc);
        add_hfield(form, "query_id", query_id);
        add_hfield(form, "query_length", query_length);
        add_hfield(form, "query_pspm", query_pspm);
        add_hfield(form, "query_offset", offset);
        var tbl = document.createElement("table");
        var row = tbl.insertRow(tbl.rows.length);
        add_label(row.insertCell(-1),"Image Type");
        add_label(row.insertCell(-1), "Error bars");
        add_label(row.insertCell(-1), "SSC");
        add_label(row.insertCell(-1), "Flip");
        add_label(row.insertCell(-1), "Width");
        add_label(row.insertCell(-1), "Height");
        row.insertCell(-1).appendChild(document.createTextNode("&nbsp;"));
        row = tbl.insertRow(tbl.rows.length);
        add_combo(row.insertCell(-1), "image_type", {"png" : "PNG (for web)", "eps" : "EPS (for publication)"}, "png");
        add_combo(row.insertCell(-1), "error_bars", {"1" : "yes", "0" : "no"}, "1");
        add_combo(row.insertCell(-1), "small_sample_correction", {"0" : "no", "1" : "yes"}, "0");
        add_combo(row.insertCell(-1), "flip", {"0" : "no", "1" : "yes"}, "0");
        add_tfield(row.insertCell(-1), "image_width", total_length, 2);
        add_tfield(row.insertCell(-1), "image_height", 10, 2);
        add_submit(row.insertCell(-1), "Download");
        form.appendChild(tbl);
        box.appendChild(form);
        cell.appendChild(box);
      }
      function create_title(title, help_link) {
        var title_ele = document.createElement("h4");
        var words = document.createTextNode(title + "&nbsp;");
        title_ele.appendChild(words);
        if (help_link) {
          var link = document.createElement("a");
          link.className = "help"
          link.href = help_link;
          var helpdiv = document.createElement("div");
          helpdiv.className = "help";
          link.appendChild(helpdiv);
          title_ele.appendChild(link);
        }
        return title_ele;
      }
      function add_hfield(form, name, value) {
        var hidden = document.createElement("input");
        hidden.type = "hidden";
        hidden.name = name;
        hidden.value = value;
        form.appendChild(hidden);
      }
      function add_label(ele,label) {
        ele.className = "downloadTd";
        var label_b = document.createElement("b");
        label_b.appendChild(document.createTextNode(label));
        ele.appendChild(label_b);
      }
      function add_tfield(ele, name, value, size) {
        ele.className = "downloadTd";
        var field = document.createElement("input");
        field.type = "text";
        field.style.marginLeft = "5px";
        field.name = name;
        field.value = value;
        field.size = size;
        ele.appendChild(field);
      }
      function add_combo(ele, name, options, selected) {
        ele.className = "downloadTd";
        var combo = document.createElement("select");
        combo.name = name;
        for (var value in options) {
          combo.style.marginLeft = "5px";
          //exclude inherited properties and undefined properties
          if (!options.hasOwnProperty(value) || options[value] === undefined) continue;
          var display = options[value];
          var option = document.createElement("option");
          option.value = value;
          if (value == selected) option.selected = true;
          option.appendChild(document.createTextNode(display));
          combo.appendChild(option);
        }
        ele.appendChild(combo);
      }
      function add_submit(ele, value) {
        ele.className = "downloadTd";
        var submit = document.createElement("input");
        submit.type = "submit";
        submit.value = value;
        ele.appendChild(submit);
      }
      function create_form(path,method) {
        method = method || "post";

        var form = document.createElement("form");
        form.method = method;
        form.action = path;
        return form;
      }
    </script>
  </xsl:template>

  <xsl:template match="background">
    <xsl:variable name="psum" select="@A + @C + @G + @T"/>
    <xsl:variable name="delta" select="0.01"/>
    <xsl:if test="($psum &gt; 1 and $psum - 1 &gt; $delta) or ($psum &lt; 1 and 1 - $psum &gt; $delta)">
      <xsl:message>Warning: background probabilities don't sum to 1 (delta <xsl:value-of select="$delta"/>).</xsl:message>
    </xsl:if>
    <xsl:value-of select="@A"/><xsl:text> </xsl:text><xsl:value-of select="@C"/>
    <xsl:text> </xsl:text><xsl:value-of select="@G"/><xsl:text> </xsl:text><xsl:value-of select="@T"/>
  </xsl:template>

  <xsl:template name="html-body">
    <body>
      <xsl:call-template name="data"/>
      <xsl:call-template name="top"/>
      <xsl:call-template name="preview"/>
      <xsl:call-template name="targets"/>
      <xsl:call-template name="queries"/>
      <xsl:call-template name="program"/>
      <xsl:call-template name="documentation"/>
      <xsl:call-template name="footer"/>
    </body>
  </xsl:template>

  <xsl:template name="data">
    <form>
      <xsl:for-each select="/tomtom/queries/query_file/query">
        <xsl:call-template name="motif-pspm">
          <xsl:with-param name="motif" select="motif"/>
        </xsl:call-template>
      </xsl:for-each>

      <xsl:for-each select="/tomtom/targets/target_file/motif">
        <xsl:call-template name="motif-pspm">
          <xsl:with-param name="motif" select="."/>
        </xsl:call-template>
      </xsl:for-each>
    </form>
  </xsl:template>

  <xsl:template name="motif-pspm">
    <xsl:param name="motif"/>
    <xsl:variable name="pspm">
      <xsl:text>&newline;</xsl:text>
      <xsl:text>letter-probability matrix: alength= 4 w= </xsl:text><xsl:value-of select="$motif/@length"/>
      <xsl:text> nsites= </xsl:text><xsl:value-of select="$motif/@nsites"/><xsl:text> E= </xsl:text>
      <xsl:value-of select="$motif/@evalue"/><xsl:text> &newline;</xsl:text>
      <xsl:for-each select="$motif/pos">
        <xsl:variable name="psum" select="@A + @C + @G + @T"/>
        <xsl:variable name="delta" select="0.01"/>
        <xsl:if test="($psum &gt; 1 and $psum - 1 &gt; $delta) or ($psum &lt; 1 and 1 - $psum &gt; $delta)">
          <xsl:message>
            Warning: Motif <xsl:value-of select="@id"/> probabilities at pos <xsl:value-of select="@i"/> don't sum to 1 
            (delta <xsl:value-of select="$delta"/>).
          </xsl:message>
        </xsl:if>
        <xsl:value-of select="@A"/><xsl:text> </xsl:text><xsl:value-of select="@C"/><xsl:text> </xsl:text>
        <xsl:value-of select="@G"/><xsl:text> </xsl:text><xsl:value-of select="@T"/>
        <xsl:text>&newline;</xsl:text>
      </xsl:for-each>
    </xsl:variable>
    <input type="hidden" id="{$motif/@id}" value="{$pspm}"/><xsl:text>&newline;</xsl:text>
  </xsl:template>

  <xsl:template name="top">
    <a name="top"/>
    <div class="pad1"><!-- TODO find out true acronym for TOMTOM -->
      <h1><img src="http://meme.nbcr.net/meme/doc/images/tomtom_logo.png" alt="TOol for Motif TO Motif comparison (TOMTOM)" /></h1>
      <p class="spaced">
        For further information on how to interpret these results or to get a 
        copy of the MEME software please access 
        <a href="http://meme.nbcr.net/">http://meme.nbcr.net</a>. 
      </p>
    </div>
    <!-- navigation -->
    <div class="pad2">
      <a class="jump" href="#query_motifs">Query Motifs</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#target_dbs">Target Databases</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#query_{/tomtom/queries/query_file[1]/query[1]/motif/@id}">Matches</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#program">Program information</a>
      <xsl:text>&nbsp;&nbsp;|&nbsp;&nbsp;</xsl:text>
      <a class="jump" href="#doc">Explanation</a>
    </div>
  </xsl:template>

  <xsl:template name="preview">
    <xsl:variable name="nmatch" select="20"/>
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Query Motifs'"/>
      <xsl:with-param name="self" select="'query_motifs'"/>
      <xsl:with-param name="next" select="'target_dbs'"/>
    </xsl:call-template>
    <div class="box">
      <table>
        <thead>
          <tr>
            <td class="pt" style="min-width:5em;"><h4>Name&nbsp;<a href="#qm_name_doc" class="help"><div class="help"/></a></h4></td>
            <xsl:if test="/tomtom/queries/query_file/query/motif/@alt">
              <td class="pt" style="min-width:8em;"><h4>Alt. Name&nbsp;<a href="#qm_alt_doc" class="help"><div class="help"/></a></h4></td>
            </xsl:if>
            <xsl:if test="/tomtom/queries/query_file/query/motif/@url">
              <td class="pt" style="min-width:8em;"><h4>Website&nbsp;<a href="#qm_website_doc" class="help"><div class="help"/></a></h4></td>
            </xsl:if>
            <td class="pt" style="min-width:9em;"><h4>Preview&nbsp;<a href="#qm_preview_doc" class="help"><div class="help"/></a></h4></td>
            <td class="pt" style="min-width:9em;"><h4>Matches&nbsp;<a href="#qm_matches_doc" class="help"><div class="help"/></a></h4></td>
            <td class="pt" style="min-width:6em;"><h4>List&nbsp;<a href="#qm_list_doc" class="help"><div class="help"/></a></h4></td>
          </tr>
        </thead>
        <tbody>
          <xsl:for-each select="/tomtom/queries/query_file/query">
            <xsl:variable name="query" select="motif"/>
            <tr>
              <td class="pt"><a href="#query_{motif/@id}"><xsl:value-of select="motif/@name"/></a></td>
              <xsl:if test="/tomtom/queries/query_file/query/motif/@alt">
                <td class="pt"><a href="#query_{motif/@id}"><xsl:value-of select="motif/@alt"/></a></td>
              </xsl:if>
              <xsl:if test="motif/@url">
                <td class="pt"><a href="{motif/@url}">link</a></td>
              </xsl:if>
              <td class="pt"><a href="#query_{motif/@id}" style="color: inherit; text-decoration:none;"><span><div id="preview_{position()}">
                    <script type="text/javascript">
                      <xsl:text>push_task(new LoadQueryTask("preview_" + </xsl:text>
                      <xsl:value-of select="position()"/>
                      <xsl:text>, &quot;</xsl:text>
                      <xsl:value-of select="motif/@name"/>
                      <xsl:text>&quot;));</xsl:text>
                    </script>
                  <xsl:for-each select="motif/pos"><xsl:call-template name="pos-consensus"/></xsl:for-each>
              </div></span></a></td>
              <td class="pt ac"><a href="#query_{motif/@id}"><xsl:value-of select="count(match)"/></a></td>
              <td class="pt">
                <xsl:for-each select="match">
                  <xsl:if test="position() -1 &lt; $nmatch">
                    <xsl:variable name="target" select="id(@target)"/>
                    <a href="#match_{$query/@id}_{$target/@id}">
                      <xsl:choose>
                        <xsl:when test="$target/@name">
                          <xsl:value-of select="$target/@name"/>
                        </xsl:when>
                        <xsl:otherwise>
                          <xsl:value-of select="$target/@id"/>
                        </xsl:otherwise>
                      </xsl:choose>
                      <xsl:if test="$target/@alt">
                        (<xsl:value-of select="$target/@alt"/>)
                      </xsl:if>
                    </a>
                    <xsl:if test="position() != last() and position() != $nmatch">
                      <xsl:text>, </xsl:text>
                    </xsl:if>
                  </xsl:if>
                </xsl:for-each>
              </td>
            </tr>
          </xsl:for-each>
        </tbody>
      </table>
    </div>
  </xsl:template>
  
  <xsl:template name="pos-consensus">
    <xsl:choose>
      <xsl:when test="@A &gt;= @C and @A &gt;= @G and @A &gt;= @T">
        <span class="A">A</span>
      </xsl:when>
      <xsl:when test="@C &gt;= @A and @C &gt;= @G and @C &gt;= @T">
        <span class="C">C</span>
      </xsl:when>
      <xsl:when test="@G &gt;= @A and @G &gt;= @C and @G &gt;= @T">
        <span class="G">G</span>
      </xsl:when>
      <xsl:otherwise>
        <span class="T">T</span>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="targets">
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Target Databases'"/>
      <xsl:with-param name="self" select="'target_dbs'"/>
      <xsl:with-param name="prev" select="'query_motifs'"/>
      <xsl:with-param name="next" select="concat('query_',/tomtom/queries/query_file[1]/query[1]/motif/@id)"/>
    </xsl:call-template>
    <div class="box">
      <table>
        <thead>
          <tr>
            <td class="pt" style="min-width:8em;"><h4>Database&nbsp;<a href="#tdb_database_doc" class="help"><div class="help"/></a></h4></td>
            <td class="pt" style="min-width:8em;"><h4>Number of Motifs&nbsp;<a href="#tdb_contains_doc" class="help"><div class="help"/></a></h4></td>
            <td class="pt" style="min-width:8em;"><h4>Motifs Matched&nbsp;<a href="#tdb_matched_doc" class="help"><div class="help"/></a></h4></td>
          </tr>
        </thead>
        <tbody>
          <xsl:for-each select="/tomtom/targets/target_file">
            <tr>
              <td class="pt"><xsl:value-of select="@name"/></td>
              <td class="pt ac"><xsl:value-of select="@loaded - @excluded"/></td>
              <td class="pt ac"><xsl:value-of select="count(motif)"/></td>
            </tr>
          </xsl:for-each>
        </tbody>
      </table>
    </div>
  </xsl:template>

  <xsl:template name="queries">
    <xsl:for-each select="/tomtom/queries/query_file/query" >
      <xsl:variable name="query_index" select="position()"/>
      <xsl:call-template name="header">
        <xsl:with-param name="title">
          <xsl:text>Matches to Query: </xsl:text><xsl:value-of select="motif/@name"/>
          <xsl:if test="motif/@alt and motif/@alt != ''">
            <xsl:text> (</xsl:text><xsl:value-of select="motif/@alt"/><xsl:text>)&newline;</xsl:text>
          </xsl:if>
        </xsl:with-param>
        <xsl:with-param name="self" select="concat('query_', motif/@id)"/>
        <xsl:with-param name="prev"><xsl:variable name="sib" select="preceding-sibling::query[1]/motif/@id" />
          <xsl:choose>
            <xsl:when test="$sib"><xsl:value-of select="concat('query_', $sib)"/></xsl:when>
            <xsl:otherwise><xsl:text>target_dbs</xsl:text></xsl:otherwise>
          </xsl:choose>
        </xsl:with-param>
        <xsl:with-param name="next"><xsl:variable name="sib" select="following-sibling::query/motif/@id" />
          <xsl:choose>
            <xsl:when test="$sib"><xsl:value-of select="concat('query_', $sib)"/></xsl:when>
            <xsl:otherwise><xsl:text>program</xsl:text></xsl:otherwise>
          </xsl:choose>
        </xsl:with-param>

      </xsl:call-template>
      <div class="box">
        <table id="table_{$query_index}" style="width:100%">
          <col style="width:20%;"/>
          <col />
        <xsl:for-each select="match">
          <xsl:call-template name="match">
            <xsl:with-param name="query_index" select="$query_index"/>
            <xsl:with-param name="match_index" select="position()"/>
          </xsl:call-template>
        </xsl:for-each>
        </table>
      </div>
    </xsl:for-each>
  </xsl:template>

  <xsl:template name="match">
    <xsl:param name="query_index"/>
    <xsl:param name="match_index"/>
    <xsl:variable name="query" select="../motif"/>
    <xsl:variable name="target" select="id(@target)"/>
    <xsl:variable name="prev_target" select="preceding-sibling::match[1]/@target"/>
    <xsl:variable name="next_target" select="following-sibling::match/@target"/>

    <xsl:variable name="overlap">
      <xsl:call-template name="logo-overlap">
        <xsl:with-param name="query" select="$query"/>
        <xsl:with-param name="target" select="$target"/>
        <xsl:with-param name="offset" select="@offset"/>
      </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="rc">
      <xsl:choose>
        <xsl:when test="@orientation = 'forward'">0</xsl:when>
        <xsl:otherwise>1</xsl:otherwise>
      </xsl:choose>
    </xsl:variable>
    <tr>
      <td style="min-width:9em;">
        <a name="match_{$query/@id}_{$target/@id}"/>
        <h4>Summary&nbsp;<a href="#summary_doc" class="help"><div class="help"/></a></h4>
      </td>
      <td style="min-width:11em;"><h4>Alignment&nbsp;<a href="#alignment_doc" class="help"><div class="help"/></a></h4></td>
    </tr>
    <tr id="tr_{$query_index}_{$match_index}">
      <td style="vertical-align:top">
        <table style="width:100%">
          <colgroup/>
          <colgroup align="left" />
          <tr>
            <th>Name</th>
            <td>
              <xsl:choose>
                <xsl:when test="$target/@url"><a href="{$target/@url}"><xsl:value-of select="$target/@name"/></a></xsl:when>
                <xsl:otherwise><xsl:value-of select="$target/@name"/></xsl:otherwise>
              </xsl:choose>
            </td>
          </tr>
          <xsl:if test="$target/@alt"><tr><th>Alt. Name</th><td><xsl:value-of select="$target/@alt"/></td></tr></xsl:if>
          <tr><th>Database</th><td><xsl:value-of select="$target/../@name"/></td></tr>
          <tr class="tspace"><th><i>p</i>-value</th><td><xsl:value-of select="@pvalue"/></td></tr>
          <tr><th><i>E</i>-value</th><td><xsl:value-of select="@evalue"/></td></tr>
          <tr><th><i>q</i>-value</th><td><xsl:value-of select="@qvalue"/></td></tr>
          <tr class="tspace"><th>Overlap</th><td><xsl:value-of select="$overlap"/></td></tr>
          <tr><th>Offset</th><td><xsl:value-of select="@offset"/></td></tr>
          <tr><th>Orientation</th><td>
              <xsl:choose><xsl:when test="@orientation = 'forward'">Normal</xsl:when><xsl:otherwise>Reverse Complement</xsl:otherwise></xsl:choose>
          </td></tr>
        </table>
      </td>
      <td>
        <xsl:variable name="logo_name"><xsl:call-template name="logo-name"/></xsl:variable>
        <xsl:variable name="target_rc">
          <xsl:choose><xsl:when test="@orientation = 'forward'">false</xsl:when><xsl:otherwise>true</xsl:otherwise></xsl:choose>
        </xsl:variable>
        <div class="logo_container" >
          <img id="m_{$query_index}_{$match_index}" src="{$logo_name}" class="logo" 
            onerror='push_task(new FixLogoTask(this.id, &quot;{$query/@name}&quot;, {$target/../@index}, &quot;{$target/@name}&quot;, {$target_rc}, {@offset}));'/>
        </div>
      </td>
    </tr>
    <tr>
      <td style="text-align:center;">
        <a id="show_{$query_index}_{$match_index}" href="javascript:show_more({$query_index}, {$match_index}, '{../motif/@id}', '{@target}', {@offset}, {$rc})">Create custom LOGO &more;</a>
      </td>
      <td style="text-align:right">
        <xsl:if test="$prev_target">[<a href="#match_{$query/@id}_{$prev_target}">Previous Match</a>]</xsl:if>
        <xsl:if test="$next_target">&nbsp;[<a href="#match_{$query/@id}_{$next_target}">Next Match</a>]</xsl:if>
        <xsl:text>&nbsp;[</xsl:text><a href="#query_{$query/@id}">Query Top</a><xsl:text>]</xsl:text>
      </td>
    </tr>
  </xsl:template>


  <xsl:template name="logo-name">
    <xsl:text>align_</xsl:text><xsl:value-of select="../motif/@name"/><xsl:text>_</xsl:text>
    <xsl:value-of select="id(@target)/../@index"/><xsl:text>_</xsl:text>
    <xsl:choose><xsl:when test="@orientation = 'forward'"><xsl:text>+</xsl:text></xsl:when>
      <xsl:otherwise><xsl:text>-</xsl:text></xsl:otherwise></xsl:choose>
    <xsl:value-of select="id(@target)/@name"/><xsl:text>.png</xsl:text>
  </xsl:template>

  <xsl:template name="logo-overlap">
    <xsl:param name="query"/>
    <xsl:param name="target"/>
    <xsl:param name="offset"/>
    <xsl:choose>
      <xsl:when test="$offset &lt; 0">
        <xsl:variable name="qlen" select="$query/@length + $offset"/>
        <xsl:choose>
          <xsl:when test="$qlen &gt; $target/@length">
            <xsl:value-of select="$target/@length"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$qlen"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$offset = 0">
        <xsl:choose>
          <xsl:when test="$query/@length &gt; $target/@length">
            <xsl:value-of select="$target/@length"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$query/@length"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <xsl:variable name="tlen" select="$target/@length - $offset"/>
        <xsl:choose>
          <xsl:when test="$query/@length &gt; $tlen">
            <xsl:value-of select="$tlen"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="$query/@length"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="program">
    <a name="program"/>
    <div class="bar">
      <div style="text-align:right;"><a href="#query_{/tomtom/queries/query_file[1]/query[last()]/motif/@id}">Previous</a><xsl:text> </xsl:text>
        <a href="#doc">Next</a><xsl:text> </xsl:text><a href="#top">Top</a></div>
      <div class="subsection">
        <a name="version"/>
        <h5>TOMTOM version</h5>
        <xsl:value-of select="/tomtom/@version"/> (Release date: <xsl:value-of select="/tomtom/@release"/>)
      </div>
      <div class="subsection">
        <a name="reference"/>
        <h5>Reference</h5>
        Shobhit Gupta, John A. Stamatoyannopoulos, Timothy L. Bailey and William Stafford Noble, "Quantifying similarity between motifs", <br />
        Genome Biology, 8:R24, 2007
      </div>
      <div class="subsection">
        <a name="command" />
        <h5>Command line summary</h5>
        <textarea rows="1" style="width:100%;" readonly="readonly">
          <xsl:value-of select="/tomtom/model/command_line"/>
        </textarea>
        <br />
        <xsl:text>Background letter frequencies (from </xsl:text>
        <xsl:choose>
          <xsl:when test="/tomtom/model/background/@from = 'first_target'">
            <xsl:text>first motif database</xsl:text>
          </xsl:when>
          <xsl:when test="/tomtom/model/background/@from = 'file'">
            <xsl:value-of select="/tomtom/model/background/@file"/>
          </xsl:when>
        </xsl:choose>
        <xsl:text>):</xsl:text><br/>
        <div style="margin-left:25px;">
          <xsl:text>A:&nbsp;</xsl:text><xsl:value-of select="format-number(/tomtom/model/background/@A, '0.000')" />
          <xsl:text>&nbsp;&nbsp; </xsl:text>
          <xsl:text>C:&nbsp;</xsl:text><xsl:value-of select="format-number(/tomtom/model/background/@C, '0.000')" />
          <xsl:text>&nbsp;&nbsp; </xsl:text>
          <xsl:text>G:&nbsp;</xsl:text><xsl:value-of select="format-number(/tomtom/model/background/@G, '0.000')" />
          <xsl:text>&nbsp;&nbsp; </xsl:text>
          <xsl:text>T:&nbsp;</xsl:text><xsl:value-of select="format-number(/tomtom/model/background/@T, '0.000')" />
        </div><br />
        <xsl:text>Result calculation took </xsl:text><xsl:value-of select="/tomtom/runtime/@seconds"/><xsl:text> seconds</xsl:text>
        <br />
      </div>      
      <a href="javascript:showHidden('model')" id="model_activator">show model parameters...</a>
      <div class="subsection" id="model_data" style="display:none;">
        <h5>Model parameters</h5>
        <xsl:text>&newline;</xsl:text>
        <textarea style="width:100%;" rows="{count(/tomtom/model/*) - 1}" readonly="readonly">
          <xsl:variable name="spaces" select="'                    '"/>
          <xsl:text>&newline;</xsl:text>
          <xsl:for-each select="/tomtom/model/*[name(.) != 'command_line']">
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
      <a href="javascript:hideShown('model')" style="display:none;" id="model_deactivator">hide model parameters...</a>
    </div>
  </xsl:template>

  <xsl:template name="documentation">
    <span class="explain">
      <xsl:call-template name="header">
        <xsl:with-param name="title" select="'Explanation of TOMTOM Results'"/>
        <xsl:with-param name="self" select="'doc'"/>
        <xsl:with-param name="prev" select="'program'"/>
      </xsl:call-template>
    
      <div class="box">
        <h4>The TOMTOM results consist of</h4>
        <ul>
          <li>The <a href="#query_doc">query motifs</a> with links to further information on matches. [<a href="#query_motifs">View</a>]</li>
          <li>The <a href="#target_doc">target databases</a> with information on the number of matches. [<a href="#target_dbs">View</a>]</li>
          <li>The <b>matches</b> between query motifs and target motifs.  [<a href="#query_{/tomtom/queries/query_file[1]/query[1]/motif/@id}">View</a>]
            <ul>
              <li>The <a href="#summary_doc">summary</a> of match information.</li>
              <li>An <a href="#alignment_doc">alignment</a> of the two motifs.</li>
            </ul>
          </li>
          <li>The <b>program</b> details including:
            <ol>
              <li>The <b>version</b> of TOMTOM and the date it was released. [<a href="#version">View</a>]</li>
              <li>The <b>reference</b> to cite if you use TOMTOM in your research. [<a href="#reference">View</a>]</li>
              <li>The <b>command line summary</b> detailing the parameters with which you ran TOMTOM. [<a href="#command">View</a>]</li>
            </ol>
          </li>
          <li>This <b>explanation</b> of how to interpret TOMTOM results.</li>
        </ul>
        <h4>Inputs</h4>
        <p>The important inputs to TOMTOM.</p>
        <a name="query_doc"/>
        <h5>Query Motifs</h5>
        <div class="doc">
          <p>The query motifs section lists the motifs that were searched for in the target databases as well as links to matches found.</p>
          <dl>
            <a name="qm_name_doc"/>
            <dt>Name</dt>
            <dd>The motif name.</dd>
            <a name="qm_alt_doc"/>
            <dt>Alt. Name</dt>
            <dd>The alternative motif name. This column may be hidden if no query motifs have an alternate name.</dd>
            <a name="qm_web_doc"/>
            <dt>Website</dt>
            <dd>A link to a website which has more information on the motif. This column may be hidden if no query motifs have website information.</dd>
            <a name="qm_preview_doc"/>
            <dt>Preview</dt>
            <dd>The motif preview. On supporting browsers this will display as a motif logo, otherwise the consensus sequence will be displayed.</dd>
            <a name="qm_matches_doc"/>
            <dt>Matches</dt>
            <dd>The number of significant matches.</dd>
            <a name="qm_list_doc"/>
            <dt>List</dt>
            <dd>Links to the first 20 matches</dd>
          </dl>
        </div>
        <a name="target_doc"/>
        <h5>Target Databases</h5>
        <div class="doc">
          <p>The target motifs section lists details on the motif databases that were specifed to search.</p>
          <dl>
            <a name="tdb_database_doc"/>
            <dt>Database</dt>
            <dd>The database name.</dd>
            <a name="tdb_contains_doc"/>
            <dt>Number of Motifs</dt>
            <dd>The number of motifs read from the motif database minus the number that had to be discarded due to conflicting ids.</dd>
            <a name="tdb_matched_doc"/>
            <dt>Motifs Matched</dt>
            <dd>The number of motifs that had a match with at least one of the query motifs.</dd>
          </dl>
        </div>
        <h4>Matches</h4>
        <p>The matches section list the group of motifs that match each query motif. For each query-target match the following information is given:</p>
        <a name="summary_doc"/>
        <h5>Summary Table</h5>
        <div class="doc">
          <p>The summary gives the important statistics about the motif match.</p>
          <dl>
            <dt>Name</dt>
            <dd>The name of the matched motif.</dd>
            <dt>Alt. Name</dt>
            <dd>The alternative name of the matched motif.</dd>
            <dt><i>p</i>-value</dt>
            <dd>The probability that the match occurred by random chance according to the null model.</dd>
            <dt><i>E</i>-value</dt>
            <dd>The expected number of false positives in the matches up to this point.</dd>
            <dt><i>q</i>-value</dt>
            <dd>The minimum <b>F</b>alse <b>D</b>iscovery <b>R</b>ate required to include the match.</dd>
            <dt>Overlap</dt>
            <dd>The number of letters that overlaped in the optimal alignment.</dd>
            <dt>Offset</dt>
            <dd>The offset of the query motif to the matched motif in the optimal alignment.</dd>
            <dt>Orientation</dt>
            <dd>The orientation of the matched motif that gave the optimal alignment.
              A value of "normal" means that the matched motif is as it appears in 
              the database otherwise the matched motif has been reverse complemented.</dd>
          </dl>
        </div>
        <a name="alignment_doc"/>
        <h5>Alignment</h5>
        <div class="doc">
          <p>The image shows the alignment of the two motifs. The matched motif is shown on the top and the query motif is shown on the bottom.</p>
        </div>
        <a name="download_logo_doc"/>
        <h5>Create custom LOGO</h5>
        <div class="doc">
          <p>By clicking the link "Create custom LOGO &more;" a form to make custom logos 
            will be displayed. The download button can then be clicked to generate a motif
            matching the selected specifications.</p>
          <dl>
            <dt>Image Type</dt>
            <dd>Two image formats, png and eps, are avaliable. The pixel based portable network graphic (png) format is commonly used on the Internet 
              and the Encapsulated PostScript (eps) format is more suitable for publications that might require scaling.</dd>
            <dt>Error bars</dt>
            <dd>Toggle error bars indicating the confidence of a motif based on the number of sites used in its creation.</dd>
            <dt>SSC</dt>
            <dd>Toggle adding pseudocounts for <b>S</b>mall <b>S</b>ample <b>C</b>orrection.</dd>
            <dt>Flip</dt>
            <dd>Toggle a full reverse complement of the alignment.</dd>
            <dt>Width</dt>
            <dd>Specify the width of the generated motif.</dd>
            <dt>Height</dt>
            <dd>Specify the height of the generated motif.</dd>
          </dl>
        </div>
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
  
