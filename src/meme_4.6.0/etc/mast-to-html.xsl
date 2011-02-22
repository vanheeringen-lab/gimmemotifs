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
  <xsl:include href="meme.css.xsl"/>
  <xsl:include href="block-diagram.xsl"/>

  <!-- useful variables -->
  <xsl:variable name="max_seq_len"><!-- unfortunately I can't remove this calculation, maybe I should do this in c and store it? -->
    <xsl:for-each select="/mast/sequences/sequence" >
      <xsl:sort select="@length" data-type="number" order="descending"/>
      <xsl:if test="position() = 1">
        <xsl:value-of select="@length" />
      </xsl:if>
    </xsl:for-each>
  </xsl:variable>
  <xsl:variable name="max_log_pvalue">
    <xsl:call-template name="ln-approx">
      <xsl:with-param name="base" select="0.0000000001"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="longest_seq_name"><!-- unfortunately I can't remove this calculcation, maybe I should do this in C and store it? -->
    <xsl:for-each select="/mast/sequences/sequence" >
      <xsl:sort select="string-length(@name)" data-type="number" order="descending"/>
      <xsl:if test="position() = 1">
        <xsl:value-of select="string-length(@name)" />
      </xsl:if>
    </xsl:for-each>
  </xsl:variable>
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
      <meta charset="UTF-8"/>
      <meta name="description" content="Motif Alignment and Search Tool (MAST) output."/>
      <title>MAST</title>
      <script type="text/javascript">
        var motifs = new Array();<xsl:text>&newline;</xsl:text>
        <xsl:for-each select="/mast/motifs/motif">
          <xsl:text>        motifs[&quot;</xsl:text><xsl:value-of select='@id'/><xsl:text>&quot;] = new Motif(&quot;</xsl:text>
          <xsl:value-of select='@name'/><xsl:text>&quot;, &quot;</xsl:text>
          <xsl:value-of select='/mast/alphabet/@type'/><xsl:text>&quot;, &quot;</xsl:text>
          <xsl:value-of select='@best_f'/><xsl:text>&quot;, &quot;</xsl:text><xsl:value-of select='@best_r' /><xsl:text>&quot;);&newline;</xsl:text>
        </xsl:for-each>
        motifs.length = <xsl:value-of select="count(/mast/motifs/motif)"/>;

        var wrap = undefined;//size to display on one line
        var wrap_timer;

        var seqmax = <xsl:value-of select="$max_seq_len"/>;

        var loadedSequences = new Array();
        //draging details
        var moving_seq;
        var moving_annobox;
        var moving_left;
        var moving_width;
        var moving_both;

        //drag needles
        var dnl = null;
        var dnr = null;
        var drag_is_rc = undefined;
        //container
        var cont = null;
        
        function mouseCoords(ev) {
          ev = ev || window.event;
          if(ev.pageX || ev.pageY){ 
	          return {x:ev.pageX, y:ev.pageY}; 
	        } 
	        return { 
	          x:ev.clientX + document.body.scrollLeft - document.body.clientLeft, 
	          y:ev.clientY + document.body.scrollTop  - document.body.clientTop 
	        };
        }

        function setup() {
          rewrap();
          window.onresize = delayed_rewrap;
        }

        function calculate_wrap() {
          var est_wrap = 0;
          var page_width;
          if (window.innerWidth) {
            page_width = window.innerWidth;
          } else if (document.body) {
            page_width = document.body.clientWidth;
          } else {
            page_width = null;
          }
          if (page_width) {
            var ruler_width = document.getElementById("ruler").offsetWidth;
            est_wrap = Math.floor(page_width / ruler_width) * 5;
          }
          if (est_wrap > 20) {
            wrap = est_wrap - 20; 
          } else {
            wrap = 100;
          }
        }

        function rewrap() {
          var previous_wrap = wrap;
          calculate_wrap();
          if (previous_wrap != wrap) {
            for (var seqid in loadedSequences) {
              //exclude inherited properties and undefined properties
              if (!loadedSequences.hasOwnProperty(seqid) || loadedSequences[seqid] === undefined) continue;
            
              var sequence = loadedSequences[seqid];
              var annobox = document.getElementById(seqid + "_annotation");
              var leftPos = parseInt(document.getElementById(seqid + "_dnl").firstChild.firstChild.nodeValue);
              var rightPos = parseInt(document.getElementById(seqid + "_dnr").firstChild.firstChild.nodeValue);
            <xsl:choose>
              <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
                <xsl:text>            var fstrand = document.getElementById(seqid + "_fstrand");&newline;</xsl:text>
                <xsl:text>            var is_rc = !fstrand.checked;&newline;</xsl:text>
              </xsl:when>
              <xsl:otherwise>
                <xsl:text>            var is_rc = undefined;&newline;</xsl:text>
              </xsl:otherwise>
            </xsl:choose>
              set_data(leftPos, rightPos - leftPos + 1, is_rc, sequence, annobox);
            }
          }
          if (wrap_timer) {
            clearTimeout(wrap_timer);
          }
          wrap_timer = setTimeout("rewrap()", 5000);
        }

        function delayed_rewrap() {
          if (wrap_timer) {
            clearTimeout(wrap_timer);
          }
          wrap_timer = setTimeout("rewrap()", 1000);
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


        function show_more(seqid) {
          if (wrap == undefined) rewrap();
          var seq_row = document.getElementById(seqid + "_blocks");
          var info_row = document.getElementById(seqid + "_info");
          var tbl = document.getElementById("tbl_sequences");
          if (info_row) {
            tbl.deleteRow(info_row.rowIndex);
            destroy_seq_handle(seqid, true);
            destroy_seq_handle(seqid, false);
            delete loadedSequences[seqid];
          } else {
            var sequence = new Sequence(seqid);
            loadedSequences[seqid] = sequence;
            info_row = tbl.insertRow(seq_row.rowIndex + 1);
            info_row.id = seqid + "_info";
            var cell = info_row.insertCell(0);
        <xsl:choose>
          <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
            cell.colSpan = 6;
          </xsl:when>
          <xsl:otherwise>
            cell.colSpan = 4;
          </xsl:otherwise>
        </xsl:choose>
            cell.style.verticalAlign = "top";
            var width = Math.min(wrap, sequence.length); 
            var start = 1;
            var is_neg_strand = false;
            //center on the first hit
            if (sequence.hits.length > 0) {
              var hit = sequence.hits[0];
              is_neg_strand = hit.is_rc;
              start = Math.max(1, hit.pos + Math.round(hit.width / 2) - Math.round(width / 2));
              if ((start + width - 1) > sequence.length) {
                start = sequence.length - width + 1;
              }
            }
            create_info_section(cell, sequence, seqid, start, width, is_neg_strand);
            var blockCont = document.getElementById(seqid + "_block_container");
            //create sequence view handles
            create_seq_handle(blockCont, seqid, true, start, seqmax);
            create_seq_handle(blockCont, seqid, false, start + width-1, seqmax);
          }
        }

        function create_info_section(cell, sequence, seqid, start, width, is_neg_strand) {
          var box = document.createElement('div');
          box.className = "box infobox";
          //generic info
          box.appendChild(document.createTextNode("Change the portion of annotated sequence by "));
          var bold_box = document.createElement('b');
          bold_box.appendChild(document.createTextNode("dragging the buttons"));
          box.appendChild(bold_box);
          box.appendChild(document.createTextNode("; hold shift to drag them individually."));
          box.appendChild(document.createElement('br'));
          box.appendChild(document.createElement('br'));
          //description
          var descTitle = document.createElement('h5');
          descTitle.appendChild(document.createTextNode("Comment:"));
          descTitle.className = "inlineTitle";
          box.appendChild(descTitle);
          box.appendChild(document.createTextNode(" " + sequence.desc));
          box.appendChild(document.createElement('br'));
          box.appendChild(document.createElement('br'));
          var cpvalueTitle = document.createElement('h5');
          cpvalueTitle.appendChild(document.createTextNode("Combined "));
          var italicp = document.createElement('i');
          italicp.appendChild(document.createTextNode("p"));
          cpvalueTitle.appendChild(italicp);
          cpvalueTitle.appendChild(document.createTextNode("-value:"));
          cpvalueTitle.className = "inlineTitle";
          box.appendChild(cpvalueTitle);
          if (sequence.cpvalue.indexOf("/") != -1) {
            var pos = sequence.cpvalue.indexOf("/");
            var forward = sequence.cpvalue.substring(0, pos);
            var reverse = sequence.cpvalue.substring(pos+1);
            var dim_forward;
            if (forward.indexOf("--") != -1) {
              dim_forward = true;
            } else if (reverse.indexOf("--") != -1) {
              dim_forward = false;
            } else {
              if ( (+forward) &gt; (+reverse)) {
                dim_forward = true;
              } else {
                dim_forward = false;
              }
            }
            var span = document.createElement('span');
            span.className = "dim";
            if (dim_forward) {// dim the forward value
              span.appendChild(document.createTextNode(forward));
              box.appendChild(document.createTextNode(" "));
              box.appendChild(span);
              box.appendChild(document.createTextNode("/"));
              box.appendChild(document.createTextNode(reverse));
            } else {// dim the reverse value
              span.appendChild(document.createTextNode(reverse));
              box.appendChild(document.createTextNode(" "));
              box.appendChild(document.createTextNode(forward));
              box.appendChild(document.createTextNode("/"));
              box.appendChild(span);
            }
          } else {
            box.appendChild(document.createTextNode(" " + sequence.cpvalue));
          }
          <xsl:if test="/mast/model/translate_dna/@value = 'y'">
          box.appendChild(document.createElement('br'));
          box.appendChild(document.createElement('br'));
          var bestframeTitle = document.createElement('h5');
          bestframeTitle.appendChild(document.createTextNode("Best frame:"));
          bestframeTitle.className = "inlineTitle";
          box.appendChild(bestframeTitle);
          if (sequence.frame.indexOf("/") != -1) {
            var pos = sequence.frame.indexOf("/");
            var forward = sequence.frame.substring(0, pos);
            var reverse = sequence.frame.substring(pos+1);
            var dim_forward;
            if (forward.indexOf("--") != -1) {
              dim_forward = true;
            } else if (reverse.indexOf("--") != -1) {
              dim_forward = false;
            } else {
              if ( (+forward) &gt; (+reverse)) {
                dim_forward = true;
              } else {
                dim_forward = false;
              }
            }
            var span = document.createElement('span');
            span.className = "dim";
            if (dim_forward) {// dim the forward value
              span.appendChild(document.createTextNode(forward));
              box.appendChild(document.createTextNode(" "));
              box.appendChild(span);
              box.appendChild(document.createTextNode("/"));
              box.appendChild(document.createTextNode(reverse));
            } else {// dim the reverse value
              span.appendChild(document.createTextNode(reverse));
              box.appendChild(document.createTextNode(" "));
              box.appendChild(document.createTextNode(forward));
              box.appendChild(document.createTextNode("/"));
              box.appendChild(span);
            }
          } else {
            box.appendChild(document.createTextNode(" " + sequence.frame));
          }
          </xsl:if>
          box.appendChild(document.createElement('br'));
          box.appendChild(document.createElement('br'));
          //sequence display
          var seqDispTitle = document.createElement('h5');
          seqDispTitle.appendChild(document.createTextNode("Annotated Sequence"));
          box.appendChild(seqDispTitle);
          <xsl:choose>
            <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
          var show_rc_only = is_neg_strand;
          var myform = document.createElement('form');
          var fstrand = document.createElement('input');
          fstrand.id = seqid + "_fstrand";
          fstrand.type = "radio";
          fstrand.name = seqid + "_strand";
          fstrand.onclick = function(evt) {
            handle_strand_change(seqid);
          };
          var rstrand = document.createElement('input');
          rstrand.id = seqid + "_rstrand";
          rstrand.type = "radio";
          rstrand.name = seqid + "_strand";
          rstrand.onclick = function(evt) {
            handle_strand_change(seqid);
          };
          fstrand.checked = !show_rc_only;
          rstrand.checked = show_rc_only;
          myform.appendChild(fstrand);
          myform.appendChild(document.createTextNode("forward strand"));
          myform.appendChild(rstrand);
          myform.appendChild(document.createTextNode("reverse strand"));
          box.appendChild(myform);
            </xsl:when>
            <xsl:otherwise>
          var show_rc_only = undefined;
            </xsl:otherwise>
          </xsl:choose>
          var annobox = document.createElement('div');
          annobox.id = seqid + "_annotation";
          set_data(start, width, show_rc_only, sequence, annobox);
          box.appendChild(annobox);
          cell.appendChild(box);
        }

        function annobox_labels() {
          var seqDispLabel = document.createElement('div');
          seqDispLabel.className = "sequence sequence_labels";
          seqDispLabel.style.height = "2.5em";
          return seqDispLabel;
        }

        function annobox_hits() {
          var seqDispHit = document.createElement('div');
          seqDispHit.className = "sequence";
          seqDispHit.style.height = "1.5em";
          return seqDispHit;
        }

        function annobox_matches() {
          var seqDispMatch = document.createElement('div');
          seqDispMatch.className = "sequence";
          seqDispMatch.style.height = "1.5em";
          return seqDispMatch;
        }

        function annobox_translations() {
          var seqDispXlate = document.createElement('div');
          seqDispXlate.className = "sequence";
          seqDispXlate.style.height = "1.5em";
          return seqDispXlate;
        }

        function annobox_sequence() {
          var seqDispSeq = document.createElement('div');
          seqDispSeq.className = "sequence";
          seqDispSeq.style.height = "1.5em";
          return seqDispSeq; 
        }
        function annobox_boundary() {
          var seqDispSeq = document.createElement('div');
          seqDispSeq.className = "sequence";
          seqDispSeq.style.height = "1.5em";
          return seqDispSeq; 
        }

        function create_seq_handle(container, seqid, isleft, pos, max) {
          var vbar = document.createElement('div');
          vbar.id = seqid + "_dn" + (isleft ? "l" : "r");
          vbar.className = "block_needle";
          //the needles sit between the sequence positions, so the left one sits at the start and the right at the end
          //this is why 1 is subtracted off the position for a left handle
          vbar.style.left = "" + ((isleft ? (pos - 1) : pos) / max * 100)+ "%";
          var label = document.createElement('div');
          label.className = "block_handle";
          label.appendChild(document.createTextNode("" + pos));
          label.style.cursor = 'pointer';
          label.title = "Drag to move the displayed range. Hold shift and drag to change " + (isleft ? "lower" : "upper") + " bound of the range.";
          vbar.appendChild(label);
          container.appendChild(vbar);
          label.onmousedown = function(evt) {
            evt = evt || window.event;
            start_drag(seqid, isleft, !(evt.shiftKey));
          };
        }

        function destroy_seq_handle(seqid, isleft) {
          var vbar = document.getElementById(seqid + "_dn" + (isleft ? "l" : "r"));
          var label = vbar.firstChild;
          label.onmousedown = null;
          vbar.parentNode.removeChild(vbar);
        }

        <xsl:if test="/mast/model/strand_handling/@value = 'separate'">
        function handle_strand_change(seqid) {
          var left = document.getElementById(seqid + "_dnl");
          var right = document.getElementById(seqid + "_dnr");
          var leftPos = parseInt(left.firstChild.firstChild.nodeValue);
          var rightPos = parseInt(right.firstChild.firstChild.nodeValue);
          var fstrand = document.getElementById(seqid + "_fstrand");
          var show_rc_only = !fstrand.checked;
          var data = loadedSequences[seqid];
          if (!data) {
            alert("Sequence not loaded?!");
          }
          var annobox = document.getElementById(seqid + "_annotation");
          set_data(leftPos, rightPos - leftPos + 1, show_rc_only, data, annobox);
        }
        </xsl:if>

        function start_drag(seqid, mouse_on_left, move_both) {
          //first record what we are moving
          moving_left = mouse_on_left;
          moving_both = move_both;
          moving_seq = loadedSequences[seqid];
          if (!moving_seq) {
            alert("Sequence not loaded?!");
          }
          //get the container for movements to be judged against
          cont = document.getElementById(seqid + "_block_container");
          //get the div's that will be moved
          dnl = document.getElementById(seqid + "_dnl");
          dnr = document.getElementById(seqid + "_dnr");
          //get the div which has all the text containers that will be updated
          moving_annobox = document.getElementById(seqid + "_annotation");
          //calculate the space between handles
          moving_width = dnr.firstChild.firstChild.nodeValue - 
              dnl.firstChild.firstChild.nodeValue;
          <xsl:text>&newline;</xsl:text>
          <xsl:choose>
            <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
              <xsl:text>          var fstrand = document.getElementById(seqid + "_fstrand");&newline;</xsl:text>
              <xsl:text>          drag_is_rc = !fstrand.checked;&newline;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>          drag_is_rc = undefined;&newline;</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
          //setup the events for draging
          document.onmousemove = handle_drag;
          document.onmouseup = finish_drag;
        }


        function calculate_width(obj) {
          return obj.clientWidth - (obj.style.paddingLeft ? obj.style.paddingLeft : 0) 
              - (obj.style.paddingRight ? obj.style.paddingRight : 0);
        }

        function calculate_drag_pos(ev) {
          var mouse = mouseCoords(ev);
          var dragable_length = calculate_width(cont);
          //I believe that the offset parent is the body
          //otherwise I would need to make this recursive
          //maybe clientLeft would work, but the explanation of
          //it is hard to understand and it apparently doesn't work
          //in firefox 2.
          var diff = mouse.x - cont.offsetLeft;
          if (diff &lt; 0) diff = 0;
          if (diff &gt; dragable_length) diff = dragable_length;
          var pos = Math.round(diff / dragable_length * (seqmax));
          return pos + 1;
        }

        function handle_drag(ev) {
          var pos = calculate_drag_pos(ev);
          update_needles(pos, moving_seq.length, false);
        }

        function finish_drag(ev) {
          document.onmousemove = null;
          document.onmouseup = null;
          var pos = calculate_drag_pos(ev);
          update_needles(pos, moving_seq.length, true);
        }

        function update_needles(pos, seqlen, updateTxt) {
          var leftPos = parseInt(dnl.firstChild.firstChild.nodeValue);
          var rightPos = parseInt(dnr.firstChild.firstChild.nodeValue);
          if (moving_both) {
            if (moving_left) {
              if ((pos + moving_width) &gt; seqlen) {
                pos = seqlen - moving_width;
              }
              leftPos = pos;
              rightPos = pos + moving_width;
            } else {
              if (pos &gt; seqlen) {
                pos = seqlen;
              } else if ((pos - moving_width) &lt; 1) {
                pos = moving_width + 1;
              }
              leftPos = pos - moving_width;
              rightPos = pos;
            }
          } else {
            if (moving_left) {
              if (pos &gt; rightPos) {
                pos = rightPos;
              }
              leftPos = pos;
            } else {
              if (pos &lt; leftPos) {
                pos = leftPos;
              } else if (pos &gt; seqlen) {
                pos = seqlen;
              }
              rightPos = pos;
            }
          }
          set_needles(seqmax, dnl, leftPos, dnr, rightPos);
          if (updateTxt) set_data(leftPos, rightPos - leftPos + 1, drag_is_rc, moving_seq, moving_annobox);
        }

        function set_needles(max, left, leftPos, right, rightPos) {
          //the needles are between the sequence positions
          //the left needle is before and the right needle after
          //think of it as an inclusive range...
          left.style.left = "" + ((leftPos -1) / max * 100)+ "%";
          left.firstChild.firstChild.nodeValue="" + leftPos;
          right.style.left = "" + (rightPos / max * 100)+ "%";
          right.firstChild.firstChild.nodeValue="" + rightPos;
        }

        function set_data(start, width, show_rc_only, data, annobox) {
          var child = annobox.firstChild;
          var line_width = Math.min(wrap, width);
          var num_per_wrap = <xsl:choose><xsl:when test="/mast/model/translate_dna/@value = 'y'">6</xsl:when><xsl:otherwise>5</xsl:otherwise></xsl:choose>;
          var end = start + width;
          for (var i = start; i &lt; end; i += line_width, line_width = Math.min(wrap, end - i)) {
            for (var j = 0; j &lt; num_per_wrap; ++j) {
              if (child) {
                while (child.firstChild) child.removeChild(child.firstChild);
              } else {
                switch (j) {
                case 0:
                  child = annobox_labels();
                  break;
                case 1:
                  child = annobox_hits();
                  break;
                case 2:
                  child = annobox_matches();
                  break;
                case 3:
      <xsl:choose>
        <xsl:when test="/mast/model/translate_dna/@value = 'y'">
                  child = annobox_translations();
                  break;
                case 4:
                  child = annobox_sequence();
                  break;
                case 5:
        </xsl:when>
        <xsl:otherwise>
                  child = annobox_sequence();
                  break;
                case 4:
        </xsl:otherwise>
      </xsl:choose>
                  child = annobox_boundary();
                  break;
                }
                annobox.appendChild(child);
              }
              switch (j) {
              case 0:
                data.append_labels(child, i, line_width, show_rc_only);
                break;
              case 1:
                data.append_hits(child, i, line_width, show_rc_only);
                break;
              case 2:
                data.append_matches(child, i, line_width, show_rc_only);
                break;
              case 3:
      <xsl:choose>
        <xsl:when test="/mast/model/translate_dna/@value = 'y'">
                data.append_translation(child, i, line_width, show_rc_only);
                break;
              case 4:
                data.append_seq(child, i, line_width, show_rc_only);
                break;
              case 5:
        </xsl:when>
        <xsl:otherwise>
                data.append_seq(child, i, line_width, show_rc_only);
                break;
              case 4:
        </xsl:otherwise>
      </xsl:choose>
                data.append_boundary(child, i, line_width, show_rc_only);
                break;
              }
              child = child.nextSibling;
            }
          }
          //clean up excess
          
          while (child) {
            var next = child.nextSibling;
            annobox.removeChild(child);
            child = next;
          }
        }

        function append_coloured_nucleotide_sequence(container, sequence) {
          var len = sequence.length;
          for (var i = 0; i &lt; len; i++) {
            var colour = "black";
            switch (sequence.charAt(i)) {
            case 'A':
              colour = "red";
              break;
            case 'C':
              colour = "blue";
              break;
            case 'G':
              colour = "orange";
              break;
            case 'T':
              colour = "green";
              break;
            }
            var letter = document.createElement('span');
            letter.style.color = colour;
            letter.appendChild(document.createTextNode(sequence.charAt(i)));
            container.appendChild(letter);
          }
        }

        function append_coloured_peptide_sequence(container, sequence) {
          var len = sequence.length;
          for (var i = 0; i &lt; len; i++) {
            var colour = "black";
            switch (sequence.charAt(i)) {
            case 'A':
            case 'C':
            case 'F':
            case 'I':
            case 'L':
            case 'V':
            case 'W':
            case 'M':
              colour = "blue";
              break;
            case 'N':
            case 'Q':
            case 'S':
            case 'T':
              colour = "green";
              break;
            case 'D':
            case 'E':
              colour = "magenta";
              break;
            case 'K':
            case 'R':
              colour = "red";
              break;
            case 'H':
              colour = "pink";
              break;
            case 'G':
              colour = "orange";
              break;
            case 'P':
              colour = "yellow";
              break;
            case 'Y':
              colour = "turquoise";
              break;
            }
            var letter = document.createElement('span');
            letter.style.color = colour;
            letter.appendChild(document.createTextNode(sequence.charAt(i)));
            container.appendChild(letter);
          }
        }

        /*
         * Finds the index of an item in a sorted array (a) that equals a key (k)
         * when compared using a comparator (c). If an exact match can not be found
         * then returns -(index+1), where index is the location that the item would be inserted.
         */
        function bsearch(a, c, k) {
          var low = 0;
          var high = a.length;
          while (low &lt; high) {
            var mid = low + Math.floor((high - low) /2);
            if (c(a[mid], k) &lt; 0) {
              low = mid + 1;
            } else {
              high = mid;
            }
          }
          if ((low &lt; a.length) &amp;&amp; (c(a[low], k) == 0)) {
            return low;
          } else {
            return -(low+1);
          }
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Start Motif Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function Motif(name, type, best_f, best_r) {
          this.name = name;
          this.is_nucleotide = (type === "nucleotide");
          this.best_f = best_f;
          if (best_r.length == 0) {
            var txt = best_f.split("");
            txt.reverse();
            best_r = txt.join("");
          }
          this.best_r = best_r;
          this.length = best_f.length;
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // End Motif Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Start Seg Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function Seg(start, segment) {
          this.start = parseInt(start);
          this.segment = segment;
        }

        function compare_seg_to_pos(seg, pos) {
          return seg.start - pos;
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // End Seg Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Start Hit Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function Hit(sequence_is_nucleotide, motif, pos, strand, pvalue, match, xlated) {
          //properties
          this.motif = motif;
          this.pos = parseInt(pos);
          this.width = 0;
          this.is_rc = (strand != "+");
          this.pvalue = pvalue;
          this.xlated = "";
          this.best = "";
          this.match = "";
          //functions
          this.find_overlap = Hit_find_overlap;
          //setup
          var seq;
          var m = match.replace(/ /g, "&nbsp;");
          if (this.is_rc) seq = this.motif.best_r;
          else seq = this.motif.best_f;
          if (sequence_is_nucleotide == this.motif.is_nucleotide) {
            this.best = seq;
            this.match = m;
            this.xlated = xlated;
            this.width = motif.length;
          } else if (sequence_is_nucleotide) {
            this.width = motif.length * 3;
            for (var i = 0; i &lt; motif.length; i++) {
              this.best += seq.charAt(i);
              this.best += "..";
              this.xlated += xlated.charAt(i);
              this.xlated += "..";
              this.match += m.charAt(i);
              this.match += "&nbsp;&nbsp;";
            }
          } else {
            throw "UNSUPPORTED_CONVERSION";
          }
        }

        function Hit_find_overlap(sequence_is_nucleotide, start, width) {
          if (this.pos &lt; start) {
            if ((this.pos + this.width) &gt; start) {
              return {hit: this, start: start, length: Math.min(width, (this.pos + this.width - start))};
            } else {
              return null;
            }
          } else if (this.pos &lt; (start + width)) {
            return {hit: this, start: this.pos, length: Math.min(this.width, (start + width - this.pos))};
          }
          return null;
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // End Hit Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Start Overlapping Hit Iterator
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function OLHIterator(sequence, start, width, is_rc) {
          //properties
          this.sequence = sequence;
          this.start = start;
          this.width = width;
          this.is_rc = is_rc;
          this.both = (is_rc === undefined);
          this.index = 0;
          this.nextO = null;
          //methods
          this.next = OLHIterator_next;
          this.has_next = OLHIterator_has_next;
          //setup
          // find the first hit which overlaps the range
          for (; this.index &lt; this.sequence.hits.length; this.index++) {
            var hit = this.sequence.hits[this.index];
            if (!this.both &amp;&amp; this.is_rc != hit.is_rc) continue;
            if ((this.nextO = hit.find_overlap(this.sequence.is_nucleotide, this.start, this.width)) != null) break;
          }
        }
        
        function OLHIterator_next() {
          if (this.nextO == null) throw "NO_NEXT_ELEMENT";
          var current = this.nextO;
          this.nextO = null;
          for (this.index = this.index + 1; this.index &lt; this.sequence.hits.length; this.index++) {
            var hit = this.sequence.hits[this.index];
            if (!this.both &amp;&amp; this.is_rc != hit.is_rc) continue;
            this.nextO = hit.find_overlap(this.sequence.is_nucleotide, this.start, this.width);
            break;
          }
          return current;
        }
        
        function OLHIterator_has_next() {
          return (this.nextO != null);
        }

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // End Overlapping Hit Iterator
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Start Sequence Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        function Sequence(seqid) {
          //properties
          this.length = parseInt(document.getElementById(seqid + "_len").value);
          this.desc = document.getElementById(seqid + "_desc").value;
          this.cpvalue = document.getElementById(seqid + "_combined_pvalue").value;
          <xsl:if test="/mast/model/translate_dna/@value = 'y'">
          this.frame = document.getElementById(seqid + "_frame").value;
          </xsl:if>
          this.is_nucleotide = (document.getElementById(seqid + "_type").value === "nucleotide");
          this.segs = new Array(); //sorted list of segs
          this.hits = new Array(); //sorted list of hits
          //functions
          this.get_overlapping_hits = Sequence_get_overlapping_hits;
          this.append_seq = Sequence_append_seq;
          this.append_translation = Sequence_append_translation;
          this.append_hits = Sequence_append_hits;
          this.append_matches = Sequence_append_matches;
          this.append_labels = Sequence_append_labels;
          this.append_boundary = Sequence_append_boundary;
          //init
          //made this parser much more permissive as a new
          //version of libxml2 broke the translate call I was
          //using to remove whitespace characters
          var mysegs = document.getElementById(seqid + "_segs");
          var tokens = mysegs.value.split(/\s+/);
          var offsetnum = "";
          var seqchunk = "";
          for (var i = 0; i &lt; tokens.length; i++) {
            var token = tokens[i];
            if (token == "") continue;
            if (token.match(/\d+/)) {
              if (offsetnum != "" &amp;&amp; seqchunk != "") {
                this.segs.push(new Seg(offsetnum, seqchunk));
              }
              offsetnum = parseInt(token, 10);
              seqchunk = "";
              continue;
            }
            seqchunk += token;
          }
          if (offsetnum != "" &amp;&amp; seqchunk != "") {
            this.segs.push(new Seg(offsetnum, seqchunk));
          }
          var myhits = document.getElementById(seqid + "_hits");
          lines = myhits.value.split(/\n/);
          for (var i in lines) {
            //exclude inherited properties and undefined properties
            if (!lines.hasOwnProperty(i) || lines[i] === undefined) continue;
            
            var line = lines[i];
            var chunks = line.split(/\t/);
            if (chunks.length != 6) continue;
            var pos = chunks[0];
            var motif = motifs[chunks[1]];
            var strand = chunks[2];
            var pvalue = chunks[3];
            var match = chunks[4];
            var xlated = chunks[5];
            this.hits.push(new Hit(this.is_nucleotide, motif, pos, strand, pvalue, match, xlated));
          }
        }

        function Sequence_get_overlapping_hits(start, width, is_rc) {
          return new OLHIterator(this, start, width, is_rc);
        }

        function Sequence_append_seq(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //make a sub container to put the sequence in
          var mycontainer = document.createElement('span');
          // find where a seg starting at start would be
          var i = bsearch(this.segs, compare_seg_to_pos, start);
          var seg;
          var sequence = "";
          if (i &lt; -1) {
            //possible partial segment, need to handle first
            i = -(i + 1);
            seg = this.segs[i-1];
            if ((seg.start + seg.segment.length) &gt; start) {
              var seg_offset = start - seg.start;
              var seg_width = Math.min(width, seg.segment.length - seg_offset);
              sequence += seg.segment.substring(seg_offset, seg_offset + seg_width);
            }
          } else if (i == -1) {
            //gap, with following segment, normal state at start of iteration
            i = 0;
          } 
          for (;i &lt; this.segs.length; i++) {
            seg = this.segs[i];
            var gap_width = Math.min(width - sequence.length, seg.start - (start + sequence.length));
            for (; gap_width &gt; 0; gap_width--) sequence += '-';
            var seg_width = Math.min(width - sequence.length, seg.segment.length);
            if (seg_width == 0) break;
            sequence += seg.segment.substring(0, seg_width);
          }
          while (sequence.length &lt; width) sequence += '-';
          // calculate which parts are aligned with hits and output them in colour
          var pos = start;
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var from, to;
            if (o.start &gt; pos) {
              from = pos - start;
              to = o.start - start;
              var subseq = sequence.substring(from, to);
              var greytxt = document.createElement('span');
              greytxt.style.color = 'grey';
              greytxt.appendChild(document.createTextNode(subseq));
              mycontainer.appendChild(greytxt);
              pos = o.start;
            }
            from = pos - start;
            to = from + o.length;
            var motifseq = sequence.substring(from, to);
            if (this.is_nucleotide) {
              append_coloured_nucleotide_sequence(mycontainer, motifseq);
            } else {
              append_coloured_peptide_sequence(mycontainer, motifseq);
            }
            pos += o.length;
          }
          if (pos &lt; (start + width)) {
            var greytxt = document.createElement('span');
            greytxt.style.color = 'grey';
            greytxt.appendChild(document.createTextNode(sequence.substring(pos - start)));
            mycontainer.appendChild(greytxt);
          }
          container.appendChild(mycontainer);
        }

        function Sequence_append_translation(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //make a sub container to put the sequence in
          var mycontainer = document.createElement('span');
          // calculate which parts are aligned with hits and output them in colour
          var pos = start;
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var space = "";
            var from, to;
            while (o.start &gt; pos) {
              space += "&nbsp;";
              ++pos;
            }
            if (space.length &gt; 0) {
              var spacer = document.createElement('span');
              spacer.appendChild(document.createTextNode(space));
              mycontainer.appendChild(spacer);
            }
            from = o.start - o.hit.pos;
            to = from + o.length;
            var motifseq = o.hit.xlated.substring(from, to);
            if (o.hit.motif.is_nucleotide) {
              append_coloured_nucleotide_sequence(mycontainer, motifseq);
            } else {
              append_coloured_peptide_sequence(mycontainer, motifseq);
            }
            pos += o.length;
          }

          container.appendChild(mycontainer);
        }

        function Sequence_append_matches(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //make a sub container to put the sequence in
          var mycontainer = document.createElement('span');
          var pos = start;
          var text = "";
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var space = "";
            var from, to;
            while (o.start &gt; pos) {
              text += "&nbsp;";
              ++pos;
            }
            from = o.start - o.hit.pos;
            to = from + o.length;
            var motifseq = o.hit.match.substring(from, to);
            text += motifseq
            pos += o.length;
          }
          mycontainer.appendChild(document.createTextNode(text));
          container.appendChild(mycontainer);
        }

        function Sequence_append_hits(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //make a sub container to put the sequence in
          var mycontainer = document.createElement('span');
          // calculate which parts are aligned with hits and output them in colour
          var pos = start;
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var space = "";
            var from, to;
            while (o.start &gt; pos) {
              space += "&nbsp;";
              ++pos;
            }
            if (space.length &gt; 0) {
              var spacer = document.createElement('span');
              spacer.appendChild(document.createTextNode(space));
              mycontainer.appendChild(spacer);
            }
            from = o.start - o.hit.pos;
            to = from + o.length;
            var motifseq = o.hit.best.substring(from, to);
            if (o.hit.motif.is_nucleotide) {
              append_coloured_nucleotide_sequence(mycontainer, motifseq);
            } else {
              append_coloured_peptide_sequence(mycontainer, motifseq);
            }
            pos += o.length;
          }

          container.appendChild(mycontainer);
        }

        function Sequence_append_labels(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //create a table (ye gods...)
          var oTable = document.createElement("table");
          var oRow = oTable.insertRow(oTable.rows.length);

          // calculate where to put the labels
          var pos = start;
          var cellindex = 0;
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var motif_center = Math.floor(o.start + o.length / 2);
            var spacer = "";
            while (pos &lt; motif_center) {
              spacer += "&nbsp;";
              pos++;
            }
            var cell1 = oRow.insertCell(cellindex++);
            cell1.appendChild(document.createTextNode(spacer));
            //add one for the div
            pos++;
            var cell2 = oRow.insertCell(cellindex++);
            //create the div that holds the label divs
            var div_container = document.createElement('div');
            div_container.className = "block_container";
            cell2.appendChild(div_container);
            //create the top label
            var top_label = "";
            <xsl:variable name="has_frame" select="/mast/model/translate_dna/@value = 'y'"/>
            <xsl:variable name="db_dna" select="/mast/sequences/database/@type = 'nucleotide'"/>
            <xsl:variable name="has_strand" select="$db_dna and /mast/model/strand_handling/@value != 'separate'"/>
            <xsl:if test="$has_strand or $has_frame">
              top_label += "(";
              <xsl:if test="$has_strand">
                if (o.hit.is_rc) {
                  top_label += "-";
                } else {
                  top_label += "+";
                }
              </xsl:if>
              <xsl:if test="$has_frame">
                switch ((o.hit.pos - 1) % 3) {
                  case 0:
                  top_label += "a";
                  break;
                  case 1:
                  top_label += "b";
                  break;
                  default:
                  top_label += "c";
                  break;
                }
              </xsl:if>
              top_label += ") ";
            </xsl:if>
            if (o.hit.motif.name.match(/^\d+$/)) {
              top_label += "Motif " + o.hit.motif.name;
            } else {
              top_label += o.hit.motif.name;
            }
            var div_top_label = document.createElement('div');
            div_top_label.className = "sequence_label sequence_label_top";
            div_top_label.style.width = "" + top_label.length + "em";
            div_top_label.style.left = "-" + (top_label.length / 2) + "em";
            div_top_label.appendChild(document.createTextNode(top_label));
            //create the bottom label
            var bottom_label = o.hit.pvalue;
            var div_bottom_label = document.createElement('div');
            div_bottom_label.className = "sequence_label sequence_label_bottom";
            div_bottom_label.style.width = "" + bottom_label.length + "em";
            div_bottom_label.style.left = "-" + (bottom_label.length / 2) + "em";
            div_bottom_label.appendChild(document.createTextNode(bottom_label));
            //add the label divs to the positioned container
            div_container.appendChild(div_top_label);
            div_container.appendChild(div_bottom_label);
          }
          container.appendChild(oTable);

        }

        function Sequence_append_boundary(container, start, width, is_rc) {
          if (start &gt; this.length) {
            alert("start: " + start + " length: " + this.length);
            throw "INDEX_OUT_OF_BOUNDS";
          }
          if ((start + width - 1) &gt; this.length) {
            alert("start: " + start + " width: " + width + " length: " + this.length);
            throw "RANGE_OUT_OF_BOUNDS";
          }
          //make a sub container to put the sequence in
          var mycontainer = document.createElement('span');
          var pos = start;
          var text = "";
          var iter = this.get_overlapping_hits(start, width, is_rc);
          while (iter.has_next()) {
            var o = iter.next();
            var end;
            while (o.start &gt; pos) {
              text += "&nbsp;";
              ++pos;
            }
            if (o.start == o.hit.pos) {
              text += "\\";
              ++pos;
            }
            end = o.start + o.length - 1;
            while (end &gt; pos) {
              text += "_";
              ++pos;
            }
            if (end == (o.hit.pos + o.hit.width -1)) {
              text += "/"
            } else {
              text += "_"
            }
            ++pos;
          }
          mycontainer.appendChild(document.createTextNode(text));
          container.appendChild(mycontainer);
        }
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // End Sequence Object
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      </script>
      <style type="text/css">
        <xsl:call-template name="meme.css" />
        .more_arrow {
          font-family:Arial,Serif;
          font-size: larger;
          font-weight: bold;
          text-decoration:none; 
        }
        div.infobox {
          background-color:#ddddff;
          margin-top: 1.6em;
          margin-bottom: 1em;
        }
        .sequence {font-family:Monospace;}
        .sequence_labels {position:relative; width:100%; height:2.5em; 
            padding:0px; margin:0px;}
        .sequence_label {font-family:Serif; position:absolute; z-index:2; 
            height:1em; text-align:center; vertical-align:middle;}
        .sequence_label_top {top:0px;}
        .sequence_label_bottom {top:1.25em;}
        .inlineTitle {display:inline;}
        .block_needle {position:absolute; z-index:4; height:30px; width:1px; 
            top:-2px; background-color:gray;}
        .block_handle {position:absolute; z-index:5; height:1.1em; width:3em; 
            top:30px; left:-1.5em; text-align:center; vertical-align:middle;
            background-color: LightGrey; border:3px outset grey;}
        table.padded-table td { padding:0px 10px; }
        table.padded-table th { padding:0px 5px; }
        td.tnum {text-align:right;}
        tr.highlight {background:#aaffaa;}
        td.dim {color: gray;}
        span.dim {color: gray;}
      </style>
    </head>
  </xsl:template>

  <xsl:template name="html-body">
    <body onload="javascript:setup()" >
      <xsl:call-template name="top"/>
      <xsl:call-template name="navigation"/>
      <xsl:call-template name="inputs"/>
      <xsl:call-template name="results"/>
      <xsl:call-template name="program"/>
      <xsl:call-template name="documentation"/>
      <xsl:call-template name="data"/>
      <xsl:call-template name="footer"/>
      <span class="sequence" id="ruler" style="visibility:hidden; white-space:nowrap;">ACGTN</span>
    </body>
  </xsl:template>

  <xsl:template name="top">
    <a name="top"/>
    <div class="pad1">
      <h1><img src="http://meme.nbcr.net/meme/images/mast.png" alt="Motif Alignment and Search Tool (MAST)" /></h1>
      <p class="spaced">
        For further information on how to interpret these results or to get a 
        copy of the MEME software please access 
        <a href="http://meme.nbcr.net/">http://meme.nbcr.net</a>. 
      </p>
    </div>
  </xsl:template>

  <xsl:template name="navigation">
    <a name="navigation"/>
    <div class="pad2">
      <a class="jump" href="#inputs">Inputs</a>
      &nbsp;&nbsp;|&nbsp;&nbsp;
      <a class="jump" href="#results">Search Results</a>
      &nbsp;&nbsp;|&nbsp;&nbsp;
      <a class="jump" href="#program">Program information</a>
      &nbsp;&nbsp;|&nbsp;&nbsp;
      <a class="jump" href="#doc">Documentation</a>
    </div>    
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

  <xsl:template name="inputs" >
    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Inputs'"/>
      <xsl:with-param name="self" select="'inputs'"/>
    </xsl:call-template>
    <div class="box">
      <a name="databases"/>
      <h4>Sequence Databases <a href="#databases_doc" class="help"><div class="help"/></a></h4>
      <div class="pad">
        <p>
          <xsl:text>The following sequence database</xsl:text>
          <xsl:choose>
            <xsl:when test="count(/mast/sequences/database) &gt; 1">
              <xsl:text>s were</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text> was</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
          <xsl:text> supplied to MAST.</xsl:text>
        </p>
        <table class="padded-table" border="0" >
          <col />
          <col />
          <col />
          <thead>
            <tr>
              <th style="text-align:left;" >Database</th>
              <th>Sequence Count</th>
              <th>Residue Count</th>
              <th>Last Modified</th>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="/mast/sequences/database">
              <tr>
                <td><xsl:value-of select="@name"/></td>
                <td class="tnum"><xsl:value-of select="@seq_count"/></td>
                <td class="tnum"><xsl:value-of select="@residue_count"/></td>
                <td><xsl:value-of select="@last_mod_date"/></td>
              </tr>
            </xsl:for-each>
          </tbody>
          <tfoot>
            <tr>
              <th style="text-align:left; padding:5px 10px;">Total</th>
              <td class="tnum"><xsl:value-of select="sum(/mast/sequences/database/@seq_count)"/></td>
              <td class="tnum"><xsl:value-of select="sum(/mast/sequences/database/@residue_count)"/></td>
              <td>&nbsp;</td>
            </tr>
          </tfoot>
        </table>
      </div>
      <a name="motifs" />
      <h4>Motifs <a href="#motifs_doc" class="help"><div class="help"/></a></h4>
      <div class="pad">
        <p>The following motifs were supplied to MAST from &quot;<xsl:value-of select="/mast/motifs/@name" />&quot; last modified on <xsl:value-of 
            select="/mast/motifs/@last_mod_date"/>.
          <xsl:if test="count(/mast/motifs/motif[@bad = 'y']) &gt; 0">
            <xsl:choose>
              <xsl:when test="/mast/model/remove_correlated/@value = 'y'">
                Motifs with a red name have been removed because they have a similarity greater 
                than <xsl:value-of select="format-number(/mast/model/max_correlation, '0.00')" /> with another motif.
              </xsl:when>
              <xsl:otherwise>
                It is recommended that motifs with a red name be removed because they have a similarity greater 
                than <xsl:value-of select="format-number(/mast/model/max_correlation, '0.00')" /> with another motif.
              </xsl:otherwise>
            </xsl:choose>
          </xsl:if>
        </p>
        <table class="padded-table" border="0" >
          <col />
          <xsl:choose>
            <xsl:when test="/mast/alphabet/@type = 'nucleotide'">
              <col />
              <col />
            </xsl:when>
            <xsl:otherwise>
              <col />
            </xsl:otherwise>
          </xsl:choose>
          <col />
          <xsl:for-each select="/mast/motifs/motif">
            <col />
          </xsl:for-each>
          <thead>
            <tr>
              <th>&nbsp;</th>
              <th>&nbsp;</th>
              <xsl:choose>
                <xsl:when test="/mast/alphabet/@type = 'nucleotide'">
                  <th colspan="2">Best possible match</th>
                </xsl:when>
                <xsl:otherwise>
                  <th>&nbsp;</th>
                </xsl:otherwise>
              </xsl:choose>
              <th colspan="{count(/mast/motifs/motif)}">Similarity</th>
            </tr>
            <tr>
              <th>Motif</th>
              <th>Width</th>
              <xsl:choose>
                <xsl:when test="/mast/alphabet/@type = 'nucleotide'">
                  <th>(+)</th>
                  <th>(-)</th>
                </xsl:when>
                <xsl:otherwise>
                  <th>Best possible match</th>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:for-each select="/mast/motifs/motif">
                <th>
                  <xsl:if test="number(@name) = NaN">
                    <xsl:text>Motif </xsl:text>
                  </xsl:if>
                  <xsl:value-of select="@name" />
                </th>
              </xsl:for-each>
            </tr>
          </thead>
          <tbody>
            <xsl:for-each select="/mast/motifs/motif">
              <xsl:variable name="motif_a_num" select="position()" />
              <xsl:variable name="motif_a" select="@id" />
              <tr>
                <xsl:choose>
                  <xsl:when test="@bad = 'y'">
                    <td style="color:red;"><xsl:value-of select="@name" /></td>
                  </xsl:when>
                  <xsl:otherwise>
                    <td><xsl:value-of select="@name" /></td>
                  </xsl:otherwise>
                </xsl:choose>
                <td><xsl:value-of select="@width" /></td>
                <td class="sequence"><xsl:call-template name="colour-sequence" ><xsl:with-param name="seq" select="@best_f" /></xsl:call-template></td>
                <xsl:if test="/mast/alphabet/@type = 'nucleotide'">
                  <td class="sequence"><xsl:call-template name="colour-sequence" ><xsl:with-param name="seq" select="@best_r" /></xsl:call-template></td>
                </xsl:if>
                <xsl:for-each select="/mast/motifs/motif">
                  <xsl:variable name="motif_b_num" select="position()" />
                  <xsl:variable name="motif_b" select="@id" />
                  <xsl:choose>
                    <xsl:when test="$motif_a != $motif_b">
                      <xsl:variable name="correlation">
                        <xsl:choose>
                          <xsl:when test="/mast/motifs/correlation[@motif_a = $motif_a and @motif_b = $motif_b]">
                            <xsl:value-of select="/mast/motifs/correlation[@motif_a = $motif_a and @motif_b = $motif_b]/@value"/>
                          </xsl:when>
                          <xsl:otherwise>
                            <xsl:value-of select="/mast/motifs/correlation[@motif_a = $motif_b and @motif_b = $motif_a]/@value"/>
                          </xsl:otherwise>
                        </xsl:choose>
                      </xsl:variable>
                      <xsl:choose>
                        <xsl:when test="$correlation &gt;= /mast/model/max_correlation">
                          <td style="color:red;"><xsl:value-of select="$correlation" /></td>
                        </xsl:when>
                        <xsl:otherwise>
                          <td><xsl:value-of select="$correlation" /></td>
                        </xsl:otherwise>
                      </xsl:choose>
                    </xsl:when>
                    <xsl:otherwise>
                      <td>-</td>
                    </xsl:otherwise>
                  </xsl:choose>
                </xsl:for-each>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>
      </div>
      <xsl:if test="/mast/motifs/nos" >
        <xsl:variable name="base_size">
          <xsl:choose>
            <xsl:when test="/mast/model/translate_dna/@value = 'y'">3</xsl:when>
            <xsl:otherwise>1</xsl:otherwise>
          </xsl:choose>
        </xsl:variable>
        <a name="nos"/>
        <h4>Nominal Order and Spacing <a href="#nos_doc" class="help"><div class="help"/></a></h4>
        <div class="pad">
          <p>The expected order and spacing of the motifs.</p>
          <table style="width:100%;">
            <tr>
              <td class="block_td">
                <div class="block_container">
                  <div class="block_rule" style="width:{((/mast/motifs/nos/@length * $base_size) div $max_seq_len) * 100}%"></div>
                  <xsl:for-each select="/mast/motifs/nos/expect">
                    <xsl:variable name="motif" select="id(@motif)" />
                    <xsl:variable name="motif_num" select="count($motif/preceding-sibling::*) + 1" />

                    <xsl:call-template name="site">
                      <xsl:with-param name="max_seq_len" select="$max_seq_len"/>
                      <xsl:with-param name="max_log_pvalue" select="$max_log_pvalue"/>
                      <xsl:with-param name="position" select="@pos * $base_size"/>
                      <xsl:with-param name="width" select="$motif/@width"/>
                      <xsl:with-param name="index" select="$motif_num"/>
                      <xsl:with-param name="name" select="$motif/@name"/>
                    </xsl:call-template>
                  </xsl:for-each>
                </div>
              </td>
            </tr>
            <tr>
              <td class="block_td" style="color: blue;">
                <div class="block_container" >
                  <xsl:call-template name="ruler">
                    <xsl:with-param name="max" select="/mast/motifs/nos/@length" />
                  </xsl:call-template>
                </div>
              </td>
            </tr>
          </table>
        </div>
      </xsl:if>
    </div>
  </xsl:template>

  <xsl:template name="results">

    <xsl:call-template name="header">
      <xsl:with-param name="title" select="'Search Results'"/>
      <xsl:with-param name="self" select="'results'"/>
    </xsl:call-template>
    <div class="box">
      <h4>Top Scoring Sequences <a href="#sequences_doc" class="help"><div class="help"/></a></h4>
      <div class="pad">
        <p>
          Each of the following <xsl:value-of select="count(/mast/sequences/sequence)"/> sequences has an <i>E</i>-value less than 
          <xsl:value-of select="/mast/model/max_seq_evalue"/>.<br />
          The motif matches shown have a position p-value less than <xsl:value-of select="/mast/model/max_weak_pvalue"/>.<br />
          <xsl:if test="count(/mast/sequences/sequence[@known_pos = 'y'])" >
          The highlighted rows have had correct positions provided for them.<br />
          </xsl:if>
          <b>Click on the arrow</b> (&more;) next to the <i>E</i>-value to view more information about a sequence.
        </p>
        <xsl:call-template name="mast_legend"/>
        <table id="tbl_sequences" style="width:100%; table-layout:fixed;" border="0">
          <col style="width:{$longest_seq_name*0.8}em;" />
          <xsl:choose>
            <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
              <col style="width:5em;" />
              <col style="width:1em;" />
              <col style="width:5em;" />
            </xsl:when>
            <xsl:otherwise>
              <col style="width:5em;" />
            </xsl:otherwise>
          </xsl:choose>
          <col style="width:3em;" />
          <col />
          <thead>
            <th>Sequence</th>
          <xsl:choose>
            <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
              <th colspan="3"><i>E</i>-value</th>
            </xsl:when>
            <xsl:otherwise>
              <th><i>E</i>-value</th>
            </xsl:otherwise>
          </xsl:choose>
            <th>&nbsp;</th>
            <th>Block Diagram</th>
          </thead>
          <tfoot>
            <tr>
              <xsl:choose>
                <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
                  <td colspan="5">&nbsp;</td>
                </xsl:when>
                <xsl:otherwise>
                  <td colspan="3">&nbsp;</td>
                </xsl:otherwise>
              </xsl:choose>
              <td class="block_td" style="color: blue;">
                <div class="block_container" >
                  <xsl:call-template name="ruler">
                    <xsl:with-param name="max" select="$max_seq_len" />
                  </xsl:call-template>
                </div>
              </td>
            </tr>
          </tfoot>
          <tbody>
            <xsl:for-each select="/mast/sequences/sequence" >
              <xsl:variable name="db" select="id(@db)"/>
              <xsl:variable name="seq_is_dna" select="$db/@type = 'nucleotide'"/>
              <xsl:variable name="mul_width">
                <xsl:choose>
                  <xsl:when test="/mast/alphabet/@type = 'amino-acid' and $seq_is_dna">3</xsl:when>
                  <xsl:otherwise>1</xsl:otherwise>
                </xsl:choose>
              </xsl:variable>
              <xsl:variable name="seq_fract" select="(@length div $max_seq_len) * 100" />
              <xsl:variable name="highlight">
                <xsl:if test="@known_pos = 'y'">
                  <xsl:text>highlight</xsl:text>
                </xsl:if>
              </xsl:variable>
              <tr id="{@id}_blocks" class="{$highlight}">
                <td>
                <xsl:choose>
                  <xsl:when test="$db/@link">
                    <xsl:variable name="link">
                      <xsl:value-of select="substring-before($db/@link, 'SEQUENCEID')" />
                      <xsl:value-of select="@name" />
                      <xsl:value-of select="substring-after($db/@link, 'SEQUENCEID')" />
                    </xsl:variable>
                    <a href="{$link}"><xsl:value-of select="@name" /></a>
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:value-of select="@name" />
                  </xsl:otherwise>
                </xsl:choose>
                </td>
                <xsl:choose>
                  <xsl:when test="/mast/model/strand_handling/@value = 'norc'">
                    <td><xsl:value-of select="./score[@strand = 'forward']/@evalue" /></td>
                  </xsl:when>
                  <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
                    <xsl:variable name="f_eval" select="./score[@strand = 'forward']/@evalue"/>
                    <xsl:variable name="r_eval" select="./score[@strand = 'reverse']/@evalue"/>
                    <xsl:variable name="f_dim"><xsl:if test="not($f_eval) or $r_eval and $f_eval &gt; $r_eval">dim</xsl:if></xsl:variable>
                    <xsl:variable name="r_dim"><xsl:if test="not($r_eval) or $f_eval and $r_eval &gt; $f_eval">dim</xsl:if></xsl:variable>
                    <td class="{$f_dim}">
                      <xsl:choose>
                        <xsl:when test="$f_eval" >
                          <xsl:value-of select="$f_eval" />
                        </xsl:when>
                        <xsl:otherwise><xsl:text>--</xsl:text></xsl:otherwise>
                      </xsl:choose>
                    </td>
                    <td>/</td>
                    <td class="{$r_dim}">
                      <xsl:choose>
                        <xsl:when test="$r_eval" >
                          <xsl:value-of select="$r_eval" />
                        </xsl:when>
                        <xsl:otherwise><xsl:text>--</xsl:text></xsl:otherwise>
                      </xsl:choose>
                    </td>
                  </xsl:when>
                  <xsl:otherwise>
                    <td><xsl:value-of select="./score/@evalue" /></td>
                  </xsl:otherwise>
                </xsl:choose>
                <td><a href="javascript:show_more('{@id}')" class="more_arrow" title="Toggle additional information">&more;</a></td>
                <td class="block_td">
                  <div class="block_container" id="{@id}_block_container" >
                    <xsl:if test="$seq_is_dna">
                      <div class="block_plus_sym" >+</div>
                      <div class="block_minus_sym" >-</div>
                    </xsl:if>
                    <div class="block_rule" style="width:{$seq_fract}%"></div>
                    <xsl:for-each select="seg/hit">
                      <xsl:variable name="motif" select="id(@motif)" />
                      <xsl:variable name="motif_num" select="count($motif/preceding-sibling::*) + 1" />
                      <xsl:variable name="strand">
                        <xsl:choose>
                          <xsl:when test="@strand = 'forward'"><xsl:text>plus</xsl:text></xsl:when>
                          <xsl:when test="@strand='reverse'"><xsl:text>minus</xsl:text></xsl:when>
                          <xsl:otherwise>
                          <xsl:value-of select="@strand" />
                          </xsl:otherwise>
                        </xsl:choose>
                      </xsl:variable>
                      <xsl:variable name="frame">
                        <xsl:if test="/mast/model/translate_dna/@value = 'y'">
                          <xsl:variable name="frame_index" select="(@pos - 1) mod 3"/>
                          <xsl:choose>
                            <xsl:when test="$frame_index = 0">
                              <xsl:text>a</xsl:text>
                            </xsl:when>
                            <xsl:when test="$frame_index = 1">
                              <xsl:text>b</xsl:text>
                            </xsl:when>
                            <xsl:otherwise>
                              <xsl:text>c</xsl:text>
                            </xsl:otherwise>
                          </xsl:choose>
                        </xsl:if>
                      </xsl:variable>
                      <xsl:call-template name="site">
                        <xsl:with-param name="max_seq_len" select="$max_seq_len"/>
                        <xsl:with-param name="max_log_pvalue" select="$max_log_pvalue"/>
                        <xsl:with-param name="position" select="@pos"/>
                        <xsl:with-param name="width" select="$motif/@width * $mul_width"/>
                        <xsl:with-param name="index" select="$motif_num"/>
                        <xsl:with-param name="strand" select="$strand"/>
                        <xsl:with-param name="pvalue" select="@pvalue"/>
                        <xsl:with-param name="name" select="$motif/@name"/>
                        <xsl:with-param name="frame" select="$frame"/>
                      </xsl:call-template>
                    </xsl:for-each>
                  </div>
                </td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>
        <xsl:call-template name="mast_legend"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template name="mast_legend">
    <div style="text-align:right">
      <xsl:for-each select="/mast/motifs/motif">
        <xsl:call-template name="legend_motif">
          <xsl:with-param name="name" select="@name"/>
          <xsl:with-param name="index" select="position()"/>
        </xsl:call-template>
      </xsl:for-each>
    </div>
  </xsl:template>

  <xsl:template name="program">
    <a name="program"/>
    <div class="bar">
      <div style="text-align:right;"><a href="#top">Top</a></div>
      <div class="subsection">
        <a name="version"/>
        <h5>MAST version</h5>
        <xsl:value-of select="/mast/@version"/> (Release date: <xsl:value-of select="/mast/@release"/>)
      </div>
      <div class="subsection">
        <a name="reference"/>
        <h5>Reference</h5>
        Timothy L. Bailey and Michael Gribskov,<br />
        "Combining evidence using p-values: application to sequence homology searches", Bioinformatics, 14(48-54), 1998.<br />
      </div>
      <div class="subsection">
        <a name="command" />
        <h5>Command line summary</h5>
        <textarea rows="1" style="width:100%;" readonly="readonly">
          <xsl:value-of select="/mast/model/command_line"/>
        </textarea>
        <br />
        <xsl:text>Background letter frequencies (from </xsl:text>
        <xsl:choose>
          <xsl:when test="/mast/alphabet/@bg_source = 'preset'">
            <xsl:text>non-redundant database</xsl:text>
          </xsl:when>
          <xsl:when test="/mast/alphabet/@bg_source = 'sequence_composition'">
            <xsl:text>sequence composition</xsl:text>
          </xsl:when>
          <xsl:when test="/mast/alphabet/@bg_source = 'file'">
            <xsl:value-of select="/mast/alphabet/@bg_file"/>
          </xsl:when>
        </xsl:choose>
        <xsl:text>):</xsl:text><br/>
        <div style="margin-left:25px;">
          <xsl:for-each select="/mast/alphabet/letter[not(@ambig) or @ambig = 'n']">
            <xsl:value-of select="@symbol"/><xsl:text>:&nbsp;</xsl:text><xsl:value-of select="format-number(@bg_value, '0.000')" />
            <xsl:if test="position() != last()"><xsl:text>&nbsp;&nbsp; </xsl:text></xsl:if>            
          </xsl:for-each>
        </div><br />
        <xsl:text>Result calculation took </xsl:text><xsl:value-of select="/mast/runtime/@seconds"/><xsl:text> seconds</xsl:text>
        <br />
      </div>      
      <a href="javascript:showHidden('model')" id="model_activator">show model parameters...</a>
      <div class="subsection" id="model_data" style="display:none;">
        <h5>Model parameters</h5>
        <xsl:text>&newline;</xsl:text>
        <textarea style="width:100%;" rows="{count(/mast/model/*) - 1}" readonly="readonly">
          <xsl:variable name="spaces" select="'                    '"/>
          <xsl:text>&newline;</xsl:text>
          <xsl:for-each select="/mast/model/*[name(.) != 'command_line']">
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
      <xsl:with-param name="title" select="'Explanation of MAST Results'"/>
      <xsl:with-param name="self" select="'doc'"/>
    </xsl:call-template>
      <div class="box">
        <h4>The MAST results consist of</h4>
        <ul>
          <li>The <a href="#input_motifs_doc"><b>inputs</b></a> to MAST including:
            <ol>
              <li>The <a href="#databases_doc"><b>sequence databases</b></a> showing the sequence 
                and residue counts. [<a href="#databases">View</a>]</li>
              <li>The <a href="#motifs_doc"><b>motifs</b></a> showing the name, width, best scoring match 
                and similarity to other motifs. [<a href="#motifs">View</a>]</li> 
              <li>
                The <a href="#nos_doc"><b>nominal order and spacing</b></a> diagram.<xsl:if test="/mast/motifs/nos"> [<a href="#nos">View</a>]</xsl:if>
              </li>
            </ol>
          </li>
          <li>The <a href="#sequences_doc"><b>search results</b></a> showing top scoring sequences with 
            tiling of all of the motifs matches shown for each of the sequences. [<a href="#results">View</a>]
          </li>
          <li>The <b>program</b> details including:
            <ol>
              <li>The <b>version</b> of MAST and the date it was released. [<a href="#version">View</a>]</li>
              <li>The <b>reference</b> to cite if you use MAST in your research. [<a href="#reference">View</a>]</li>
              <li>The <b>command line summary</b> detailing the parameters with which you ran MAST. [<a href="#command">View</a>]</li>
            </ol>
          </li>
          <li>This <b>explanation</b> of how to interpret MAST results.</li>
        </ul>
        <a name="input_motifs_doc"/>
        <h4>Inputs</h4>
        <p>MAST received the following inputs.</p>
        <a name="databases_doc"/>
        <h5>Sequence Databases</h5>
        <div class="doc">
          <p>This table summarises the sequence databases specified to MAST.</p>
          <dl>
            <dt>Database</dt>
            <dd>The name of the database file.</dd>
            <dt>Sequence Count</dt>
            <dd>The number of sequences in the database.</dd>
            <dt>Residue Count</dt>
            <dd>The number of residues in the database.</dd>
          </dl>
        </div>
        <a name="motifs_doc"/>
        <h5>Motifs</h5>
        <div class="doc">
          <p>Summary of the motifs specified to MAST.</p>
          <dl>
            <dt>Name</dt>
            <dd>The name of the motif. If the motif has been removed or removal is recommended to avoid highly similar motifs 
              then it will be displayed in red text.</dd>
            <dt>Width</dt>
            <dd>The width of the motif. No gaps are allowed in motifs supplied to MAST as it only works for motifs of a fixed width.</dd>
            <dt>Best possible match</dt>
            <dd>The sequence that would achieve the best possible match score and its reverse complement for nucleotide motifs.</dd>
            <dt>Similarity</dt>
            <dd>
              MAST computes the pairwise correlations between each pair of motifs. The correlation between two motifs is the 
              maximum sum of Pearson's correlation coefficients for aligned columns divided by the width of the shorter motif. 
              The maximum is found by trying all alignments of the two motifs. Motifs with correlations below 0.60 have little 
              effect on the accuracy of the combined scores. Pairs of motifs with higher correlations should be removed from 
              the query. Correlations above the supplied threshold are shown in red text.
            </dd>
          </dl>
        </div>
        <a name="nos_doc"/>
        <h5>Nominal Order and Spacing</h5>
        <div class="doc">
          <p>This diagram shows the normal spacing of the motifs specified to MAST.</p>
        </div>
        <a name="sequences_doc"/>
        <h4>Search Results</h4>
        <p>MAST provides the following motif search results.</p>
        <h5>Top Scoring Sequences</h5>
        <div class="doc">
          <p>
            This table summarises the top scoring sequences with a <a href="#evalue_doc">Sequence <i>E</i>-value</a>
            better than the threshold (default 10). The sequences are sorted by the <a href="#evalue_doc">Sequence 
            <i>E</i>-value</a> from most to least significant.
          </p>
          <dl>
            <dt>Sequence</dt>
            <dd>The name of the sequence. This maybe be linked to search a sequence database for the sequence name.</dd>
            <dt><i>E</i>-value</dt>
            <dd>
              The <a href="#evalue_doc"><i>E</i>-value</a> of the sequence. For DNA only; if strands were scored seperately 
              then there will be two <i>E</i>-values for the sequence seperated by a "/". The score for the provided sequence 
              will be first and the score for the reverse-complement will be second.
            </dd>
            <dt>&more;</dt>
            <dd>
              Click on this to show <a href="#additional_seq_doc">additional information</a> about the sequence such as a 
              description, combined p-value and the annotated sequence.
            </dd>
            <dt>Block Diagram</dt>
            <dd>
              The block diagram shows the best non-overlapping tiling of motif matches on the sequence.
              <ul style="padding-left: 1em; margin-top:5px;">
                <li>The length of the line shows the length of a sequence relative to all the other sequences.</li>
                <li>A block is shown where the <a href="#pos_pvalue_doc" >positional <i>p</i>-value</a>
                  of a motif is less (more significant) than the significance threshold which is 0.0001 by default.</li>
                <li>If a significant motif match (as specified above) overlaps other significant motif matches then 
                  it is only displayed as a block if its <a href="#pos_pvalue_doc" >positional <i>p</i>-value</a> 
                  is less (more significant) then the product of the <a href="#pos_pvalue_doc" >positional 
                  <i>p</i>-values</a> of the significant matches that it overlaps.</li>
                <li>The position of a block shows where a motif has matched the sequence.</li>
                <li>The width of a block shows the width of the motif relative to the length of the sequence.</li>
                <li>The colour and border of a block identifies the matching motif as in the legend.</li>
                <li>The height of a block gives an indication of the significance of the match as 
                  taller blocks are more significant. The height is calculated to be proportional 
                  to the negative logarithm of the <a href="#pos_pvalue_doc" >positional <i>p</i>-value</a>, 
                  truncated at the height for a <i>p</i>-value of 1e-10.</li>
                  <li>Hovering the mouse cursor over the block causes the display of the motif name 
                    and other details in the hovering text.</li>
                <li>DNA only; blocks displayed above the line are a match on the given DNA, whereas blocks 
                  displayed below the line are matches to the reverse-complement of the given DNA.</li>
                <li>DNA only; when strands are scored separately then blocks may overlap on opposing strands.</li>
              </ul>
            </dd>
          </dl>
        </div>
        <a name="additional_seq_doc"/>
        <h5>Additional Sequence Information</h5>
        <div class="doc">
          <p>
            Clicking on the &more; link expands a box below the sequence with additional information and adds two dragable buttons 
            below the block diagram.
          </p>
          <dl>
            <dt>Description</dt>
            <dd>The description appearing after the identifier in the fasta file used to specify the sequence.</dd>
            <dt>Combined <i>p</i>-value</dt>
            <dd>The <a href="#combined_pvalue_doc">combined <i>p</i>-value</a> of the sequence. DNA only; if strands were scored 
              seperately then there will be two <i>p</i>-values for the sequence seperated by a "/". The score for the provided sequence 
              will be first and the score for the reverse-complement will be second.</dd>
            <a name="anno_doc"/>
            <dt>Annotated Sequence</dt>
            <dd>
              The annotated sequence shows a portion of the sequence with the matching motif sequences displayed above. 
              The displayed portion of the sequence can be modified by sliding the two buttons below the sequence block diagram
              so that the portion you want to see is between the two needles attached to the buttons. By default the two buttons
              move together but you can drag one individually by holding shift before you start the drag. If the strands were
              scored seperately then they can't be both displayed at once due to overlaps and so a radio button offers the choice
              of strand to display.
            </dd>
          </dl>
        </div>
        <h4>Scoring</h4>
        <p>MAST scores sequences using the following measures.</p>
        <a name="score_doc" />
        <h5>Position score calculation</h5>
        <div class="doc">
          <p>
            The score for the match of a position in a sequence to a motif is computed by by summing the appropriate entry 
            from each column of the position-dependent scoring matrix that represents the motif. Sequences shorter than 
            one or more of the motifs are skipped.
          </p>
        </div>
        <a name="pos_pvalue_doc" />
        <h5>Position <i>p</i>-value</h5>
        <div class="doc">
          <p>
            The position p-value of a match is the probability of a single random subsequence of the length of the motif 
            scoring at least as well as the observed match.
          </p>
        </div>
        <a name="seq_pvalue_doc" />
        <h5>Sequence <i>p</i>-value</h5>
        <div class="doc">
          <p>
            The sequence p-value of a score is defined as the probability of a random sequence of the same length containing 
            some match with as good or better a score.
          </p>
        </div>
        <a name="combined_pvalue_doc" />
        <h5>Combined <i>p</i>-value</h5>
        <div class="doc">
          <p>The combined p-value of a sequence measures the strength of the match of the sequence to all the motifs and is calculated by</p>
          <ol>
            <li>finding the score of the single best match of each motif to the sequence (best matches may overlap),</li>
            <li>calculating the sequence p-value of each score,</li>
            <li>forming the product of the p-values,</li>
            <li>taking the p-value of the product.</li>
          </ol>
        </div>
        <a name="evalue_doc" />
        <h5>Sequence <i>E</i>-value</h5>
        <div class="doc">
          <p>
            The E-value of a sequence is the expected number of sequences in a random database of the same size that would match 
            the motifs as well as the sequence does and is equal to the combined p-value of the sequence times the number of 
            sequences in the database.
          </p>
        </div>
      </div>
    </span>
  </xsl:template>

  <xsl:template name="data" >
    <xsl:text>&newline;</xsl:text>
    <xsl:comment >This data is used by the javascript to create the detailed views</xsl:comment><xsl:text>&newline;</xsl:text>
    <form>
      <xsl:for-each select="/mast/sequences/sequence">
        <xsl:variable name="segs">
          <xsl:text>&newline;</xsl:text>
          <xsl:for-each select="seg"> 
            <xsl:value-of select="@start"/><xsl:text>&tab;</xsl:text>
            <xsl:value-of select="translate(data, ' &#10;&#9;', '')" />
            <xsl:text>&newline;</xsl:text>
          </xsl:for-each>
        </xsl:variable>
        <xsl:variable name="hits">
          <xsl:text>&newline;</xsl:text>
          <xsl:for-each select=".//hit"> 
            <xsl:value-of select="@pos" /><xsl:text>&tab;</xsl:text> 
            <xsl:value-of select="@motif" /><xsl:text>&tab;</xsl:text>
            <xsl:choose>
              <xsl:when test="not (@strand) or @strand = 'forward'">
                <xsl:text>+</xsl:text>
              </xsl:when>
              <xsl:otherwise>
                <xsl:text>-</xsl:text>
              </xsl:otherwise>
            </xsl:choose><xsl:text>&tab;</xsl:text> 
            <xsl:value-of select="@pvalue"/><xsl:text>&tab;</xsl:text> 
            <xsl:value-of select="@match"/><xsl:text>&tab;</xsl:text>
            <xsl:value-of select="@translation"/><xsl:text>&newline;</xsl:text> 
          </xsl:for-each>
        </xsl:variable>
        <xsl:variable name="combined_pvalue">
          <xsl:choose>
            <xsl:when test="/mast/model/strand_handling/@value = 'norc'">
              <xsl:value-of select="./score[@strand = 'forward']/@combined_pvalue" />
            </xsl:when>
            <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
              <xsl:choose>
                <xsl:when test="./score[@strand = 'forward']/@combined_pvalue">
                  <xsl:value-of select="./score[@strand = 'forward']/@combined_pvalue" />
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>--&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</xsl:text>
                </xsl:otherwise>
              </xsl:choose>
              <xsl:text>&nbsp;&nbsp;/&nbsp;&nbsp;</xsl:text>
              <xsl:choose>
                <xsl:when test="./score[@strand = 'reverse']/@combined_pvalue">
                  <xsl:value-of select="./score[@strand = 'reverse']/@combined_pvalue" />
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>--</xsl:text>
                </xsl:otherwise>
              </xsl:choose>
            </xsl:when>
            <xsl:otherwise>
              <xsl:value-of select="./score/@combined_pvalue" />
            </xsl:otherwise>
          </xsl:choose>
        </xsl:variable>
        <xsl:if test="/mast/model/translate_dna/@value = 'y'">
          <xsl:variable name="best_frame">
            <xsl:choose>
              <xsl:when test="/mast/model/strand_handling/@value = 'norc'">
                <xsl:value-of select="./score[@strand = 'forward']/@frame" />
              </xsl:when>
              <xsl:when test="/mast/model/strand_handling/@value = 'separate'">
                <xsl:choose>
                  <xsl:when test="./score[@strand = 'forward']/@frame">
                    <xsl:value-of select="./score[@strand = 'forward']/@frame" />
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:text>--&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</xsl:text>
                  </xsl:otherwise>
                </xsl:choose>
                <xsl:text>&nbsp;&nbsp;/&nbsp;&nbsp;</xsl:text>
                <xsl:choose>
                  <xsl:when test="./score[@strand = 'reverse']/@frame">
                    <xsl:value-of select="./score[@strand = 'reverse']/@frame" />
                  </xsl:when>
                  <xsl:otherwise>
                    <xsl:text>--</xsl:text>
                  </xsl:otherwise>
                </xsl:choose>
              </xsl:when>
              <xsl:otherwise>
                <xsl:value-of select="./score/@frame" />
              </xsl:otherwise>
            </xsl:choose>
          </xsl:variable>
          <input type="hidden" id="{@id}_frame" value="{$best_frame}" /><xsl:text>&newline;</xsl:text>
        </xsl:if>
        <input type="hidden" id="{@id}_len" value="{@length}" /><xsl:text>&newline;</xsl:text>
        <input type="hidden" id="{@id}_desc" value="{@comment}" /><xsl:text>&newline;</xsl:text>
        <input type="hidden" id="{@id}_combined_pvalue" value="{$combined_pvalue}" /><xsl:text>&newline;</xsl:text>
        <input type="hidden" id="{@id}_type" value="{id(@db)/@type}" /><xsl:text>&newline;</xsl:text>
        <input type="hidden" id="{@id}_segs" value="{$segs}" /><xsl:text>&newline;</xsl:text>
        <input type="hidden" id="{@id}_hits" value="{$hits}" /><xsl:text>&newline;</xsl:text>
      </xsl:for-each>
    </form>
  </xsl:template>

  <xsl:template name="colour-sequence">
    <xsl:param name="seq"/>
    <xsl:variable name="colour">
      <xsl:call-template name="pick-colour">
        <xsl:with-param name="sym" select="substring($seq, 1, 1)"/>
      </xsl:call-template>
    </xsl:variable>
    <span style="color:{$colour}"><xsl:value-of select="substring($seq, 1, 1)"/></span>
    <xsl:if test="string-length($seq) &gt; 1">
      <xsl:call-template name="colour-sequence">
        <xsl:with-param name="seq" select="substring($seq, 2)"/>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>

  <xsl:template name="pick-colour">
    <xsl:param name="sym"/>
    <xsl:param name="isnuc" select="/mast/alphabet/@type = 'nucleotide'" />
    <xsl:choose>
      <xsl:when test="$isnuc">
        <xsl:choose>
          <xsl:when test="$sym = 'A'">
            <xsl:text>red</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'C'">
            <xsl:text>blue</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'G'">
            <xsl:text>orange</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'T'">
            <xsl:text>green</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>black</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <xsl:choose>
          <xsl:when test="contains('ACFILVWM', $sym)">
            <xsl:text>blue</xsl:text>
          </xsl:when>
          <xsl:when test="contains('NQST', $sym)">
            <xsl:text>green</xsl:text>
          </xsl:when>
          <xsl:when test="contains('DE', $sym)">
            <xsl:text>magenta</xsl:text>
          </xsl:when>
          <xsl:when test="contains('KR', $sym)">
            <xsl:text>red</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'H'">
            <xsl:text>pink</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'G'">
            <xsl:text>orange</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'P'">
            <xsl:text>yellow</xsl:text>
          </xsl:when>
          <xsl:when test="$sym = 'Y'">
            <xsl:text>turquoise</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>black</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:otherwise>
    </xsl:choose>
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
</xsl:stylesheet>
