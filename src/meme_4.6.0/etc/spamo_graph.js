      var SpamoEnum = {"INTEGER" : 1, "DECIMAL" : 2, "STRING" : 3, "CMD" : 4, "COMMENT" : 5, "EOF" : 6, "ERROR" : 7};

      var SpamoConstants = {
        "size_title" : 18,
        "font_title" : "normal 18px Times, sans-serif",
        "size_subtitle" : 14,
        "font_subtitle" : "normal 14px Times, sans-serif",
        //font for the axis descriptions
        "size_axis" : 12,
        "font_axis" : "normal 12px Times, sans-serif",
        //font for the numbers on the axis
        "size_label" : 8,
        "font_label" : "normal 8px Times, sans-serif",
        "size_tic" : 2,
        "spacer_tic" : 2,
        "spacer_side" : 5,
        "spacer_top" : 5,
        "spacer_bottom" : 5,
        "spacer_title" : 0,
        "spacer_subtitle" : 5,
        "spacer_axis" : 0,
        "spacer_base" : 3,
        "spacer_mid" : 10,
        "spacer_border" : 5,
        "y_axis_label_start" : "# of Sequences with Spacing (total = ",
        "y_axis_label_end" : ")",
        "x_axis_label" : "Spacings from Primary to Secondary Motif (Gap)",
        "same_strand_label" : "Same Strand",
        "oppo_strand_label" : "Opposite Strand",
        "upstream_label" : "Upstream",
        "downstream_label" : "Downstream",
      };

      var SpamoMeasureBox;

      /*
       * Measures the character heights in string
       * by painting them over the top of each other and finding the
       * bounding box. Returns the height of the total string
       * as well as calculating the vertical centering offset.
       */
      function string_height(text, font, size) {
        var i, x, y, box;
        x = size / 2;
        y = size + (size / 2);
        if (!SpamoMeasureBox) SpamoMeasureBox = document.createElement("canvas");
        SpamoMeasureBox.width = size * 2;
        SpamoMeasureBox.height = size * 2;
        ctx = SpamoMeasureBox.getContext('2d');
        ctx.font = font;
        for (i = 0; i < text.length; ++i) {
          ctx.fillText(text.charAt(i), x, y);
        }
        box = bbox(ctx, size * 2, size * 2, true);

        return {height: box.height, vcenteroffset: (box.height / 2)}
      }

      /*
       * Calculates a tight bounding box on the contents of a canvas
       */
      function bbox(ctx, cwidth, cheight, skip_width_calculation, skip_height_calculation) {
        if (typeof skip_width_calculation != 'boolean') skip_width_calculation = false;
        if (typeof skip_height_calculation != 'boolean') skip_height_calculation = false;
        if (skip_width_calculation && skip_height_calculation) throw "bbox: can't skip calculation of both width and height";
        var data = ctx.getImageData(0, 0, cwidth, cheight).data;
        var r = 0, c = 0;// r: row, c: column
        var top_line = -1, bottom_line = -1, left_line = -1, right_line = -1;
        var box_width = 0, box_height = 0;
        if (!skip_height_calculation) {
          // Find the top-most line with a non-white pixel
          for (r = 0; r < cheight; r++) {
            for (c = 0; c < cwidth; c++) {
              if (data[r * cwidth * 4 + c * 4 + 3]) {
                top_line = r;
                break;
              }
            }
            if (top_line != -1) break;
          }
          
          //find the last line with a non-white pixel
          if (top_line != -1) {
            for (r = cheight-1; r >= top_line; r--) {
              for(c = 0; c < cwidth; c++) {
                if(data[r * cwidth * 4 + c * 4 + 3]) {
                  bottom_line = r;
                  break;
                }
              }
              if (bottom_line != -1) break;
            }
            box_height = bottom_line - top_line + 1;
          }
        }

        if (!skip_width_calculation) {
          // Find the left-most line with a non-white pixel
          for (c = 0; c < cwidth; c++) {
            for (r = 0; r < cheight; r++) {
              if (data[r * cwidth * 4 + c * 4 + 3]) {
                left_line = c;
                break;
              }
            }
            if (left_line != -1) break;
          }

          //find the right most line with a non-white pixel
          if (left_line != -1) {
            for (c = cwidth-1; c >= left_line; c--) {
              for(r = 0; r < cheight; r++) {
                if(data[r * cwidth * 4 + c * 4 + 3]) {
                  right_line = c;
                  break;
                }
              }
              if (right_line != -1) break;
            }
            box_width = right_line - left_line + 1;
          }
        }
        //return the bounds
        if (!skip_height_calculation && !skip_width_calculation) {
          return {top: top_line, bottom: bottom_line, left: left_line, right: right_line, width: box_width, height: box_height};
        } else if (!skip_height_calculation) {
          return {top: top_line, bottom: bottom_line, height: box_height};
        } else {
          return {left: left_line, right: right_line, width: box_width};
        }
      }


      function SpamoHistogram(text) {
        //initialise graph state to default
        this.graph_title = "";
        this.graph_margin = 150;
        this.graph_mlength = 10;
        this.graph_binsize = 1;
        this.graph_binmax = 50;
        this.graph_sequences = 0;
        this.graph_significant = 0.05;
        this.graph_has_border = false;
        this.graph_has_subtitle = false;
        this.graph_has_y_axis_label = true;
        this.graph_has_x_axis_label = true;
        this.graph_has_stream_label = false;
        this.graph_has_strand_label = false;
        

        var tokeniser = new SpamoTokeniser(text);
        var token = tokeniser.next();
        var stack = [];
        while (token.type != SpamoEnum.EOF) {
          switch (token.type) {
            case SpamoEnum.ERROR:
              die(token.value);
            case SpamoEnum.INTEGER:
            case SpamoEnum.DECIMAL:
            case SpamoEnum.STRING:
              stack.push(token);
              break;
            case SpamoEnum.CMD:
              if (token.value == "title") {
                var title = stack.pop();
                if (title === undefined || title.type != SpamoEnum.STRING) die("Command 'title' expects string on stack.");
                this.graph_title = title.value;
              } else if (token.value == "margin") {
                var margin = stack.pop();
                if (margin === undefined || margin.type != SpamoEnum.INTEGER) die("Command 'margin' expects integer on stack.");
                if (margin.value <= 0) die("Command 'margin' expects an integer 1 or larger.");
                this.graph_margin = margin.value;
              } else if (token.value == "motif-length") {
                var mlen = stack.pop();
                if (mlen === undefined || mlen.type != SpamoEnum.INTEGER) die("Command 'motif-length' expects integer on stack.");
                if (mlen.value <= 0) die("Command 'motif-length' expects an integer 1 or larger.");
                this.graph_mlength = mlen.value;
              } else if (token.value == "bin-size") {
                var bin = stack.pop();
                if (bin === undefined || bin.type != SpamoEnum.INTEGER) die("Command 'bin-size' expects integer on stack.");
                if (bin.value <= 0) die("Command 'bin-size' expects an integer 1 or larger.");
                this.graph_binsize = bin.value;
              } else if (token.value == "bin-max") {
                var bin_max = stack.pop();
                if (bin_max === undefined || bin_max.type != SpamoEnum.INTEGER) die("Command 'bin-max' expects integer on stack.");
                if (bin_max.value < 0) die("Command 'bin-max' expects a positive integer.");
                this.graph_binmax = bin_max.value;
              } else if (token.value == "sequences") {
                var sequences = stack.pop();
                if (sequences === undefined || bin_max.type != SpamoEnum.INTEGER) die("Command 'sequences' expects integer on stack.");
                if (sequences.value < 0) die("Command 'sequences' expects a positive integer.");
                this.graph_sequences = sequences.value;
              } else if (token.value == "threshold") {
                var thresh = stack.pop();
                if (thresh === undefined || (thresh.type != SpamoEnum.INTEGER && thresh.type != SpamoEnum.DECIMAL)) 
                  die("Command 'threshold' expects a number on stack.");
                if (thresh.value < 0 || thresh.value >= 1) die("Command 'threshold' expects 0 <= value < 1.");
                this.graph_significant = thresh.value;
              } else if (token.value == "show-border") {
                this.graph_has_border = true;
              } else if (token.value == "show-subtitle") {
                this.graph_has_subtitle = true;
              } else if (token.value == "hide-y-axis-label") {
                this.graph_has_y_axis_label = false;
              } else if (token.value == "hide-x-axis-label") {
                this.graph_has_x_axis_label = false;
              } else if (token.value == "show-stream-labels") {
                this.graph_has_stream_label = true;
              } else if (token.value == "show-strand-labels") {
                this.graph_has_strand_label = true;
              } else {
                die("Unknown command '"+token.value+"'");
              }
          }
          token = tokeniser.next();
        }
        //calculate bin count to read in the bin values
        var possible_count = this.graph_margin - this.graph_mlength + 1;
        var bin_count = Math.floor(possible_count / this.graph_binsize) + (possible_count % this.graph_binsize > 1 ? 1 : 0);
        var i;
        
        //get the bins
        this.bins_opposite_right = SpamoHistogram_bins(stack, bin_count);
        this.bins_opposite_left = SpamoHistogram_bins(stack, bin_count);
        this.bins_same_right = SpamoHistogram_bins(stack, bin_count);
        this.bins_same_left = SpamoHistogram_bins(stack, bin_count);

        //check the stack
        if (stack.length != 0) die("Everything loaded but stack not empty!");

        //calculate variables
        this.graph_has_title = this.graph_title.length > 0;
        this.graph_subtitle = "margin: " + this.graph_margin + "   bin size: " + this.graph_binsize;
        this.y_axis_max = (Math.floor(this.graph_binmax / 5) + 1) * 5;
        this.full_bin_count = Math.floor((this.graph_margin - this.graph_mlength + 1) / this.graph_binsize);
      }

      function SpamoHistogram_bins(stack, bin_count) {
        bins = new Array(bin_count);
        for (i = bin_count - 1; i >= 0; --i) {
          var pvalue = stack.pop();
          if (pvalue === undefined || (pvalue.type != SpamoEnum.INTEGER && pvalue.type != SpamoEnum.DECIMAL) 
              || pvalue.value < 0 || pvalue.value > 1) die("Expected pvalue on stack.");
          var count = stack.pop();
          if (count === undefined || (count.type != SpamoEnum.INTEGER && count.type != SpamoEnum.DECIMAL)
              || count.value < 0) die("Expected count on stack.");
          bins[i] = new SpamoBin(count.value, pvalue.value); 
        }
        return bins;
      }

      function SpamoBin(count, pvalue) {
        this.count = count;
        this.pvalue = pvalue;
      }

      function SpamoTokeniser(text) {
        this.text = text;
        this.pos = 0;
        this.next = SpamoTokeniser_next;
      }

      function SpamoTokeniser_next() {
        var whitespace = /\s/;
        //skip over whitespace
        for (;this.pos < this.text.length; this.pos++) {
          if (!this.text.charAt(this.pos).match(whitespace)) break;
        }
        //check for eof
        if (this.pos >= this.text.length) return new SpamoToken(SpamoEnum.EOF, null);
        //check for a comment
        if (this.text.charAt(this.pos) == '%') {
          this.pos++;
          //find the next new line
          var eol = this.text.indexOf("\n", this.pos);
          var line;
          if (eol == -1) {
            //last line in file!
            line = this.text.substring(this.pos);
            this.pos = this.text.length;
          } else {
            line = this.text.substring(this.pos, eol);
            this.pos = eol + 1;
          }
          return new SpamoToken(SpamoEnum.COMMENT, line);
        }
        //check for a string
        if (this.text.charAt(this.pos) == '(') {
          var start = ++this.pos;
          var eos;
          while (true) {
            //find the end of the string
            eos = this.text.indexOf(")", this.pos);
            if (eos == -1) {
              //uhoh, unterminated string!
              return new SpamoToken(SpamoEnum.ERROR, "Unterminated string at " + this.pos);
            }
            this.pos = eos+1;
            if (this.text.charAt(eos-1) != "\\") break;
          }
          return new SpamoToken(SpamoEnum.STRING, this.text.substring(start, eos));
        }
        //read the next word
        {
          var start = this.pos;
          for (;this.pos < this.text.length; this.pos++) {
            if (this.text.charAt(this.pos).match(whitespace)) break;
          }
          var word = this.text.substring(start, this.pos);
          //check for integer
          if (word.match(/^\d+$/)) {
            return new SpamoToken(SpamoEnum.INTEGER, parseInt(word, 10));
          }
          //check for decimal
          var num = word * 1;
          if (!isNaN(num)) {
            return new SpamoToken(SpamoEnum.DECIMAL, num);
          }
          //assume must be command
          return new SpamoToken(SpamoEnum.CMD, word);
        } 
      }

      function SpamoToken(type, value) {
        this.type = type;
        this.value = value;
      }

      function SpamoMetrics(spamo, canvas, ctx) {
        var i;

        this.width = canvas.width;
        this.height = canvas.height;

        //count axis numbers width
        this.count_width = 0;
        ctx.font = SpamoConstants.font_label;
        for (i = 0; i <= spamo.y_axis_max; i+= 5) {
          var width = ctx.measureText("" + i).width;
          if (width > this.count_width) {
            this.count_width = width;
          }
        }
        //spacing axis numbers width
        this.spacing_width = 0;
        ctx.font = SpamoConstants.font_label;
        for (i = -spamo.graph_margin; i <= spamo.graph_margin; ++i) {
          var width = ctx.measureText("" + i).width;
          if (width > this.spacing_width) {
            this.spacing_width = width;
          }
        }
        //calculate the height of a quadrant
        this.quad_height = this.height;
        if (spamo.graph_has_border) {
          this.quad_height -= 2 * SpamoConstants.spacer_border;
        }
        this.quad_height -= SpamoConstants.spacer_top;
        if (spamo.graph_has_title) {
          this.quad_height -= SpamoConstants.size_title;
        }
        if (spamo.graph_has_title && spamo.graph_has_subtitle) {
          this.quad_height -= SpamoConstants.spacer_title;
        }
        if (spamo.graph_has_subtitle) {
          this.quad_height -= SpamoConstants.size_subtitle;
        }
        if (spamo.graph_has_title || spamo.graph_has_subtitle) {
          this.quad_height -= SpamoConstants.spacer_subtitle;
        }
        if (spamo.graph_has_stream_label) {
          this.quad_height -= SpamoConstants.size_axis;
        }
        this.quad_height -= SpamoConstants.spacer_base * 2;
        this.quad_height -= this.spacing_width;
        if (spamo.graph_has_x_axis_label) {
          this.quad_height -= SpamoConstants.spacer_axis + SpamoConstants.size_axis;
        }
        this.quad_height -= SpamoConstants.spacer_bottom;
        this.quad_height /= 2;

        //calculate the max width of a quadrant
        this.quad_max_width = this.width;
        if (spamo.graph_has_border) {
          this.quad_max_width -= 2 * SpamoConstants.spacer_border;
        }
        this.quad_max_width -= 2 * SpamoConstants.spacer_side;
        if (spamo.graph_has_y_axis_label) {
          this.quad_max_width -= SpamoConstants.size_axis + SpamoConstants.spacer_axis; 
        }
        this.quad_max_width -= this.count_width + SpamoConstants.spacer_tic + SpamoConstants.spacer_mid;
        if (spamo.graph_has_strand_label) {
          this.quad_max_width -= SpamoConstants.size_axis + SpamoConstants.spacer_axis;
        }
        this.quad_max_width /= 2;

        //calculate the width of a graph bar the size of a single unit
        this.unit_width = this.quad_max_width / spamo.graph_margin;
        if (this.unit_width > 1) this.unit_width = Math.floor(this.unit_width);

        //calculate the width of the first possible bar
        this.partial_bar_width = ((spamo.graph_margin - spamo.graph_mlength + 1) % spamo.graph_binsize) * this.unit_width;

        //calculate the width of a full bar
        this.full_bar_width = this.unit_width * spamo.graph_binsize;

        //calculate the true width of a quadrant
        this.quad_width = this.unit_width * spamo.graph_margin;

        //add the excess onto the central spacer
        this.spacer_mid = SpamoConstants.spacer_mid + (this.quad_max_width - this.quad_width) * 2;

        //calculate the left quad
        this.left_quad = 0;
        if (spamo.graph_has_border) {
          this.left_quad += SpamoConstants.spacer_border;
        }
        this.left_quad += SpamoConstants.spacer_side;
        if (spamo.graph_has_y_axis_label) {
          this.left_quad += SpamoConstants.size_axis + SpamoConstants.spacer_axis;
        }
        this.left_quad += this.count_width + SpamoConstants.spacer_tic + SpamoConstants.size_tic;

        //calculate the right quad
        this.right_quad = this.left_quad + this.quad_width + this.spacer_mid;

        //calculate the bottom quad
        this.bottom_quad = 0;
        if (spamo.graph_has_border) {
          this.bottom_quad += SpamoConstants.spacer_border;
        }
        this.bottom_quad += SpamoConstants.spacer_bottom;
        if (spamo.graph_has_x_axis_label) {
          this.bottom_quad += SpamoConstants.size_axis + SpamoConstants.spacer_axis;
        }

        //calculate the top quad
        this.top_quad = this.bottom_quad + this.quad_height + (2 * SpamoConstants.spacer_base) + this.spacing_width;
        
        
      }

      function calculate_spamo_bar_height(spamo, metrics, counts) {
        return (counts / spamo.y_axis_max) * metrics.quad_height;
      }

      function die(msg) {
        throw msg;
      }

      function draw_spamo_border(spamo, metrics, ctx) {
        var bdr = SpamoConstants.spacer_border;
        var bdr2 = 2 * bdr;
        //add on 0.5 so the rectangle is in the center of a pixel and so is only 1 pixel wide
        ctx.strokeRect(bdr + 0.5, bdr + 0.5, metrics.width - bdr2, metrics.height - bdr2); 
      }

      function draw_spamo_title(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_title;
        var title_width = ctx.measureText(spamo.graph_title).width;
        var x = (metrics.width / 2) - (title_width / 2);
        var y = metrics.height - SpamoConstants.spacer_top - SpamoConstants.size_title;
        if (spamo.graph_has_border) {
          y -= SpamoConstants.spacer_border;
        }
        ctx.fillText(spamo.graph_title, x, metrics.height - y);
      }

      function draw_spamo_subtitle(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_subtitle;
        var subtitle_width = ctx.measureText(spamo.graph_subtitle).width;
        var x = (metrics.width / 2) - (subtitle_width / 2);
        var y = metrics.height - SpamoConstants.spacer_top - SpamoConstants.size_subtitle;
        if (spamo.graph_has_border) {
          y -= SpamoConstants.spacer_border;
        }
        if (spamo.graph_has_title) {
          y -= (SpamoConstants.size_title + SpamoConstants.spacer_title);
        }
        ctx.fillText(spamo.graph_subtitle, x, metrics.height - y);
      }

      function draw_spamo_y_axis_label(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_axis;
        var txt = SpamoConstants.y_axis_label_start + spamo.graph_sequences + SpamoConstants.y_axis_label_end;
        var axislabel_width = ctx.measureText(txt).width;
        var x = SpamoConstants.spacer_side + SpamoConstants.size_axis;
        if (spamo.graph_has_border) {
          x += SpamoConstants.spacer_border;
        }
        var y = (metrics.height / 2) - (axislabel_width / 2);
        ctx.save();
        ctx.translate(x, metrics.height - y);
        ctx.rotate(-(Math.PI / 2));
        ctx.fillText(txt, 0, 0);
        ctx.restore();
      }

      function draw_spamo_x_axis_label(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_axis;
        var axislabel_width = ctx.measureText(SpamoConstants.x_axis_label).width;
        var x = (metrics.width / 2) - (axislabel_width / 2);
        var y = SpamoConstants.spacer_bottom;
        if (spamo.graph_has_border) {
          y += SpamoConstants.spacer_border;
        }
        ctx.fillText(SpamoConstants.x_axis_label, x, metrics.height - y);
      }

      function draw_spamo_strand_labels(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_axis;
        var upstream_width = ctx.measureText(SpamoConstants.upstream_label).width;
        var downstream_width = ctx.measureText(SpamoConstants.downstream_label).width;
        var y = metrics.top_quad + metrics.quad_height;
        var x = metrics.left_quad + (metrics.quad_width / 2) - (upstream_width / 2);
        ctx.fillText(SpamoConstants.upstream_label, x, metrics.height - y);
        x = metrics.right_quad + (metrics.quad_width / 2) - (downstream_width / 2);
        ctx.fillText(SpamoConstants.downstream_label, x, metrics.height - y);
      }

      function draw_spamo_stream_labels(spamo, metrics, ctx) {
        ctx.font = SpamoConstants.font_axis;
        var same_strand_width = ctx.measureText(SpamoConstants.same_strand_label).width;
        var oppo_strand_width = ctx.measureText(SpamoConstants.oppo_strand_label).width;
        var x = metrics.right_quad + metrics.quad_width + SpamoConstants.spacer_axis + SpamoConstants.size_axis;
        var y = metrics.top_quad + (metrics.quad_height / 2) - (same_strand_width / 2);
        ctx.save();
        ctx.translate(x, metrics.height - y);
        ctx.rotate(-(Math.PI / 2));
        ctx.fillText(SpamoConstants.same_strand_label, 0, 0);
        ctx.restore();
        y = metrics.bottom_quad + (metrics.quad_height / 2) - (oppo_strand_width / 2);
        ctx.save();
        ctx.translate(x, metrics.height - y);
        ctx.rotate(-(Math.PI / 2));
        ctx.fillText(SpamoConstants.oppo_strand_label, 0, 0);
        ctx.restore();
      }

      function draw_spamo_y_axis_part(spamo, metrics, ctx, start_y, height) {
        ctx.font = SpamoConstants.font_label;
        ctx.beginPath();
        ctx.moveTo(metrics.left_quad - 0.5, metrics.height - start_y);
        ctx.lineTo(metrics.left_quad - 0.5, metrics.height - (start_y + height));
        ctx.stroke();
        y_inc = -(height / spamo.y_axis_max * 5);
        ctx.save();
        ctx.translate(metrics.left_quad, metrics.height - start_y);
        for (i = 0; i <= spamo.y_axis_max; i += 5) {
          var lbl = "" + i;
          var lbl_width = ctx.measureText(lbl).width;
          var lbl_heightbox = string_height(lbl, SpamoConstants.font_label, SpamoConstants.size_label);
          var x = -(lbl_width + SpamoConstants.size_tic);
          var y = lbl_heightbox.vcenteroffset;
          //draw text
          ctx.fillText(lbl, x, y);

          //draw tic
          ctx.beginPath();
          ctx.moveTo(0, 0.5);
          ctx.lineTo(-SpamoConstants.size_tic, 0.5);
          ctx.stroke();

          ctx.translate(0, y_inc);
        }
        ctx.restore();
      }

      function draw_spamo_y_axis(spamo, metrics, ctx) {
        draw_spamo_y_axis_part(spamo, metrics, ctx, metrics.top_quad, metrics.quad_height);  
        draw_spamo_y_axis_part(spamo, metrics, ctx, metrics.bottom_quad + metrics.quad_height, -metrics.quad_height);  
      }

      function draw_spamo_divider_part(spamo, metrics, ctx, x, y, len, dir) {
        var offset;
        ctx.beginPath();
        ctx.moveTo(x, y);
        offset = 1;
        while (offset <= len) {
          ctx.lineTo(x, y + offset * dir);
          offset += 2;
          ctx.moveTo(x, y + offset * dir);
          offset += 1;
        }

        ctx.stroke();
      }

      function draw_spamo_divider(spamo, metrics, ctx) {
        draw_spamo_divider_part(spamo, metrics, ctx, metrics.left_quad + metrics.quad_width + 0.5, 
            metrics.height - metrics.top_quad, metrics.quad_height, -1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.right_quad - 0.5, 
            metrics.height - metrics.top_quad, metrics.quad_height, -1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.left_quad + metrics.quad_width + 0.5, 
            metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.right_quad - 0.5, 
            metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);
      }

      function draw_spamo_motif_boundary(spamo, metrics, ctx) {
        if (spamo.graph_mlength == 1) return;
        var offset = metrics.unit_width * (spamo.graph_mlength - 1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.left_quad + offset - 0.5,
            metrics.height - metrics.top_quad, metrics.quad_height, -1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.right_quad + metrics.quad_width - offset + 0.5, 
            metrics.height - metrics.top_quad, metrics.quad_height, -1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.left_quad + offset - 0.5, 
            metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);

        draw_spamo_divider_part(spamo, metrics, ctx, metrics.right_quad + metrics.quad_width - offset + 0.5, 
            metrics.height - (metrics.bottom_quad + metrics.quad_height), metrics.quad_height, 1);
      }

      function draw_spamo_ender(spamo, metrics, ctx) {
        var left = metrics.right_quad + metrics.quad_width + 0.5;
        ctx.save();
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(left, metrics.height - metrics.top_quad);
        ctx.lineTo(left, metrics.height - (metrics.top_quad + metrics.quad_height));
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(left, metrics.height - (metrics.bottom_quad + metrics.quad_height));
        ctx.lineTo(left, metrics.height - metrics.bottom_quad);
        ctx.stroke();

        ctx.restore();
      }

      function draw_spamo_x_axis_num(spamo, metrics, ctx, num) {
        ctx.font = SpamoConstants.font_label;
        var txt = "" + num;
        var width = ctx.measureText(txt).width;
        var vcenter = string_height(txt, SpamoConstants.font_label, SpamoConstants.size_label).vcenteroffset;
        ctx.save();
        ctx.rotate(-(Math.PI / 2));
        ctx.fillText(txt, -width - 0.5, vcenter);
        ctx.restore();
      }

      function draw_spamo_x_axis(spamo, metrics, ctx) {
        var x = metrics.left_quad + metrics.quad_width;
        var y = metrics.bottom_quad + metrics.quad_height + SpamoConstants.spacer_base + metrics.spacing_width;
        var skip = Math.ceil(SpamoConstants.size_label / metrics.unit_width);
        if (skip == 0) skip = 1;
        if (skip % 10 > 0) {
          skip -= (skip % 10) + 10;
        }
        ctx.save();
        ctx.translate(x - 0.5, metrics.height - y);
        for (var num = 0; num < spamo.graph_margin; num++) {
          if (num % skip == 0) {
            draw_spamo_x_axis_num(spamo, metrics, ctx, num);
          }
          ctx.translate(-metrics.unit_width, 0);
        }
        ctx.restore();
        x = metrics.right_quad;
        ctx.save();
        ctx.translate(x + 0.5, metrics.height - y);
        for (var num = 0; num < spamo.graph_margin; num++) {
          if (num % skip == 0) {
            draw_spamo_x_axis_num(spamo, metrics, ctx, num);
          }
          ctx.translate(metrics.unit_width, 0);
        }
        ctx.restore();

      }


      function draw_spamo_bar(spamo, metrics, ctx, bar_width, bin) {
        var bar_height = calculate_spamo_bar_height(spamo, metrics, bin.count);
        if (bin.pvalue <= spamo.graph_significant) {
          ctx.strokeStyle = "rgb(205, 0, 0)";
        } else {
          ctx.strokeStyle = "rgb(205, 205, 205)";
        }
        //this is easier to think of as lines than rectangles
        ctx.lineWidth = bar_width;
        ctx.beginPath();
        ctx.moveTo(bar_width / 2, 0);
        ctx.lineTo(bar_width / 2, -bar_height);
        ctx.stroke();
      }

      function draw_spamo_bars_part(spamo, metrics, ctx, bins) {
        var i;
        //skip unused
        if (spamo.graph_mlength > 1) ctx.translate(metrics.unit_width * (spamo.graph_mlength - 1), 0);
        i = bins.length -1;
        //draw partial bar
        if (metrics.partial_bar_width > 0) {
          draw_spamo_bar(spamo, metrics, ctx, metrics.partial_bar_width, bins[i--]);
          ctx.translate(metrics.partial_bar_width, 0);
        }
        //draw full bar
        for (; i >= 0; i--) {
          draw_spamo_bar(spamo, metrics, ctx, metrics.full_bar_width, bins[i]);
          ctx.translate(metrics.full_bar_width, 0);
        }

      }

      function draw_spamo_bars(spamo, metrics, ctx) {
        // bottom right
        ctx.save();
        ctx.translate(metrics.right_quad + metrics.quad_width, metrics.height - (metrics.bottom_quad + metrics.quad_height));
        ctx.scale(-1, -1);
        draw_spamo_bars_part(spamo, metrics, ctx, spamo.bins_opposite_right);
        ctx.restore();
        // bottom left
        ctx.save();
        ctx.translate(metrics.left_quad, metrics.height - (metrics.bottom_quad + metrics.quad_height));
        ctx.scale(1, -1);
        draw_spamo_bars_part(spamo, metrics, ctx, spamo.bins_opposite_left);
        ctx.restore();
        // top right
        ctx.save();
        ctx.translate(metrics.right_quad + metrics.quad_width, metrics.height - metrics.top_quad);
        ctx.scale(-1, 1);
        draw_spamo_bars_part(spamo, metrics, ctx, spamo.bins_same_right);
        ctx.restore();
        // top left
        ctx.save();
        ctx.translate(metrics.left_quad, metrics.height - metrics.top_quad);
        draw_spamo_bars_part(spamo, metrics, ctx, spamo.bins_same_left);
        ctx.restore();
      }


      function draw_spamo_on_canvas(spamo, canvas) {
        var ctx = canvas.getContext('2d');
        var metrics = new SpamoMetrics(spamo, canvas, ctx);
        
        if (spamo.graph_has_border) draw_spamo_border(spamo, metrics, ctx);
        if (spamo.graph_has_title) draw_spamo_title(spamo, metrics, ctx);
        if (spamo.graph_has_subtitle) draw_spamo_subtitle(spamo, metrics, ctx);
        if (spamo.graph_has_y_axis_label) draw_spamo_y_axis_label(spamo, metrics, ctx);
        if (spamo.graph_has_x_axis_label) draw_spamo_x_axis_label(spamo, metrics, ctx);
        if (spamo.graph_has_strand_label) draw_spamo_strand_labels(spamo, metrics, ctx);
        if (spamo.graph_has_stream_label) draw_spamo_stream_labels(spamo, metrics, ctx);
        draw_spamo_y_axis(spamo, metrics, ctx);
        draw_spamo_x_axis(spamo, metrics, ctx);
        draw_spamo_bars(spamo, metrics, ctx);
        draw_spamo_ender(spamo, metrics, ctx);
        draw_spamo_divider(spamo, metrics, ctx);
        draw_spamo_motif_boundary(spamo, metrics, ctx);
      }
