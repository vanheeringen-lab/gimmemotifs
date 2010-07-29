    //======================================================================
    // start Alphabet object
    //======================================================================
    function Alphabet(alphabet, bg) {
      //variable prototype
      this.freqs = new Array();
      this.alphabet = new Array();
      this.letter_count = 0;
      //method prototype
      this.get_ic = Alphabet_get_ic;
      this.get_size = Alphabet_get_size;
      this.get_index = Alphabet_get_index;
      this.get_letter = Alphabet_get_letter;
      this.get_colour = Alphabet_get_colour;
      this.get_bg_freq = Alphabet_get_bg_freq;
      this.is_nucleotide = Alphabet_is_nucleotide;
      this.is_ambig = Alphabet_is_ambig;
      this.toString = Alphabet_to_string;
      //construct
      var is_letter = /^\w$/;
      var is_prob = /^((1(\.0+)?)|(0(\.\d+)?))$/;
      for (var pos = 0; pos < alphabet.length; pos++) {
        var letter = alphabet.charAt(pos);
        if (is_letter.test(letter)) {
          this.alphabet[this.letter_count] = letter.toUpperCase();
          this.freqs[this.letter_count] = -1;
          this.letter_count++;
        }
      }
      if (!(bg === undefined)) {
        var parts = bg.split(/\s+/);
        for (var i = 0, pos = 0; (i + 1) < parts.length; i += 2) {
          var letter = parts[i];
          var freq = parts[i+1];
          if (is_letter.test(letter) && is_prob.test(freq)) {
            letter = letter.toUpperCase();          //find the letter it matches
            for (;pos < this.letter_count; pos++) {
              if (this.alphabet[pos] == letter) break;
            }
            if (pos >= this.letter_count) throw "NOT_IN_ALPHABET";
            this.freqs[pos] = (+freq);
          }
        }
      }
    }


    function Alphabet_get_ic() {
      if (this.is_nucleotide()) {
        return 2;
      } else {
        return Math.log(20) / Math.LN2;
      }
    }

    function Alphabet_get_size() {
      return this.letter_count;
    }

    function Alphabet_get_letter(alph_index) {
      if (alph_index < 0 || alph_index >= this.letter_count) {
        throw "BAD_ALPHABET_INDEX";
      }
      return this.alphabet[alph_index];
    }

    function Alphabet_get_bg_freq(alph_index) {
      if (alph_index < 0 || alph_index >= this.letter_count) {
        throw "BAD_ALPHABET_INDEX";
      }
      if (this.freqs[alph_index] == -1) throw "BG_FREQ_NOT_SET";
      return this.freqs[alph_index];
    }

    function Alphabet_get_colour(alph_index) {
      var red = "rgb(204,0,0)";
      var blue = "rgb(0,0,204)";
      var orange = "rgb(255,179,0)";
      var green = "rgb(0,127,0)";
      if (alph_index < 0 || alph_index >= this.letter_count) {
        throw "BAD_ALPHABET_INDEX";
      }
      if (this.is_nucleotide()) {
        switch (this.alphabet[alph_index]) {
          case "A":
            return red;
          case "C":
            return blue;
          case "G":
            return orange;
          case "T":
            return green;
        }
      } else {
        switch (this.alphabet[alph_index]) {
          case "A":
          case "C":
          case "F":
          case "I":
          case "L":
          case "V":
          case "W":
          case "M":
            return blue;
          case "N":
          case "Q":
          case "S":
          case "T":
            return green;
          case "D":
          case "E":
            return "magenta";
          case "K":
          case "R":
            return red;
          case "H":
            return "pink";
          case "G":
            return orange;
          case "P":
            return "yellow";
          case "Y":
            return "turquoise";
        }
      }
      return "black";
    }

    function Alphabet_is_ambig(alph_index) {
      if (alph_index < 0 || alph_index >= this.letter_count) {
        throw "BAD_ALPHABET_INDEX";
      }
      if (this.is_nucleotide()) {
        return ("ACGT".indexOf(this.alphabet[alph_index]) == -1);
      } else {
        return ("ACDEFGHIKLMNPQRSTVWY".indexOf(this.alphabet[alph_index]) == -1);
      }
    }

    function Alphabet_get_index(letter) {
      for (i = 0; i < this.letter_count; i++) {
        if (this.alphabet[i] == letter.toUpperCase()) return i;
      }
      throw "UNKNOWN_LETTER";
    }

    function Alphabet_is_nucleotide() {
      //TODO basic method, make better
      if (this.letter_count < 20) return true;
      return false;
    }

    function Alphabet_to_string() {
      return (this.is_nucleotide() ? "Nucleotide" : "Protein") + " Alphabet " + (this.alphabet.join(""));
    }

    //======================================================================
    // end Alphabet object
    //======================================================================

    //======================================================================
    // start Symbol object
    //======================================================================
    function Symbol(alph_index, scale, alphabet) {
      //variable prototype
      this.symbol = alphabet.get_letter(alph_index);
      this.scale = scale;
      this.colour = alphabet.get_colour(alph_index);
      //function prototype
      this.get_symbol = Symbol_get_symbol;
      this.get_scale = Symbol_get_scale;
      this.get_colour = Symbol_get_colour;
      this.toString = Symbol_to_string;
    }

    function Symbol_get_symbol() {
      return this.symbol;
    }

    function Symbol_get_scale() {
      return this.scale;
    }

    function Symbol_get_colour() {
      return this.colour;
    }

    function Symbol_to_string() {
      return this.symbol + " " + (Math.round(this.scale*1000)/10) + "%";
    }

    function compare_symbol(sym1, sym2) {
      if (sym1.get_scale() < sym2.get_scale()) {
        return -1;
      } else if (sym1.get_scale() > sym2.get_scale()) {
        return 1;
      } else {
        return 0;
      }
    }
    //======================================================================
    // end Symbol object
    //======================================================================

    //======================================================================
    // start Pspm object
    //======================================================================
    function Pspm(pssm) {
      //variable prototype
      this.alph_length = 0;
      this.motif_length = 0;
      this.pspm = new Array();
      //function prototype
      this.copy = Pspm_copy;
      this.reverse_complement = Pspm_reverse_complement;
      this.get_stack = Pspm_get_stack;
      this.get_stack_ic = Pspm_get_stack_ic;
      this.get_motif_length = Pspm_get_motif_length;
      this.get_alph_length = Pspm_get_alph_length;
      this.toString = Pspm_to_string;
      //construct
      var pspm_header = /letter-probability matrix:\s+alength=\s+(\d+)\s+w=\s+(\d+)\s+/;
      var is_prob = /^((1(\.0+)?)|(0(\.\d+)?))$/;
      var is_empty = /^\s*$/;
      var lines = pssm.split(/\s*\n\s*/);
      var read_pssm = false;
      var line_num = 0;
      var col_num = 0;
      for (line_index in lines) {
        var line = lines[line_index];
        if (is_empty.test(line)) {
          continue;
        }
        if (!read_pssm) {
          var header_match = pspm_header.exec(line);
          if (header_match != null) {
            read_pssm = true;
            this.alph_length = (+header_match[1]);
            this.motif_length = (+header_match[2]);
            this.pspm = new Array(this.motif_length);
          }
          continue;
        }
        if (line_num >= this.motif_length) {
          throw "TOO_MANY_ROWS";
        }
        this.pspm[line_num] = new Array(this.alph_length);
        col_num = 0;
        var parts = line.split(/\s+/);
        for (part_index in parts) {
          var prob = parts[part_index];
          if (!is_prob.test(prob)) continue;
          if (col_num >= this.alph_length) {
            throw "TOO_MANY_COLS";
          }
          this.pspm[line_num][col_num] = (+prob);
          col_num++;
        }
        if (col_num != this.alph_length) {
          throw "TOO_FEW_COLS";
        }
      line_num++;
      }
      if (line_num != this.motif_length) {
        throw "TOO_FEW_ROWS";
      }
    }

    function Clone() {}

    function Pspm_copy() {
      Clone.prototype = this;
      var clone = new Clone();
      //so far only a shallow copy, need to copy everything
      clone.alph_length = (0+this.alph_length);
      clone.motif_length = (0+this.motif_length);
      clone.pspm = new Array(this.motif_length);
      for (row = 0; row < this.motif_length; row++) {
        clone.pspm[row] = new Array(this.alph_length);
        for (col = 0; col < this.alph_length; col++) {
          clone.pspm[row][col] = (0+this.pspm[row][col]);
        }
      }
      return clone;
    }

    function Pspm_reverse_complement(alphabet) {
      if (this.alph_length != alphabet.get_size()) {
        throw "ALPHABET_MISMATCH";
      }
      if (!alphabet.is_nucleotide()) {
        throw "NO_PROTEIN_RC";
      }
      //reverse
      var x = 0;
      var y = this.motif_length-1;
      while (x < y) {
        var temp = this.pspm[x];
        this.pspm[x] = this.pspm[y];
        this.pspm[y] = temp;
        x++;
        y--;
      }
      //complement
      var a_index = alphabet.get_index("A");
      var c_index = alphabet.get_index("C");
      var g_index = alphabet.get_index("G");
      var t_index = alphabet.get_index("T");
      for (i = 0; i < this.motif_length; i++) {
        var row = this.pspm[i];
        //swap A and T
        var temp = row[a_index];
        row[a_index] = row[t_index];
        row[t_index] = temp;
        //swap C and G
        temp = row[c_index];
        row[c_index] = row[g_index];
        row[g_index] = temp;
      }
      //note that ambigs are ignored because they don't effect motifs
      return this; //allow function chaining...
    }

    function Pspm_get_stack(position, alphabet) {
      if (this.alph_length != alphabet.get_size()) {
        throw "ALPHABET_MISMATCH";
      }
      var row = this.pspm[position];
      var stack_ic = this.get_stack_ic(position, alphabet);
      var alphabet_ic = alphabet.get_ic();
      var stack = new Array();
      for (i = 0; i < this.alph_length; i++) {
        if (alphabet.is_ambig(i)) continue;
        var sym = new Symbol(i, row[i]*stack_ic/alphabet_ic, alphabet);
        if (sym.get_scale() <= 0) continue;
        stack.push(sym);
      }
      stack.sort(compare_symbol);
      return stack;
    }

    function Pspm_get_stack_ic(position, alphabet) {
      if (this.alph_length != alphabet.get_size()) {
        throw "ALPHABET_MISMATCH";
      }
      var row = this.pspm[position];
      var H = 0;
      for (var i = 0; i < this.alph_length; i++) {
        if (alphabet.is_ambig(i)) continue;
        if (row[i] == 0) continue;
        H -= (row[i] * (Math.log(row[i]) / Math.LN2));
      }
      if (alphabet.is_nucleotide()) {
        return 2 - H;
      } else {
        return (Math.log(20)/Math.LN2) - H;
      }
    }

    function Pspm_get_motif_length() {
      return this.motif_length;
    }

    function Pspm_get_alph_length() {
      return this.alph_length;
    }

    function Pspm_to_string() {
      var str = "";
      for (row_index in this.pspm) {
        var row = this.pspm[row_index];
        str += row.join("\t") + "\n";
      }
      return str;
    }
    //======================================================================
    // end Pspm object
    //======================================================================

    //======================================================================
    // start TextMetrics
    //======================================================================
    function TextMetrics(ctx, canvas_width, canvas_height, letter) {
      //variable prototypes
      this.ascent = 0;
      this.descent = 0;
      this.width = 0;
      //function prototypes
      this.get_height = TextMetrics_get_height;
      //construct
      var baseline = Math.round(canvas_height/2);
      ctx.clearRect(0, 0, canvas_width, canvas_height);
      ctx.fillText(letter, 0, baseline);
      var data = ctx.getImageData(0, 0, canvas_width, canvas_height).data;
      var r = 0, c = 0;
      var start_line = -1, finish_line = -1;
      // Find the first line with a non-white pixel
      for(r = 0; r <= baseline; r++) {
        for(c = 0; c < canvas_width; c++) {
          if(data[r * canvas_width * 4 + c * 4 + 3]) {
            start_line = r;
            break;
          }
        }
        if (start_line != -1) break;
      }
      if (start_line == -1) throw "LETTER_INVISIBLE";
      
      //find the last line with a non-white pixel
      for (r = canvas_height-1; r >= start_line; r--) {
        for(c = 0; c < canvas_width; c++) {
          if(data[r * canvas_width * 4 + c * 4 + 3]) {
            finish_line = r;
            break;
          }
        }
        if (finish_line != -1) break;
      }
      ctx.clearRect(0, 0, canvas_width, canvas_height);
      this.ascent = baseline - start_line;
      this.descent = (finish_line > baseline ? finish_line - baseline : 0);
      this.width = ctx.measureText(letter).width;
    }

    function TextMetrics_get_height() {
      return this.ascent + this.descent;
    }

    //======================================================================
    // end TextMetrics
    //======================================================================

    //======================================================================
    // start LogoMetrics object
    //======================================================================

    function LogoMetrics(ctx, canvas_width, canvas_height, alphabet, positions, use_preferred_and_scale) {
      if (use_preferred_and_scale === undefined) use_preferred_and_scale = false;
      //variable prototypes
      this.canvas_width = canvas_width;
      this.canvas_height = canvas_height;
      this.scale_x = 1;
      this.scale_y = 1;
      this.pad_top = 5;
      this.pad_left = 10;
      this.pad_right = 5;
      this.pad_bottom = 0;
      this.y_label = "bits"
      this.y_label_height = 11;
      this.y_label_font = "bold " + this.y_label_height + "px sans-serif";
      this.y_label_spacer = 3;
      this.y_num_height = 12;
      this.y_num_width = 0;
      this.y_num_font = "bold " + this.y_num_height + "px sans-serif";
      this.y_tic_width = 5;
      this.stack_pad_left = 0;
      this.stack_preferred_width = 25;
      this.stack_preferred_height = 200;
      this.stack_font = "bold 20px sans-serif";
      this.stack_height = 0;
      this.stack_width = 0;
      this.stacks_pad_right = 10;
      this.x_num_above = 2;
      this.x_num_height = 14;
      this.x_num_width = 0;
      this.x_num_font = "bold " + this.x_num_height + "px sans-serif";
      this.fine_txt_height = 7;
      this.fine_txt_above = 2;
      this.fine_txt_font = "italic " + this.fine_txt_height + "px sans-serif";
      this.letter_metrics = new Array();
      //function prototypes
      this.get_letter_metrics = LogoMetrics_get_letter_metrics;
      //calculate the width of the y axis numbers
      ctx.font = this.y_num_font;
      for (var i = 0; i <= 2; i++) {
        this.y_num_width = Math.max(this.y_num_width, ctx.measureText("" + i).width);
      }
      //calculate the width of the x axis numbers (but they are rotated so it becomes height)
      ctx.font = this.x_num_font;
      for (var i = 1; i <= positions; i++) {
        this.x_num_width = Math.max(this.x_num_width, ctx.measureText("" + i).width);
      }
      //calculate how much vertical space is left over for drawing the stack
      var extra_height = this.pad_bottom + this.fine_txt_height + this.fine_txt_above + 
          Math.max(this.y_num_height/2,
            this.fine_txt_height + this.fine_txt_above + this.x_num_height + this.x_num_above) +
          this.pad_top + this.y_num_height/2;
      if (use_preferred_and_scale) {
        this.stack_height = this.stack_preferred_height;
        this.scale_y = this.canvas_height / (extra_height + this.stack_height);
      } else {
        this.stack_height = this.canvas_height - extra_height;
      }
      if (this.stack_height <= 0) {
        throw "TOO_SMALL_FOR_STACK";
      }
      //calculate how much horizontal space is left over for drawing the stack
      var extra_width = this.pad_left + this.y_label_height + this.y_label_spacer + 
          this.y_num_width + this.y_tic_width + this.stacks_pad_right + this.pad_right;
      if (use_preferred_and_scale) {
        this.stack_width = this.stack_preferred_width;
        this.scale_x = this.canvas_width / (extra_width + (this.stack_width * positions));
      } else {
        this.stack_width = ((this.canvas_width - extra_width) / positions) - this.stack_pad_left;
      }
      //calculate the height of the letters, used for scaling
      ctx.font = this.stack_font;
      for (var i = 0; i < alphabet.get_size(); i++) {
        if (alphabet.is_ambig(i)) continue;
        this.letter_metrics[alphabet.get_letter(i)] =
            new TextMetrics(ctx, canvas_width, canvas_height, alphabet.get_letter(i));
      }
    }

    function LogoMetrics_get_letter_metrics(letter) {
      var metrics = this.letter_metrics[letter];
      if (metrics == null) {
        throw "LETTER_NOT_CALCULATED";
      }
      return metrics;
    }

    //======================================================================
    // end LogoMetrics object
    //======================================================================
    
    //found this trick at http://talideon.com/weblog/2005/02/detecting-broken-images-js.cfm
    function image_ok(img) {
      // During the onload event, IE correctly identifies any images that
      // weren't downloaded as not complete. Others should too. Gecko-based
      // browsers act like NS4 in that they report this incorrectly.
      if (!img.complete) {
        return false;
      }

      // However, they do have two very useful properties: naturalWidth and
      // naturalHeight. These give the true size of the image. If it failed
      // to load, either of these should be zero.
      if (typeof img.naturalWidth != "undefined" && img.naturalWidth == 0) {
        return false;
      }

      // No other way of checking: assume it's ok.
      return true;
    }
      
    function supports_text(ctx) {
      if (!ctx.fillText) return false;
      if (!ctx.measureText) return false;
      return true;
    }

    //draws the scale, returns the width
    function draw_scale(ctx, metrics, alphabet_ic) {
      var tic_height = metrics.stack_height / alphabet_ic;

      ctx.save();

      ctx.translate(metrics.y_label_height, metrics.y_num_height/2);

      //draw the axis label
      ctx.save();
      ctx.font = metrics.y_label_font;
      ctx.translate(0, metrics.stack_height/2);
      ctx.save();
      ctx.rotate(-(Math.PI / 2));
      ctx.textAlign = "center";
      ctx.fillText("bits", 0, 0);
      ctx.restore();
      ctx.restore();

      ctx.translate(metrics.y_label_spacer + metrics.y_num_width, 0);

      //draw the axis tics
      ctx.save();
      ctx.translate(0, metrics.stack_height);
      ctx.font = metrics.y_num_font;
      ctx.textAlign = "right";
      ctx.textBaseline = "middle"
      for (var i = 0; i <= alphabet_ic; i++) {
        //draw the number
        ctx.fillText("" + i, 0, 0);
        //draw the tic
        ctx.beginPath();
        ctx.moveTo(0, 0);
        ctx.lineTo(metrics.y_tic_width, 0);
        ctx.stroke();
        //prepare for next tic
        ctx.translate(0, -tic_height);
      }
      ctx.restore();

      ctx.translate(metrics.y_tic_width, 0);

      ctx.beginPath();
      ctx.moveTo(0, 0);
      ctx.lineTo(0, metrics.stack_height);
      ctx.stroke();

      ctx.restore();
    }

    function draw_stack_num(ctx, metrics, row_index) {
      ctx.save();
      ctx.font = metrics.x_num_font;
      ctx.translate(0, metrics.stack_height + metrics.x_num_above);
      ctx.save();
      ctx.rotate(-(Math.PI / 2));
      ctx.textBaseline = "middle"
      ctx.textAlign = "right"
      ctx.fillText("" + (row_index + 1), 0, 0);
      ctx.restore();
      ctx.restore();
    }

    //draws a stack of symbols
    function draw_stack(ctx, metrics, symbols) {
      var lpad = 2;
      var pos = metrics.stack_height;
      for (var i in symbols) {
        var sym = symbols[i];
        var sym_height = metrics.stack_height*sym.get_scale();
        var letter = metrics.get_letter_metrics(sym.get_symbol());
        //attempting to draw something smaller than a pixel causes display corruption
        if (sym_height >= 1) {
          //it's better to see the letter than to pad it
          var pad = lpad;
          if (sym_height - pad < 5) {
            pad = Math.max(0, sym_height - 5);
          }
          //draw letter
          ctx.save();
          ctx.translate(0, pos);
          ctx.translate(0, -(pad/2));
          ctx.fillStyle = sym.get_colour();
          ctx.textAlign = "center";
          ctx.save();
          ctx.scale(metrics.stack_width/letter.width, (sym_height - pad)/letter.get_height());
          ctx.translate(0, -letter.descent);
          ctx.fillText(sym.get_symbol(), 0, 0);
          ctx.restore();
          ctx.restore();
        }

        pos = pos - sym_height;
      }
    }

    function draw_logo(pspm, alphabet, fine_txt, image_id) {
      var image = document.getElementById(image_id);
      if (!image) {
        alert("Can't find specified image id (" +  image_id + ")");
        return;
      }
      //if the image has loaded then there is no reason to use the canvas
      if (image_ok(image)) return;
      //the image has failed to load so replace it with a canvas if we can
      var canvas = document.createElement("canvas");
      //check for canvas support before attempting anything
      if (!canvas.getContext) return;
      var ctx = canvas.getContext('2d');
      //check for html5 text drawing support
      if (!supports_text(ctx)) return;
      //size the canvas the same as the image (this means I'll have to set it in html)
      canvas.width = image.width;
      canvas.height = image.height;
      canvas.id = image.id;
      canvas.alt = image.alt;
      canvas.style.display = image.style.display;
      //replace the image with the canvas
      image.parentNode.replaceChild(canvas, image);
      //assume that the user wants the canvas scaled equally so calculate what the best width for this image should be
      var metrics = new LogoMetrics(ctx, canvas.width, canvas.height, alphabet, pspm.get_motif_length(), true);
      ctx.save();
      ctx.scale(metrics.scale_x, metrics.scale_y);
      ctx.save();
      //create margin
      ctx.translate(metrics.pad_left, metrics.pad_top);
      //draw scale
      draw_scale(ctx, metrics, alphabet.get_ic());
      ctx.translate(metrics.y_label_height + metrics.y_label_spacer + 
          metrics.y_num_width + metrics.y_tic_width, 0);
      //draw letters
      ctx.font = metrics.stack_font;
      ctx.translate(metrics.stack_width/2, metrics.y_num_height / 2);
      for (var row_index = 0; row_index < pspm.get_motif_length(); row_index++) {
        ctx.translate(metrics.stack_pad_left,0);
        draw_stack_num(ctx, metrics, row_index);
        draw_stack(ctx, metrics, pspm.get_stack(row_index, alphabet));
        ctx.translate(metrics.stack_width, 0);
      }
      ctx.translate(metrics.stacks_pad_right - metrics.stack_width/2, metrics.stack_height +
          metrics.x_num_above + metrics.x_num_width + 
          metrics.fine_txt_above + metrics.fine_txt_height);
      ctx.font = metrics.fine_txt_font;
      ctx.textAlign = "right";
      ctx.fillText(fine_txt, 0,0);
      ctx.restore();
      ctx.restore();
    }

    //example code for creating a logo
    //function do_load() {
    //  var alphabet = new Alphabet(document.getElementById("alphabet").value, document.getElementById("bgfreq").value);
    //  var pspm = new Pspm(document.getElementById("pspm1").value);
    //  var pspm_rc = pspm.copy().reverse_complement(alphabet);
    //  draw_logo(pspm, alphabet, "Motif 1", "logo1");
    //  draw_logo(pspm_rc, alphabet, "Motif 1 Reverse Complement", "logo_rc1");
    //}

