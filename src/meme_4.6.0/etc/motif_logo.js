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
      var green = "rgb(0,128,0)";
      var yellow = "rgb(255,255,0)";
      var purple = "rgb(204,0,204)";
      var magenta = "rgb(255,0,255)";
      var pink = "rgb(255,204,204)";
      var turquoise = "rgb(51,230,204)";
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
            return magenta;
          case "K":
          case "R":
            return red;
          case "H":
            return pink;
          case "G":
            return orange;
          case "P":
            return yellow;
          case "Y":
            return turquoise;
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
    function Pspm(pssm, name, ltrim, rtrim) {
      if (ltrim === undefined) ltrim = 0;
      if (rtrim === undefined) rtrim = 0;
      //variable prototype
      this.alph_length = 0;
      this.motif_length = 0;
      this.pspm = new Array();
      this.name = (typeof name == "string" ? name : "");
      this.nsites = 0;
      this.evalue = 0;
      this.ltrim = ltrim;
      this.rtrim = rtrim;
      //function prototype
      this.copy = Pspm_copy;
      this.reverse_complement = Pspm_reverse_complement;
      this.get_stack = Pspm_get_stack;
      this.get_stack_ic = Pspm_get_stack_ic;
      this.get_motif_length = Pspm_get_motif_length;
      this.get_alph_length = Pspm_get_alph_length;
      this.get_left_trim = Pspm_get_left_trim;
      this.get_right_trim = Pspm_get_right_trim;
      this.toString = Pspm_to_string;
      //construct
      var pspm_header = /letter-probability matrix:\s+alength=\s+(\d+)\s+w=\s+(\d+)(\s+nsites=\s+(\S+))?(\s+E=\s+(\S+))?\s*/;
      var is_empty = /^\s*$/;
      var lines = pssm.split(/\s*\n\s*/);
      var read_pssm = false;
      var line_num = 0;
      var col_num = 0;
      for (line_index in lines) {
        //exclude inherited properties and undefined properties
        if (!lines.hasOwnProperty(line_index) || lines[line_index] === undefined) continue;

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
            if (header_match[4]) this.nsites = parseFloat(header_match[4]);//not always an integer
            if (header_match[6]) this.evalue = parseFloat(header_match[6]);
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
          //exclude inherited properties and undefined properties
          if (!parts.hasOwnProperty(part_index) || parts[part_index] === undefined) continue;
          
          var prob = parts[part_index];
          if (col_num >= this.alph_length) {
            throw "TOO_MANY_COLS";
          }
          this.pspm[line_num][col_num] = (+prob);
          //check the probability is within bounds
          if (this.pspm[line_num][col_num] > 1 || this.pspm[line_num][col_num] < 0) {
            throw "NUM_NOT_PROB";
          }
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
      clone.name = (""+this.name);
      clone.nsites = (0+this.nsites);
      clone.evalue = (0+this.evalue);
      clone.ltrim = (0+this.ltrim);
      clone.rtrim = (0+this.rtrim);
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
      //swap triming
      var temp_trim = this.ltrim;
      this.ltrim = this.rtrim;
      this.rtrim = temp_trim;
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
      return alphabet.get_ic() - H;
    }

    function Pspm_get_error(alphabet) {
      var asize;
      if (this.nsites == 0) return 0;
      if (alphabet.is_nucleotide()) {
        asize = 4;
      } else {
        asize = 20;
      }
      return (asize-1) / (2 * Math.log(2)*this.nsites);
    }

    function Pspm_get_motif_length() {
      return this.motif_length;
    }

    function Pspm_get_alph_length() {
      return this.alph_length;
    }

    function Pspm_get_left_trim() {
      return this.ltrim;
    }

    function Pspm_get_right_trim() {
      return this.rtrim;
    }

    function Pspm_to_string() {
      var str = "";
      for (row_index in this.pspm) {
        //exclude inherited properties and undefined properties
        if (!this.pspm.hasOwnProperty(row_index) || this.pspm[row_index] === undefined) continue;
        
        var row = this.pspm[row_index];
        str += row.join("\t") + "\n";
      }
      return str;
    }
    //======================================================================
    // end Pspm object
    //======================================================================
    
    //======================================================================
    // start Logo object
    //======================================================================
    function Logo(alphabet, fine_text) {
      this.alphabet = alphabet;
      this.fine_text = fine_text;
      this.pspm_list = [];
      this.pspm_column = [];
      this.rows = 0;
      this.columns = 0;

      //functions
      this.add_pspm = Logo_add_pspm;
      this.get_columns = Logo_get_columns;
      this.get_rows = Logo_get_rows;
      this.get_pspm = Logo_get_pspm;
      this.get_offset = Logo_get_offset;
    }

    function Logo_add_pspm(pspm, column) {
      if (column === undefined) column = 0;
      else if (column < 0) throw "COLUMN_OUT_OF_BOUNDS";
      this.pspm_list[this.rows] = pspm;
      this.pspm_column[this.rows] = column;
      this.rows++;
      var col = column + pspm.get_motif_length();
      if (col > this.columns) this.columns = col;
    }

    function Logo_get_columns() {
      return this.columns;
    }

    function Logo_get_rows() {
      return this.rows;
    }

    function Logo_get_pspm(row_index) {
      if (row_index < 0 || row_index >= this.rows) throw "INDEX_OUT_OF_BOUNDS";
      return this.pspm_list[row_index];
    }

    function Logo_get_offset(row_index) {
      if (row_index < 0 || row_index >= this.rows) throw "INDEX_OUT_OF_BOUNDS";
      return this.pspm_column[row_index];
    }

    //======================================================================
    // end Logo object
    //======================================================================

    //======================================================================
    // start RasterizedAlphabet
    //======================================================================

    // Rasterize Alphabet
    // 1) Measure width of text at default font for all symbols in alphabet
    // 2) sort in width ascending
    // 3) Drop the top and bottom 10% (designed to ignore outliers like 'W' and 'I')
    // 4) Calculate the average as the maximum scaling factor (designed to stop I becoming a rectangular blob).
    // 5) Assume scale of zero would result in width of zero, interpolate scale required to make perfect width font
    // 6) Draw text onto temp canvas at calculated scale
    // 7) Find bounds of drawn text
    // 8) Paint on to another canvas at the desired height (but only scaling width to fit if larger).
    function RasterizedAlphabet(alphabet, font, target_width) {
      //variable prototypes
      this.lookup = []; //a map of letter to index
      this.rasters = []; //a list of rasters
      this.dimensions = []; //a list of dimensions

      //function prototypes
      this.draw = RasterizedAlphabet_draw;

      //construct
      var default_size = 60; // size of square to assume as the default width
      var safety_pad = 20; // pixels to pad around so we don't miss the edges
      // create a canvas to do our rasterizing on
      var canvas = document.createElement("canvas");
      // assume the default font would fit in a canvas of 100 by 100
      canvas.width = default_size + 2 * safety_pad;
      canvas.height = default_size + 2 * safety_pad;
      // check for canvas support before attempting anything
      if (!canvas.getContext) throw "NO_CANVAS_SUPPORT";
      var ctx = canvas.getContext('2d');
      // check for html5 text drawing support
      if (!supports_text(ctx)) throw "NO_CANVAS_TEXT_SUPPORT";
      // calculate the middle
      var middle = Math.round(canvas.width / 2);
      // calculate the baseline
      var baseline = Math.round(canvas.height - safety_pad);
      // list of widths
      var widths = [];
      var count = 0;
      var letters = [];
      //now measure each letter in the alphabet
      for (var i = 0; i < alphabet.get_size(); ++i) {
        if (alphabet.is_ambig(i)) continue; //skip ambigs as they're never rendered
        var letter = alphabet.get_letter(i);
        letters.push(letter);
        var pos = count++;
        this.lookup[letter] = pos;
        //clear the canvas
        canvas.width = canvas.width;
        // get the context and prepare to draw our width test
        var ctx = canvas.getContext('2d');
        ctx.font = font;
        ctx.fillStyle = alphabet.get_colour(i);
        ctx.textAlign = "center";
        ctx.translate(middle, baseline);
        // draw the test text
        ctx.fillText(letter, 0, 0);
        //measure
        var size = RasterizedAlphabet_measure(ctx, canvas.width, canvas.height);
        if (size.width == 0) throw "INVISIBLE_LETTER"; //maybe the fill was white on white?
        widths.push(size.width);
        this.dimensions[pos] = size;
      }
      //sort the widths
      widths.sort(function(a,b) {return a - b;});
      //drop 10% of the items off each end
      var tenpercent = Math.floor(widths.length / 10);
      for (var i = 0; i < tenpercent; ++i) {
        widths.pop();
        widths.shift();
      }
      //calculate average width
      var avg_width = 0;
      for (var i = 0; i < widths.length; ++i) avg_width += widths[i];
      avg_width /= widths.length;
      // calculate scales
      for (var i = 0; i < this.dimensions.length; ++i) {
        var size = this.dimensions[i];
        // calculate scale
        var scale = target_width / Math.max(avg_width, size.width);
        // estimate scaled height
        var target_height = size.height * scale;
        // create an approprately sized canvas
        var raster = document.createElement("canvas");
        raster.width = target_width; // if it goes over the edge too bad...
        raster.height = target_height + safety_pad * 2;
        // calculate the middle
        middle = Math.round(raster.width / 2);
        // calculate the baseline
        baseline = Math.round(raster.height - safety_pad);
        // get the context and prepare to draw the rasterized text
        ctx = raster.getContext('2d');
        ctx.font = font;
        ctx.fillStyle = alphabet.get_colour(i);
        ctx.textAlign = "center";
        ctx.translate(middle, baseline);
        ctx.save();
        ctx.scale(scale, scale);
        // draw the rasterized text
        ctx.fillText(letters[i], 0, 0);
        ctx.restore();
        this.rasters[i] = raster;
        this.dimensions[i] = RasterizedAlphabet_measure(ctx, raster.width, raster.height);
      }
    }

    function RasterizedAlphabet_measure(ctx, cwidth, cheight) {
      var data = ctx.getImageData(0, 0, cwidth, cheight).data;
      var r = 0, c = 0;// r: row, c: column
      var top_line = -1, bottom_line = -1, left_line = -1, right_line = -1;
      var txt_width = 0, txt_height = 0;
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
        txt_height = bottom_line - top_line + 1;
      }

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
        txt_width = right_line - left_line + 1;
      }

      //return the bounds
      return {bound_top: top_line, bound_bottom: bottom_line, bound_left: left_line, bound_right: right_line, width: txt_width, height: txt_height};
    }

    function RasterizedAlphabet_draw(ctx, letter, dx, dy, dWidth, dHeight) {
      var index = this.lookup[letter];
      var raster = this.rasters[index];
      var size = this.dimensions[index];
      ctx.drawImage(raster, 0, size.bound_top -1, raster.width, size.height+1, dx, dy, dWidth, dHeight);
    }

    //======================================================================
    // end RasterizedAlphabet
    //======================================================================

    //======================================================================
    // start LogoMetrics object
    //======================================================================

    function LogoMetrics(ctx, canvas_width, canvas_height, logo_columns, logo_rows, allow_space_for_names) {
      if (allow_space_for_names === undefined) allow_space_for_names = false;
      //variable prototypes
      this.canvas_width = canvas_width;
      this.canvas_height = canvas_height;
      this.scale_x = 1;
      this.scale_y = 1;
      this.pad_top = 5;
      this.pad_left = 10;
      this.pad_right = 5;
      this.pad_bottom = 0;
      this.pad_middle = 20;
      this.name_height = 14;
      this.name_font = "bold " + this.name_height + "px Times, sans-serif";
      this.name_spacer = 0;
      this.y_label = "bits"
      this.y_label_height = 12;
      this.y_label_font = "bold " + this.y_label_height + "px Helvetica, sans-serif";
      this.y_label_spacer = 3;
      this.y_num_height = 12;
      this.y_num_width = 0;
      this.y_num_font = "bold " + this.y_num_height + "px Helvetica, sans-serif";
      this.y_tic_width = 5;
      this.stack_pad_left = 0;
      this.stack_font = "bold 25px Helvetica, sans-serif";
      this.stack_height = 90;
      this.stack_width = 26;
      this.stacks_pad_right = 5;
      this.x_num_above = 2;
      this.x_num_height = 12;
      this.x_num_width = 0;
      this.x_num_font = "bold " + this.x_num_height + "px Helvetica, sans-serif";
      this.fine_txt_height = 6;
      this.fine_txt_above = 2;
      this.fine_txt_font = "normal " + this.fine_txt_height + "px Helvetica, sans-serif";
      this.letter_metrics = new Array();
      this.summed_width = 0;
      this.summed_height = 0;
      //function prototypes
      //none
      //calculate the width of the y axis numbers
      ctx.font = this.y_num_font;
      for (var i = 0; i <= 2; i++) {
        this.y_num_width = Math.max(this.y_num_width, ctx.measureText("" + i).width);
      }
      //calculate the width of the x axis numbers (but they are rotated so it becomes height)
      ctx.font = this.x_num_font;
      for (var i = 1; i <= logo_columns; i++) {
        this.x_num_width = Math.max(this.x_num_width, ctx.measureText("" + i).width);
      }
      
      //calculate how much vertical space we want to draw this
      //first we add the padding at the top and bottom since that's always there
      this.summed_height += this.pad_top + this.pad_bottom;
      //all except the last row have the same amount of space allocated to them
      if (logo_rows > 1) {
        var row_height = this.stack_height + this.pad_middle;
        if (allow_space_for_names) {
          row_height += this.name_height;
          //the label is allowed to overlap into the spacer
          row_height += Math.max(this.y_num_height/2, this.name_spacer); 
          //the label is allowed to overlap the space used by the other label
          row_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above); 
        } else {
          row_height += this.y_num_height/2; 
          //the label is allowed to overlap the space used by the other label
          row_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above); 
        }
        this.summed_height += row_height * (logo_rows - 1);
      }
      //the last row has the name and fine text below it but no padding
      this.summed_height += this.stack_height + this.y_num_height/2;
      if (allow_space_for_names) {
        this.summed_height += this.fine_txt_height + this.fine_txt_above + this.name_height;
        this.summed_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above + this.name_spacer);
      } else {
        this.summed_height += Math.max(this.y_num_height/2, this.x_num_height + this.x_num_above + this.fine_txt_height + this.fine_txt_above);
      }

      //calculate how much horizontal space we want to draw this
      //first add the padding at the left and right since that's always there
      this.summed_width += this.pad_left + this.pad_right;
      //add on the space for the y-axis label
      this.summed_width += this.y_label_height + this.y_label_spacer;
      //add on the space for the y-axis
      this.summed_width += this.y_num_width + this.y_tic_width;
      //add on the space for the stacks
      this.summed_width += (this.stack_pad_left + this.stack_width) * logo_columns;
      //add on the padding after the stacks (an offset from the fine text)
      this.summed_width += this.stacks_pad_right;

      //calculate scaling factors
      this.scale_y = this.canvas_height / this.summed_height;
      this.scale_x = this.canvas_width / this.summed_width;

      //maintain aspect ratio
      if (this.scale_y > this.scale_x) {
        this.scale_y = this.scale_x;
      } else {
        this.scale_x = this.scale_y;
      }


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
      ctx.lineWidth = 1.5;
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
      ctx.textBaseline = "middle";
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
      ctx.translate(metrics.stack_width / 2, metrics.stack_height + metrics.x_num_above);
      ctx.save();
      ctx.rotate(-(Math.PI / 2));
      ctx.textBaseline = "middle"
      ctx.textAlign = "right"
      ctx.fillText("" + (row_index + 1), 0, 0);
      ctx.restore();
      ctx.restore();
    }

    function draw_stack(ctx, metrics, symbols, raster) {
      var preferred_pad = 0;
      var sym_min = 5;

      ctx.save();//1
      ctx.translate(0, metrics.stack_height);
      for (var i in symbols) {
        //exclude inherited properties and undefined properties
        if (!symbols.hasOwnProperty(i) || symbols[i] === undefined) continue;
        
        var sym = symbols[i];
        var sym_height = metrics.stack_height * sym.get_scale();
        
        var pad = preferred_pad;
        if (sym_height - pad < sym_min) {
          pad = Math.min(pad, Math.max(0, sym_height - sym_min));
        }
        sym_height -= pad;

        //translate to the correct position
        ctx.translate(0, -(pad/2 + sym_height));
        //draw
        raster.draw(ctx, sym.get_symbol(), 0, 0, metrics.stack_width, sym_height);
        //translate past the padding
        ctx.translate(0, -(pad/2));
      }
      ctx.restore();//1
    }

    //draws a stack of symbols
    function draw_stack_old(ctx, metrics, symbols) {
      var lpad = 2;
      var sym_min = 5;
      var pos = metrics.stack_height;
      for (var i in symbols) {
        //exclude inherited properties and undefined properties
        if (!symbols.hasOwnProperty(i) || symbols[i] === undefined) continue;
        
        var sym = symbols[i];
        var sym_height = metrics.stack_height*sym.get_scale();
        var letter = metrics.get_letter_metrics(sym.get_symbol());
        //attempting to draw something smaller than a pixel causes display corruption
        if (sym_height >= 1) {
          //it's better to see the letter than to pad it
          var pad = lpad;
          if (sym_height - pad < sym_min) {
            pad = Math.min(pad, Math.max(0, sym_height - sym_min));
          }
          //move to the correct drawing position
          ctx.save();//s1
          ctx.translate(0, pos);
          //create a clipping rectangle to ensure the letter doesn't overlap when it's distorted
          ctx.save();//s2
          //ctx.beginPath(); //disabled clipping because after the improvements in the text metrics it looks better without
          //ctx.moveTo(-metrics.stack_width/2,0);
          //ctx.lineTo(metrics.stack_width/2, 0);
          //ctx.lineTo(metrics.stack_width/2, -sym_height);
          //ctx.lineTo(-metrics.stack_width/2, -sym_height);
          //ctx.lineTo(-metrics.stack_width/2,0);
          //ctx.clip();
          //now draw
          ctx.translate(0, -(pad/2));
          ctx.translate(0, -letter.get_descent(sym_height - pad));
          ctx.fillStyle = sym.get_colour();
          ctx.textAlign = "center";
          ctx.save();//s3
          ctx.scale(letter.wscale, letter.get_hscale(sym_height - pad));
          ctx.fillText(sym.get_symbol(), 0, 0);
          ctx.restore();//s3

          ctx.restore();//s2
          ctx.restore();//s1
        }

        pos = pos - sym_height;
      }
    }
    
    function draw_dashed_line(ctx, pattern, start, x1, y1, x2, y2) {
        var x, y, len, i;
        var dx = x2 - x1;
        var dy = y2 - y1;
        var tlen = Math.pow(dx*dx + dy*dy, 0.5);
        var theta = Math.atan2(dy,dx);
        var mulx = Math.cos(theta);
        var muly = Math.sin(theta);
        var lx = [];
        var ly = [];
        for (i = 0; i < pattern; ++i) {
          lx.push(pattern[i] * mulx);
          ly.push(pattern[i] * muly);
        }
        i = start;
        x = x1;
        y = y1;
        len = 0;
        ctx.beginPath();
        while (len + pattern[i] < tlen) {
          ctx.moveTo(x, y);
          x += lx[i];
          y += ly[i];
          ctx.lineTo(x, y);
          len += pattern[i];
          i = (i + 1) % pattern.length;
          x += lx[i];
          y += ly[i];
          len += pattern[i];
          i = (i + 1) % pattern.length;
        }
        if (len < tlen) {
          ctx.moveTo(x, y);
          x += mulx * (tlen - len);
          y += muly * (tlen - len);
          ctx.lineTo(x, y);
        }
        ctx.stroke();
    }

    function draw_trim_background(ctx, metrics, pspm, offset) {
      var lwidth = metrics.stack_width * pspm.get_left_trim();
      var rwidth = metrics.stack_width * pspm.get_right_trim();
      var mwidth = metrics.stack_width * pspm.get_motif_length();
      var rstart = mwidth - rwidth;
      ctx.save();//s8
      ctx.translate(offset * metrics.stack_width, 0);
      ctx.fillStyle = "rgb(240, 240, 240)";
      if (pspm.get_left_trim() > 0) ctx.fillRect(0, 0, lwidth, metrics.stack_height);
      if (pspm.get_right_trim() > 0) ctx.fillRect(rstart, 0, rwidth, metrics.stack_height);
      ctx.fillStyle = "rgb(51, 51, 51)";
      if (pspm.get_left_trim() > 0) draw_dashed_line(ctx, [3], 0, lwidth-0.5, 0, lwidth-0.5,  metrics.stack_height);
      if (pspm.get_right_trim() > 0) draw_dashed_line(ctx, [3], 0, rstart+0.5, 0, rstart+0.5,  metrics.stack_height);
      ctx.restore();//s8
    }

    function draw_logo_on_canvas(logo, canvas, show_names, scale) {
      var draw_name = (typeof show_names == "boolean" ? show_names : (logo.get_rows() > 1));
      var cwidth = canvas.width;
      var cheight = canvas.height;
      //need a minimum 46 x 120 canvas to draw the font size checks on
      if (canvas.width < 46) canvas.width = 46;
      if (canvas.height < 120) canvas.height = 120;
      var ctx = canvas.getContext('2d');
      //assume that the user wants the canvas scaled equally so calculate what the best width for this image should be
      var metrics = new LogoMetrics(ctx, canvas.width, canvas.height, logo.get_columns(), logo.get_rows(), draw_name);
      ctx.save();//s1
      if (typeof scale == "number") {
        //resize the canvas to fit the scaled logo
        cwidth = metrics.summed_width * scale;
        cheight = metrics.summed_height * scale;
      } else {
        if (cwidth == 0 || cheight == 0 || scale == 0) {
          throw "CANVAS_MUST_HAVE_DIMENSIONS";
        }
        scale = Math.min(cwidth / metrics.summed_width, cheight / metrics.summed_height);
      }
      var raster = new RasterizedAlphabet(logo.alphabet, metrics.stack_font, metrics.stack_width * scale * 2);
      if (cwidth != canvas.width || cheight != canvas.height) {
        canvas.width = cwidth;
        canvas.height = cheight;
        //as the canvas has been resized the context is now out of date
        ctx = canvas.getContext('2d');
      }
      ctx.scale(scale, scale);
      ctx.save();//s2
      ctx.save();//s7
      //create margin
      ctx.translate(metrics.pad_left, metrics.pad_top);
      for (var pspm_i = 0; pspm_i < logo.get_rows(); ++pspm_i) {
        var pspm = logo.get_pspm(pspm_i);
        var offset = logo.get_offset(pspm_i);
        //optionally draw name if this isn't the last row or is the only row 
        if (draw_name && (logo.get_rows() == 1 || pspm_i != (logo.get_rows()-1))) {
          ctx.save();//s4
          ctx.translate(metrics.summed_width/2, metrics.name_height);
          ctx.font = metrics.name_font;
          ctx.textAlign = "center";
          ctx.fillText(pspm.name, 0, 0);
          ctx.restore();//s4
          ctx.translate(0, metrics.name_height + Math.min(0, metrics.name_spacer - metrics.y_num_height/2));
        }
        //draw scale
        draw_scale(ctx, metrics, logo.alphabet.get_ic());
        ctx.save();//s5
        //translate across past the scale
        ctx.translate(metrics.y_label_height + metrics.y_label_spacer + 
            metrics.y_num_width + metrics.y_tic_width, 0);
        //draw the trimming background
        if (pspm.get_left_trim() > 0 || pspm.get_right_trim() > 0) {
          draw_trim_background(ctx, metrics, pspm, offset);
        }
        //draw letters
        ctx.translate(0, metrics.y_num_height / 2);
        for (var col_index = 0; col_index < logo.get_columns(); col_index++) {
          ctx.translate(metrics.stack_pad_left,0);
          if (col_index >= offset && col_index < (offset + pspm.get_motif_length())) {
            var motif_position = col_index - offset;
            draw_stack_num(ctx, metrics, motif_position);
            draw_stack(ctx, metrics, pspm.get_stack(motif_position, logo.alphabet), raster);
          }
          ctx.translate(metrics.stack_width, 0);
        }
        ctx.restore();//s5
        ////optionally draw name if this is the last row but isn't the only row 
        if (draw_name && (logo.get_rows() != 1 && pspm_i == (logo.get_rows()-1))) {
          //translate vertically past the stack and axis's        
          ctx.translate(0, metrics.y_num_height/2 + metrics.stack_height + 
              Math.max(metrics.y_num_height/2, metrics.x_num_above + metrics.x_num_width + metrics.name_spacer));

          ctx.save();//s6
          ctx.translate(metrics.summed_width/2, metrics.name_height);
          ctx.font = metrics.name_font;
          ctx.textAlign = "center";
          ctx.fillText(pspm.name, 0, 0);
          ctx.restore();//s6
          ctx.translate(0, metrics.name_height);
        } else {
          //translate vertically past the stack and axis's        
          ctx.translate(0, metrics.y_num_height/2 + metrics.stack_height + Math.max(metrics.y_num_height/2, metrics.x_num_above + metrics.x_num_width));
        }
        //if not the last row then add middle padding
        if (pspm_i != (logo.get_rows() -1)) {
          ctx.translate(0, metrics.pad_middle);
        }
      }
      ctx.restore();//s7
      ctx.translate(metrics.summed_width - metrics.pad_right, metrics.summed_height - metrics.pad_bottom);
      ctx.font = metrics.fine_txt_font;
      ctx.textAlign = "right";
      ctx.fillText(logo.fine_text, 0,0);
      ctx.restore();//s2
      ctx.restore();//s1
    }

    function create_canvas(c_width, c_height, c_id, c_title, c_display) {
      var canvas = document.createElement("canvas");
      //check for canvas support before attempting anything
      if (!canvas.getContext) return null;
      var ctx = canvas.getContext('2d');
      //check for html5 text drawing support
      if (!supports_text(ctx)) return null;
      //size the canvas
      canvas.width = c_width;
      canvas.height = c_height;
      canvas.id = c_id;
      canvas.title = c_title;
      canvas.style.display = c_display;
      return canvas;
    }

    function logo_1(alphabet, fine_text, pspm) {
      var logo = new Logo(alphabet, fine_text);
      logo.add_pspm(pspm);
      return logo;
    }
    
    function logo_2(alphabet, fine_text, target, query, query_offset) {
      var logo = new Logo(alphabet, fine_text);
      if (query_offset < 0) {
        logo.add_pspm(target, -query_offset);
        logo.add_pspm(query);
      } else {
        logo.add_pspm(target);
        logo.add_pspm(query, query_offset);
      }      
      return logo;
    }

    /*
     * Specifies an alternate source for an image.
     * If the image with the image_id specified has
     * not loaded then a generated logo will be used 
     * to replace it.
     *
     * Note that the image must either have dimensions
     * or a scale must be set.
     */
    function alternate_logo(logo, image_id, scale) {
      var image = document.getElementById(image_id);
      if (!image) {
        alert("Can't find specified image id (" +  image_id + ")");
        return;
      }
      //if the image has loaded then there is no reason to use the canvas
      if (image_ok(image)) return;
      //the image has failed to load so replace it with a canvas if we can.
      var canvas = create_canvas(image.width, image.height, image_id, image.title, image.style.display);
      if (canvas == null) return;
      //draw the logo on the canvas
      draw_logo_on_canvas(logo, canvas, undefined, scale);
      //replace the image with the canvas
      image.parentNode.replaceChild(canvas, image);
    }

    /*
     * Specifes that the element with the specified id
     * should be replaced with a generated logo.
     */
    function replace_logo(logo, replace_id, scale, title_txt, display_style) {
      var element = document.getElementById(replace_id);
      if (!replace_id) {
        alert("Can't find specified id (" + replace_id + ")");
        return;
      }
      //found the element!
      var canvas = create_canvas(50, 120, replace_id, title_txt, display_style);
      if (canvas == null) return;
      //draw the logo on the canvas
      draw_logo_on_canvas(logo, canvas, undefined, scale);
      //replace the element with the canvas
      element.parentNode.replaceChild(canvas, element);
    }

