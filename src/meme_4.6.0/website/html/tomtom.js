// javascript for tomtom.cgi

var codes = {A : 1, C : 2, G : 4, T : 8, U : 8, R : 5, Y : 10, M : 3, K : 12, W : 9, S : 6, B : 14, D : 13, H : 11, V : 7, N : 15};
var preview_delay_timer = null;
var cached_alphabet = null;
var logo_error_count = 0;
var logo_warning_count = 0;
var log_box = null;
var vis_ids = [];

// setup_form
//
// Sets up the visible portions of the form to match the choices made
// in the radio buttons.
//
function setup_form() {
  //get the logging text box
  log_box = document.getElementById('motif_log');
  //check that the motif format is synchronized with the div
  if (document.getElementById('type_s').checked) motif_format('consensus');
  else if (document.getElementById('type_f').checked) motif_format('matrix');
  else if (document.getElementById('type_m').checked) motif_format('meme');
  else {
    var inline_motif = document.getElementById('type_i');
    if (inline_motif != null && inline_motif.checked) motif_format('inline');
  }
  update_preview();
  show_db_upload();
}

//
// reset_form
// 
// Sets the form back to the default state.
//
function reset_form() {
  // if there was a inline motif then select it by default
  if (document.getElementById('type_i') == null) motif_format('consensus');
  else motif_format('inline');
  document.getElementById('motif_s').value = '';
  document.getElementById('motif_f').value = '';
  update_preview();
  show_hide('', 'queue_div');
  show_hide('', 'db_upload_div');
}

//
// check_form
//
// Checks the content of the form for errors and reports any it finds to the user.
//
function check() {
  //check the motif
  if (document.getElementById('type_m').checked) {
    var motif_file = document.getElementById('motif_m').value;
    if (motif_file == null || motif_file.match(/^\s*$/)) {
      alert("Please enter a meme motif file for upload.");
      return false;
    }
    if (!document.getElementById('m_all').checked) {
      var motif_filter = document.getElementById('m_filter').value;
      if (!motif_filter.match(/^\s*$/)) {
        var names = motif_filter.split(/\s+/);
        for (var i = 0; i < names.length; ++i) {
          if (!names[i].match(/^[a-zA-Z0-9:_\.][a-zA-Z0-9:_\.-]*$/)) {
            alert("A name in the motif filter contains unsupported characters." + 
                " Motif names can not begin with a dash and may only have " +
                "alphanumeric characters plus ':', '_', '.' and '-'.");
            return false;
          }
        }
      }
    }
  } else {
    update_preview();
    if (logo_error_count > 0) {
      alert("Errors were found while checking the motif. Please read the error log for more detail.");
      return false;
    }
    if (logo_warning_count > 0) {
      var ignore_warnings = confirm("Warnings were generated while checking the motif. Do you want to continue anyway?");
      if (!ignore_warnings) {
        preview_tab(false);
        return false;
      }
    }
  }
  //check the database
  var database = document.getElementById('target_db').value;
  var db_file = document.getElementById("db_file").value;
  if (database == null || database.match(/^\s*$/)) {
    alert("Please select a database to search.");
    return false;
  } else if (database === "upload" && db_file == "") {
    alert("Please select a motif database file to upload.");
    return false;
  }
  //check the email
  if (need_email()) {
    var email = document.getElementById('email').value;
    var email_confirm = document.getElementById('email_confirm').value;
    if (email.indexOf('@') == -1) {
      alert("Please enter an email address.");
      return false;
    } else if (email != email_confirm) {
      alert("Please correct your email address as it doesn't match the confirmation address.");
      return false;
    }
  }
  //check the threshold
  var ttype = document.getElementById('thresh_type').value;
  var threshold = parseNum(document.getElementById('thresh').value);
  if (ttype == "qvalue" && (isNaN(threshold) || threshold <= 0 || threshold >= 1)) {
    alert("Please enter a number between 0 and 1 for the q-value threshold.");
    return false;
  } else if (ttype == "evalue" && (isNaN(threshold) || threshold <= 0)) {
    alert("Please enter a number larger than zero for the E-value threshold.");
    return false;
  } 
  return true;
}

//
// motif_format
//
// Sets the data entry items visible for the selected format.
//
function motif_format(format) {
  document.getElementById('consensus').style.display = 'none';
  document.getElementById('matrix').style.display = 'none';
  document.getElementById('meme').style.display = 'none';
  // check the inline one only if it's avaliable
  if (document.getElementById('type_i') != null){
    document.getElementById('inline').style.display = 'none';
  }
  var preview = document.getElementById('preview_div');
  if (format == 'meme') preview.style.display = 'none';
  else preview.style.display = 'block';
  document.getElementById(format).style.display = 'block'; 
  update_email_visibility();
  update_preview();
}

//
// preview_tab
//
function preview_tab(is_logo) {
  if (is_logo) {
    document.getElementById('preview_canvas_tab').className = 'tab active_tab';
    document.getElementById('preview_log_tab').className = 'tab';
    show_hide('preview_canvas_div', 'preview_log_div');
  } else {
    document.getElementById('preview_canvas_tab').className = 'tab';
    document.getElementById('preview_log_tab').className = 'tab active_tab';
    show_hide('preview_log_div', 'preview_canvas_div');
  }
}

function show_db_upload() {
  var value = document.getElementById('target_db').value;
  
  if (value == 'upload') {
    show_hide('db_upload_div', '');  
    show_hide('queue_div', '');
  } else {
    show_hide('', 'db_upload_div');  
    if (!need_email()) show_hide('', 'queue_div');
  }
}

//
// show_hide
//
// Generic show/hide function for making toggleable sections.
//
function show_hide(show_id, hide_id) {
  if (hide_id != '') document.getElementById(hide_id).style.display = 'none';
  if (show_id != '') document.getElementById(show_id).style.display = 'block';
}

//
// toggle
//
// Generic toggle visibility of a previously hidden element.
//
function show_toggle(id) {
  var vis = vis_ids[id];
  if (vis) {
    document.getElementById(id).style.display = 'none';
    delete vis_ids[id];
  } else {
    document.getElementById(id).style.display = 'block';
    vis_ids[id] = true;
  }
}

//
// example
//
// Sets the count matrix to the example motif.
//
function example_col_matrix() {
  document.getElementById('motif_f').value =  
      ' 3  3 19  0  1  0  2 26  5\n' +
      ' 8  0  0  1  0  1 23  1 15\n' +
      '14  0  9 27 26  4  3  0  4\n' + 
      ' 3 25  0  0  1 23  0  1  4';
  update_preview();
}

function example_row_matrix() {
  document.getElementById('motif_f').value =  
      ' 3  8 14  3\n' +
      ' 3  0  0 25\n' +
      '19  0  9  0\n' +
      ' 0  1 27  0\n' +
      ' 1  0 26  1\n' +
      ' 0  1  4 23\n' +
      ' 2 23  3  0\n' +
      '26  1  0  1\n' +
      ' 5 15  4  4\n';
  update_preview();
}
//
// addMsg
//
// Adds a message to the log box
function addMsg(msg, header, indent) {
  if (header === undefined) header = '';
  if (indent === undefined) {
    if (header == '') indent = '';
    else indent = '  ';
  }
  var out = '';
  var lines = msg.split(/\n/);
  out = header + lines[0];
  for (var i = 1; i < lines.length; ++i) {
    out += "\n" + indent + lines[i];
  }
  if (log_box == null) return;
  if (log_box.value == '') {
    log_box.value = out;
  } else {
    log_box.value += "\n" + out;
  }
}

//
// addError
//
// Adds an error to the log box.
//
function addError(msg) {
  addMsg(msg, "Bad Input: ");
  logo_error_count++;
}

//
// addWarning
//
// Adds a warning to the log box.
//
function addWarning(msg) {
  addMsg(msg, "Warning: ");
  logo_warning_count++;
}

//
// parseNum
//
// Parses a floating point number not allowing for anything else but white space.
//
function parseNum(text) {
  if (text.match(/^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$/)) {
    return parseFloat(text);
  } else {
    return Number.NaN;
  }
}


//
// getAlphabetLetter
//
// Get the text for the letter's bg and convert it to a number
// substituting defaults when errors occur.
//
function getAlphabetLetter(base, letter) {
  var text = document.getElementById(base + letter).value;
  var bg = parseNum(text);
  if (isNaN(bg)) {
    addWarning("Background " + letter + " is not a number.");
    bg = 0.25;
  } else if (bg < 0) {
    addWarning("Background " + letter + " is < 0.");
    bg = 0.25;
  }
  return bg;
}

//
// getAlphabet
//
// Get the bg for each of the alphabet's letters
// and normalise it into a 0 order markov model.
//
function getAlphabet(base) {
  var bgA = getAlphabetLetter(base, 'A'); // provide default
  var bgC = getAlphabetLetter(base, 'C'); // provide default
  var bgG = getAlphabetLetter(base, 'G'); // provide default
  var bgT = getAlphabetLetter(base, 'T'); // provide default
  var bgSum = bgA + bgC + bgG + bgT;
  bgA /= bgSum;
  bgC /= bgSum;
  bgG /= bgSum;
  bgT /= bgSum;
  return new Alphabet("ACGT", "A " + bgA + " C " + bgC + " G " + bgG + " T " + bgT);  
}

//
// getSites
//
// Get the sites.
function getSites(id) {
  var text = document.getElementById(id).value;
  var num = parseNum(text);
  var sites = num;
  if (isNaN(num)) {
    addWarning("Sites is not a number.");
    sites = 10;
  } else if (Math.round(num) != num) {
    addWarning("Sites is a decimal number.");
    sites = Math.round(num);
  }
  if (sites < 1) {
    addWarning("Sites is < 1.");
    sites = 10;
  }
  return sites;
}

//
// getPseudocount
//
// Get the pseudocount
//
function getPseudocount(id) {
  var text = document.getElementById(id).value;
  var num = parseNum(text);
  var count = num;
  if (isNaN(num)) {
    addWarning("Pseudocount is not a number.");
    count = 1;
  } else if (num < 0) {
    addWarning("Pseudocount is < 0.");
    count = 1;
  }
  return count;
}

function seq_parse_error(seq, error_pos, error_msg) {
  var msg = "";
  var offset = 0;
  if (error_pos.length > 0) {
    msg += "\n" + seq + "\n";
    for (var i = 0; i < error_pos.length; ++i) {
      if (error_pos[i] < offset) continue;
      for (; offset < error_pos[i]; ++offset) msg += " ";
      msg += "^";
      offset++;
    }
    msg += " " + error_msg;
  }
  return msg;
}

//
// consensus_logo
//
// Creates a logo object from the consensus data
//
function consensus_logo() {
  var sites = getSites('s_sites');
  var pseudocount = getPseudocount('s_pseudo')
  var alpha = getAlphabet('s_bg_');
  var seq = document.getElementById('motif_s').value;
  if (seq.length == 0) {
    addError("No IUPAC motif.");
    return null;
  }
  //==========================================================================
  // convert the consensus to a series of numbers
  //==========================================================================
  var nums = [];
  {
    //status variables
    var in_bracket = false;
    var bracket_value;
    var bracket_offset;
    //problems
    var unrecognized = [];
    var unexpected_open = [];
    var unexpected_close = [];
    var nocontent = [];
    
    for (var i = 0; i < seq.length; i++) {
      var letter = seq.charAt(i);
      if (!in_bracket) {
        if (letter == '[') {
          in_bracket = true;
          bracket_value = 0; //binary 0000
          bracket_offset = i;
        } else {
          var num = codes[letter.toUpperCase()];
          if (num != undefined) {
            nums.push(num);
          } else if (letter == ']') {
            unexpected_close.push(i);
          } else {
            unrecognized.push(i);
          }
        }
      } else {
        if (letter == ']') {
          in_bracket = false;
          if (bracket_value != 0) {
            nums.push(bracket_value);
          } else {
            nocontent.push(bracket_offset);
          }
        } else {
          var num = codes[letter.toUpperCase()];
          if (num != undefined) {
            bracket_value |= num;
          } else if (letter == '[') {
            unexpected_open.push(i);
          } else {
            unrecognized.push(i);
          }
        }
      }
    }
    if (in_bracket) { //auto close brackets
      if (bracket_value != 0) {
        nums.push(bracket_value);
      } else {
        nocontent.push(bracket_offset);
      }
    }
    var problem_count = unrecognized.length + unexpected_open.length + unexpected_close.length + nocontent.length + (in_bracket ? 1 : 0);
    if (problem_count > 0) {
      var msg = "" + problem_count + " IUPAC problem" + (problem_count > 1 ? "s" : "")
      msg += seq_parse_error(seq, unrecognized, "Unrecognized");
      msg += seq_parse_error(seq, unexpected_open, "Unexpected '['");
      msg += seq_parse_error(seq, unexpected_close, "Unexpected ']'");
      if (in_bracket) {
        msg += "\n" + seq + "\n";
        for (offset = 0; offset < bracket_offset; ++offset) msg += " ";
        msg += "^ Missing ']'";
      }
      msg += seq_parse_error(seq, nocontent, "No content in [...]");
      addWarning(msg);
    }
  }
  //==========================================================================
  // check for content to continue with
  //==========================================================================
  if (nums.length == 0) {
    addError("No recognised IUPAC codes.");
    return null;
  }
  //==========================================================================
  // convert the numbers into a pspm
  //==========================================================================
  var pspm_text = "letter-probability matrix: alength= 4 w= " + nums.length + " nsites= " + sites + " E= 0";
  {
    for (var i = 0; i < nums.length; ++i) {
      var freqs = [0, 0, 0, 0];
      var num = nums[i];
      //count how many letters are set to get the divisor for sites
      // at the same time distribute the pseudocount
      var count = 0;
      for (var a = 0; a < 4; ++a) {
        if (1<<a & num) ++count;
        freqs[a] = alpha.get_bg_freq(a) * pseudocount;
      }
      var site_part = sites / count;
      for (var a = 0; a < 4; ++a) {
        if (1<<a & num) freqs[a] += site_part;
        //normalise
        freqs[a] /= (sites + pseudocount);
      }
      // create row
      pspm_text += "\n" + freqs.join(" ");
    }
  }
  // create a Pspm object
  var pspm = new Pspm(pspm_text, "query");
  // return the logo
  return logo_1(alpha, "MEME", pspm) 
}


//
// matrix_transpose
//
// Take an array representing a matrix and it's dimensions and 
// return an array representing the transpose.
//
// There's probably some way of doing this in place but 
// (with the exception of square matricies) I can't think 
// of it so I'm making a copy.
//
// So a matrix like this:
//
// a b c d
// e f g h
//
// Would become like this:
//
// a e
// b f
// c g
// d h
function matrix_transpose(matrix_array, matrix_width, matrix_height) {
  // test that the dimensions we were given make sense
  if (matrix_array.length != matrix_width * matrix_height) throw "BAD_DIMENSIONS";
  // apparently this is a slow way of making an array
  var copy_array = new Array(matrix_width * matrix_height);
  var copy_width = matrix_height;
  var copy_height = matrix_width;
  for (var x = 0; x < matrix_width; ++x) {
    for (var y = 0; y < matrix_height; ++y) {
      var copy_x = y;
      var copy_y = x;
      copy_array[copy_x + copy_y * copy_width] = matrix_array[x + y * matrix_width];
    }
  }
  return copy_array;
}

//
// matrix_preview
// 
// Creates a motif preview of a count or probability matrix
//
function matrix_logo() {
  var orientation = document.getElementById('f_orient').value;
  var sites = getSites('f_sites');
  var pseudocount = getPseudocount('f_pseudo');
  var alpha = getAlphabet('f_bg_');
  var matrix = document.getElementById('motif_f').value;
  //==========================================================================
  // convert the matrix into a pspm
  //==========================================================================
  var pspm_text = '';
  {
    var matrix_array = [];
    var lines = matrix.split(/\n/);
    var height = 0;
    var width = (orientation == "row" ? 4 : -1);
    var total = 0;
    //=========================================================================
    // read the matrix into an array
    //=========================================================================
    for (var line_i = 0; line_i < lines.length; ++line_i) {
      var line = lines[line_i].replace(/^\s+|\s+$/g, '');//trim
      //skip empty lines
      if (line.match(/^\s*$/)) continue;
      //break into numbers
      var nums = line.split(/\s+/);
      if (width == -1) {
        // use the first width as the expected width
        width = nums.length;
      } else if (nums.length != width) { 
        // the dimensions don't match the expected so error
        if (height == 0) { 
          addError("A row oriented matrix must have 4 entries per row.");
        } else {
          addError("Row " + (height+1) + " of the matrix does not have " + width + " entries as expected from the previous rows.");
        }
        return null;
      }
      // count this row
      height++;
      for (var num_i = 0; num_i < nums.length; ++num_i) {
        var num = parseNum(nums[num_i]);
        if (isNaN(num)) {
          addWarning("Matrix value at row: " + (height) + " column: " + (num_i + 1) + " is not a number. Substituting 0.");
          num = 0;
        } else if (num < 0) {
          addWarning("Matrix value at row: " + (height) + " column: " + (num_i + 1) + " is < 0. Substituting 0.");
          num = 0;
        }
        matrix_array.push(num);
        total += num;
      }
    }
    // check the dimensions
    if (orientation == "col" && height != 4) {
      addError("A column oriented matrix must have 4 entries per column.");
      return null;
    }
    if (width == -1) {
      addError("No matrix.");
      return null;
    }
    if (orientation == "auto" && width != 4 && height != 4) {
      addError("Matrix width or height must be 4.");
      return null;
    }
    //=========================================================================
    // correct the orientation of the matrix
    //=========================================================================
    var transposed = false;
    if (orientation == "auto") {
      if (width == 4 && height == 4) {
        //calculate the variance of rows and cols and choose the smaller
        var avg = total / 4;
        var row_variance = 0;
        var col_variance = 0;
        for (var i = 0; i < 4; ++i) {
          var row_sum = 0;
          var col_sum = 0;
          for (var j = 0; j < 4; ++j) {
            row_sum += matrix_array[i*width + j];
            col_sum += matrix_array[j*width + i];
          }
          row_variance += Math.pow(row_sum - avg, 2);
          col_variance += Math.pow(col_sum - avg, 2);
        }
        row_variance /= 4;
        col_variance /= 4;
        // if one has a variance of 0 then it's not really a guess
        if (row_variance == col_variance) {
          addWarning("No difference in variance of summed rows or summed columns. Assuming a row matrix.");
        } else if (row_variance != 0 && col_variance != 0) {
          addWarning("Neither rows or columns sum to a constant. Guessing orientation based on lowest variance.");
        }
        if (row_variance > col_variance) { //guess col matrix as it has lower variance
          //transpose
          matrix_array = matrix_transpose(matrix_array, width, height);
          transposed = true;
        }
      } else if (height == 4) { // must be column matrix
        //transpose
        matrix_array = matrix_transpose(matrix_array, width, height);
        height = width;
        width = 4;
        transposed = true;
      }
    } else if (orientation == "col") {
      //transpose
      matrix_array = matrix_transpose(matrix_array, width, height);
      height = width;
      width = 4;
      transposed = true;
    }
    //=========================================================================
    // convert into a position specific probability matrix
    //=========================================================================
    // we should have a row matrix now
    // calculate the average row sum
    var row_avg = total / height;
    // calculate the actual row sums and the row variance
    var row_sums = [];
    var variance = 0;
    var max_diff = 0;
    //for each row
    for (var y = 0; y < height; ++y) {
      //sum the row
      var row_sum = 0;
      for (var x = 0; x < 4; ++x) {
        row_sum += matrix_array[y*width + x];
      }
      if (row_sum == 0) {
        if (transposed) {
          addError("Column " + y + " in the column matrix summed to zero indicating that all nucleotides had zero probability!");
        } else {
          addError("Row " + y + " in the row matrix summed to zero indicating that all nucleotides had zero probability!");
        }
        return null; //can't preview
      }      
      var diff = row_sum - row_avg;
      max_diff = Math.max(max_diff, Math.abs(diff));
      row_sums.push(row_sum);
      variance += Math.pow(diff, 2);
    }
    variance /= height;

    // check the variance
    if (max_diff > 0.1) {
      var msg = "Not all the " + (transposed ? "columns" : "rows") + " summed to the same value. Sums are:\n" + row_sums.join(" ");
      addWarning(msg);
    }

    // detect probability matrix
    if (Math.round(row_avg) > 1) { //it's probably a count matrix so use row_avg as sites
      sites = Math.round(row_avg);
      if (sites != row_avg) {
        addWarning("The number of sites as derived from the count matrix was not a whole number. Rounding.");
      }
    }

    pspm_text = "letter-probability matrix: alength= 4 w= " + height + " nsites= " + sites + " E= 0";
    // normalise each row individually
    for (var y = 0; y < height; ++y) {
      //normalise
      for (var x = 0; x < width; ++x) matrix_array[y*width + x] /= row_sums[y];
      //convert to counts
      for (var x = 0; x < width; ++x) matrix_array[y*width + x] *= sites;
      //add pseudocounts
      for (var x = 0; x < width; ++x) matrix_array[y*width + x] += (alpha.get_bg_freq(x) * pseudocount);
      //renormalise
      for (var x = 0; x < width; ++x) matrix_array[y*width + x] /= (sites + pseudocount);
      //output
      pspm_text += "\n" + matrix_array[y*width] + " " + matrix_array[y*width + 1] + " " + matrix_array[y*width + 2] + " " + matrix_array[y*width + 3];
    }
  }

  // create a Pspm object
  var pspm = new Pspm(pspm_text, "query");
  // draw the logo
  return logo_1(alpha, "MEME", pspm);
}

//
// inline_logo
//
// Creates a logo of a inline motif (in minimal meme format)
//
function inline_logo() {
  var motifs = document.getElementById('motif_i').value;
  var bg;
  var pspm_text;
  // now search the motif for the background and the first pspm
  var line_i = 0;
  var lines = motifs.split(/\n/);
  // find the background
  for (; line_i < lines.length; line_i++) {
    if (lines[line_i].match(/Background letter frequencies/)) {
      bg = lines[line_i+1].replace(/^\s+|\s+$/g, '');//trim
      break;
    }
  }
  // find the first pspm
  var pspm_regex = /letter-probability matrix: alength= 4 w= (\d+)/;
  for (; line_i < lines.length; line_i++) {
    var hit = pspm_regex.exec(lines[line_i])
    if (hit) {
      var width = hit[1];
      pspm_text = lines[line_i];
      for (var i = 1; i <= width; i++) {
        pspm_text += "\n" + lines[line_i + i];
      }
      break;
    }
  }
  var alphabet = new Alphabet("ACGT", bg);
  var pspm = new Pspm(pspm_text, "query");
  // draw the logo
  return logo_1(alphabet, "MEME", pspm);
}


//
// update_preview
//
// Creates a preview of the motif for the sequence, matrix and inline formats. 
// It can't do meme because it never has access to the data (unless I do something 
// very fancy with ajax)
//
function update_preview() {
  // ensure that update_preview isn't called twice for no reason
  if (preview_delay_timer != null) clearTimeout(preview_delay_timer);
  // get the preview canvas
  var canvas = document.getElementById('motif_preview');
  // clear the canvas
  canvas.width = canvas.width;
  // clear the log
  if (log_box != null) log_box.value = '';
  // reset the error count
  logo_error_count = 0;
  logo_warning_count = 0;
  // check for canvas logo support
  var preview_support = false;
  if (canvas.getContext) {
    var ctx = canvas.getContext('2d');
    if (supports_text(ctx)) {
      preview_support = true;
    }
  }
  // disable preview if unsupported
  if (!preview_support){
    addMsg("Info: Preview unsupported.");
    document.getElementById('preview_canvas_tab').style.display = 'none';
  }
  // attempt to get a new logo
  var logo = null;
  if (document.getElementById('type_s').checked) {
    logo = consensus_logo();
  } else if (document.getElementById('type_f').checked) {
    logo = matrix_logo();
  } else if (document.getElementById('type_m').checked) {
    // no way to preview meme files
  } else {
    var inline_motif = document.getElementById('type_i');
    if (inline_motif != null && inline_motif.checked) {
      logo = inline_logo();
    }
  }
  // give a success message
  if (logo_warning_count == 0 && logo_error_count == 0) addMsg("Info: Motif ok.");
  // draw the preview
  if (logo != null) {
    if (preview_support) draw_logo_on_canvas(logo, canvas);
  }
  // switch to the log on error
  if (logo_error_count > 0) {
    if (preview_support) document.getElementById('preview_canvas_tab').style.display = 'none';
    preview_tab(false);
    log_box.style.background = "#FFE4E1";
  } else {
    if (preview_support) document.getElementById('preview_canvas_tab').style.display = 'inline';
    if (logo_warning_count > 0) {
      log_box.style.background = "#EEE8AA";
    } else  {
      log_box.style.background = "#E4ECEC";
    }
  }
  // switch to the preview when no errors or warnings
  if (preview_support && logo_warning_count == 0 && logo_error_count == 0) preview_tab(true);
}

//
// delayed_preview
//
// Waits a short delay for more changes and then updates the preview.
//
function delayed_preview() {
  if (preview_delay_timer != null) clearTimeout(preview_delay_timer);
  preview_delay_timer = setTimeout("update_preview()", 500);
}

//
// need_email
//
// Returns true if an email is needed.
//
function need_email() {
  if (document.getElementById('type_m').checked) {
    if (document.getElementById('m_all').checked) {
      return true;
    } else if (document.getElementById('m_filter').value.match(/\S\s+\S/)) {
      return true;
    }
  } else {
    var inline_radio = document.getElementById('type_i');
    if (inline_radio != null && inline_radio.checked) {
      if (document.getElementById('i_count').value > 1) return true;
    }
  }
  if (document.getElementById('target_db').value == 'upload') {
    return true;
  }
  return false;
}

//
// update_email_visibility
//
// If the email is required then make the field visible
//
function update_email_visibility() {
  if (need_email()) {
    show_hide('queue_div', '');
  } else {
    show_hide('', 'queue_div');
  }
}
