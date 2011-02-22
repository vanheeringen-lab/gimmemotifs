
var fasta_db_url;

function check(form) {

  var method = form.search_methods
  var alphabet = form.alphabet;
  var motif = form.motifs;
  var sequence_source = form.sequence_source
  var uploaddb = form.upload_db;
  var database = form.database;
  var email = form.address;

  // If GLAM is checked they need to specify the alphabet
  if (method[2].checked) {
    if (! alphabet.value || alphabet.value == ',') {
      alert("Please select a sequence alphabet.");
      alphabet.focus();
      return false;
    }
  }
  if (!motif.value) {
    alert("Please select a motif file to upload.");
    motif.focus();
    return false;
  }
  // If sequence source is upload
  if (sequence_source[0].checked) {
    if (! uploaddb.value) {
      alert("Please select a sequence database to upload.");
      uploaddb.focus();
      return false;
    }
  }
  // If sequence source is supported database
  else if (sequence_source[1].checked) {
    if (! database.value) {
      alert("Please select a sequence database to search.");
      database.focus();
      return false;
    }
  }
  else {
   alert("Please select a sequence database source.");
  }
  if(! email.value) {
    alert("Please fill in the \"e-mail address\" field.");
    email.focus();
    return false;
  }
  return true;

}

function isSearchMethodChecked() {
  for (i = 0; i < document.forms[0].search_methods.length; i++) {
    if (document.forms[0].search_methods[i].checked) {
      return true;
    }
  }
  return false;
}

function isSequenceDatabaseSet() {
  if ((document.forms[0].sequence_source[0].checked && document.forms[0].upload_db.value) 
    || (document.forms[0].sequence_source[1].checked && document.forms[0].database.value != "")){
    return true;
  }
  else {
    return false;
  }
}

function getStage() {

  var stage = 1;
  if(isSearchMethodChecked()) {
    stage = 2;
  }
  if (document.forms[0].motifs.value) {
    stage = 3;
  }
  if (isSequenceDatabaseSet()) {
    stage = 4;
  }
  if (document.forms[0].address.value) {
    stage = 5;
  }

  return stage;
}

function setStage(i) {
  switch(i) {
    case 1:
      document.getElementById('submit_button').style.display = "None";
      document.getElementById('step1').style.display = "list-item";
      document.getElementById('step2').style.display = "None";
      document.getElementById('step3').style.display = "None";
      document.getElementById('step4').style.display = "None";
      document.getElementById('step5').style.display = "None";
      document.getElementById('progress_step1').style.color = "red";
      document.getElementById('progress_line1').style.background = "gray";
      document.getElementById('progress_step2').style.color = "gray";
      document.getElementById('progress_line2').style.background = "gray";
      document.getElementById('progress_step3').style.color = "gray";
      document.getElementById('progress_line3').style.background= "gray";
      document.getElementById('progress_step4').style.color = "gray";
      document.getElementById('progress_line4').style.background= "gray";
      break;
    case 2:
      document.getElementById('submit_button').style.display = "None";
      document.getElementById('step1').style.display = "list-item";
      document.getElementById('step2').style.display = "list-item";
      document.getElementById('step3').style.display = "None";
      document.getElementById('step4').style.display = "None";
      document.getElementById('step5').style.display = "None";
      document.getElementById('progress_step1').style.color = "black";
      document.getElementById('progress_line1').style.background = "#02656A";
      document.getElementById('progress_step2').style.color = "red";
      document.getElementById('progress_line2').color = "gray";
      document.getElementById('progress_step3').style.color = "gray";
      document.getElementById('progress_line3').color = "gray";
      document.getElementById('progress_step4').style.color = "gray";
      break;
    case 3:
      document.getElementById('submit_button').style.display = "None";
      document.getElementById('step1').style.display = "list-item";
      document.getElementById('step2').style.display = "list-item";
      document.getElementById('step3').style.display = "list-item";
      document.getElementById('step4').style.display = "None";
      document.getElementById('step5').style.display = "None";
      document.getElementById('progress_step1').style.color = "black";
      document.getElementById('progress_line1').style.background = "#02656A";
      document.getElementById('progress_step2').style.color = "black";
      document.getElementById('progress_line2').style.background = "#02656A";
      document.getElementById('progress_step3').style.color = "red";
      document.getElementById('progress_line3').color = "gray";
      document.getElementById('progress_step4').style.color = "gray";
      if (document.forms[0].sequence_source[0].checked) {
        document.getElementById('uploaded_database').style.display = "inline";
        document.getElementById('supported_database').style.display = "None";
      }
      else if (document.forms[0].sequence_source[1].checked) {
        document.getElementById('uploaded_database').style.display = "None";
        document.getElementById('supported_database').style.display = "inline";
      }
      else {
        document.getElementById('uploaded_database').style.display = "None";
        document.getElementById('supported_database').style.display = "None";
      }
      break;
    case 4:
      document.getElementById('submit_button').style.display = "None";
      document.getElementById('step1').style.display = "list-item";
      document.getElementById('step2').style.display = "list-item";
      document.getElementById('step3').style.display = "list-item";
      document.getElementById('step4').style.display = "list-item";
      document.getElementById('step5').style.display = "list-item";
      document.getElementById('progress_step1').style.color = "black";
      document.getElementById('progress_line1').style.background = "#02656A";
      document.getElementById('progress_step2').style.color = "black";
      document.getElementById('progress_line2').style.background = "#02656A";
      document.getElementById('progress_step3').style.color = "black";
      document.getElementById('progress_line3').style.background = "#02656A";
      document.getElementById('progress_step4').style.color = "black";
      document.getElementById('progress_line4').style.background = "#02656A";
      document.getElementById('progress_step5').style.color = "red";
      if (document.forms[0].sequence_source[0].checked) {
        document.getElementById('uploaded_database').style.display = "inline";
        document.getElementById('supported_database').style.display = "None";
      }
      else if (document.forms[0].sequence_source[1].checked) {
        document.getElementById('uploaded_database').style.display = "None";
        document.getElementById('supported_database').style.display = "inline";
      }
      else {
        document.getElementById('uploaded_database').style.display = "inline";
        document.getElementById('supported_database').style.display = "None";
      }
      break;
    case 5:
      document.getElementById('submit_button').style.display = "inline";
      document.getElementById('step1').style.display = "list-item";
      document.getElementById('step2').style.display = "list-item";
      document.getElementById('step3').style.display = "list-item";
      document.getElementById('step4').style.display = "list-item";
      document.getElementById('step5').style.display = "list-item";
      document.getElementById('progress_step1').style.color = "black";
      document.getElementById('progress_line1').style.background = "#02656A";
      document.getElementById('progress_step2').style.color = "black";
      document.getElementById('progress_line2').style.background = "#02656A";
      document.getElementById('progress_step3').style.color = "black";
      document.getElementById('progress_line3').style.background = "#02656A";
      document.getElementById('progress_step4').style.color = "black";
      document.getElementById('progress_line4').style.background = "#02656A";
      document.getElementById('progress_step5').style.color = "black";
      if (document.forms[0].sequence_source[0].checked) {
        document.getElementById('uploaded_database').style.display = "inline";
        document.getElementById('supported_database').style.display = "None";
      }
      else if (document.forms[0].sequence_source[1].checked) {
        document.getElementById('uploaded_database').style.display = "None";
        document.getElementById('supported_database').style.display = "inline";
      }
      else {
        document.getElementById('uploaded_database').style.display = "inline";
        document.getElementById('supported_database').style.display = "None";
      }
      break;
  }
}

function OnLoad(url) {
  fasta_db_url = url;
  setStage(getStage());
  ChooseSearchMethod();
  ChangeMotifFile();
  ChooseSequenceSource();
}

function ChangeEmail() {

  setStage(5);
  return true;
}

function ChooseSearchMethod() {

  // List of search methods
  var searchActions = [ 'cgi-bin/fimo.cgi', 'cgi-bin/mast.cgi', 'cgi-bin/glam2scan.cgi'];
  var action_field = document.getElementById('target_action');

  // Only display content relevent to the selected search type 
  for (i = 0; i < document.forms[0].search_methods.length; i++) {
    var search_method = document.forms[0].search_methods[i];
    var search_parameters_div = document.getElementById(search_method.value + '_parameters');
    var motif_text_div = document.getElementById(search_method.value + '_motif');
    var motif_sample_div = document.getElementById(search_method.value + '_motif_sample');
    var db_doc_div = document.getElementById(search_method.value + '_db_doc');
    if (search_method.checked) {
      search_parameters_div.style.display = "block";
      motif_text_div.style.display = "inline";
      motif_sample_div.style.display = "inline";
      db_doc_div.style.display = "inline";
      document.forms[0].action = searchActions[i];
      action_field.value = searchActions[i];
      // MAST is short seqs only
      do_catlist_update(i==1);
      if (getStage() < 3) {
        setStage(2);
      }
    }
    else {
     search_parameters_div.style.display = "None";
     motif_text_div.style.display = "None";
     motif_sample_div.style.display = "None";
     db_doc_div.style.display = "None";
    }
  }
}

function ChangeMotifFile() {
  if (document.forms[0].motifs.value) {
    if (getStage() < 4) {
      setStage(3);
    }
  }
}

function ChooseSequenceSource(ctrl) {
  if (ctrl) {
    if (ctrl.value == "supported_database") {
      document.getElementById('supported_database').style.display = "inline";
      document.getElementById('uploaded_database').style.display = "None";
      document.getElementById('upload_db').value = "";
    }
    else {
      document.getElementById('supported_database').style.display = "None";
      document.getElementById('uploaded_database').style.display = "inline";
      document.getElementById('database').value = "";
    }
  }
}

function ChangeSequenceDatabase() {
  if (isSequenceDatabaseSet()) {
    if (getStage() < 5) setStage(4);
  } else {
    setStage(3);
  }
}

function ResetForm() {
  setStage(1);
}

function is_mast() {
  var mast_radio = document.getElementById('mast_search_method');
  if (!mast_radio) throw "MAST_RADIO_BUTTON_MISSING";
  return (mast_radio.checked);
}

function do_dblist_update(catid) {
  var url = fasta_db_url;
  url += "&mode=xml&catid=~category~";
  if (is_mast()) url += "&short_only=1";
  setStage(3); 
  send_dblist_request(url, catid, "database");//this method is in selectdb.js
}

function do_catlist_update(is_short_only) {
  var url = fasta_db_url;
  url += "&mode=categories_xml";
  if (is_short_only) url += "&short_only=1";
  send_catlist_request(url, "category"); //this method is in selectdb.js
  send_dblist_request("", "", "database"); //this method is in selectdb.js
}
