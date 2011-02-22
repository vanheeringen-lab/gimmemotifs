function on_page_show() {
  if (console != undefined) console.log("on_page_show");
  //setup_motif_input();
}

function show_db_upload() {
  var value = document.getElementById('target_db').value;
  if (value == 'upload') {
    show_hide('db_upload_div', '');  
  } else {
    show_hide('', 'db_upload_div');  
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
