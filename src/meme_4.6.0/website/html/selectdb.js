//This file contains AJAX code to allow selection of a different category of databases

function get_request_obj() {
  if (window.XMLHttpRequest) {
    // code for IE7+, Firefox, Chrome, Opera, Safari
    return new XMLHttpRequest();
  }
  if (window.ActiveXObject) {
    // code for IE6, IE5
    return new ActiveXObject("Microsoft.XMLHTTP");
  }
  throw "NO_REQUEST_OBJ";
}

/*
 * Extracts the categories out of a request
 * and sets the passed list to contain those
 * categories and an empty item
 */
function send_catlist_request(url, list_id, default_index) {
  var list_obj = document.getElementById(list_id);
  if (list_obj === undefined) {
    alert("List id \"" + list_id + "\" not found.");
    return;
  }
  //send the request
  send_list_request(list_obj, url, default_index);
}

function send_dblist_request(url_pattern, category_index, list_id, default_index) {
  var list_obj = document.getElementById(list_id);
  if (list_obj === undefined) {
    alert("List id \"" + list_id + "\" not found.");
    return;
  }
  //exit if no category to query
  if (category_index == "") {
    list_obj.disabled = true;
    clear_list(list_obj);
    add_option(list_obj, "", "");
    return;
  }
  //create the url to query
  var url = url_pattern.replace(/~category~/i, category_index);
  //send the request
  send_list_request(list_obj, url, default_index);
}

function send_list_request(list_obj, url, default_index) {
  //clear the existing nodes
  list_obj.disabled = true;
  clear_list(list_obj);
  //add an status line
  add_option(list_obj, "", "loading");
  //now send the request
  var request_obj = get_request_obj();
  //setup the response handler
  request_obj.onreadystatechange=function() {
    switch (request_obj.readyState) {
    case 0: //not initilized
    case 1: //setup
    case 2: //sent
    case 3: //in progress
      break;
    case 4: //complete
      parse_list_response(list_obj, request_obj, default_index);
    }
  }
  //send an asynchronos request for the url
  request_obj.open("GET", url, true);
  request_obj.send(null);
}

function parse_list_response(list_obj, request_obj, default_index) {
  if (request_obj.status != 200) {
    alert("AJAX failed");
    return;
  }
  var xmlDoc = request_obj.responseXML.documentElement;
  //clear the existing nodes
  clear_list(list_obj);
  //add an empty line
  add_option(list_obj, "", "");
  //parse the list
  var items = xmlDoc.getElementsByTagName("item");
  for (var i = 0; i < items.length; i++) {
    var item = items[i];
    var item_value = item.getAttribute("index");
    var item_text = item.getAttribute("name");
    add_option(list_obj, item_value, item_text);
  }
  if (!(default_index === undefined))list_obj.selectedIndex = default_index;
  list_obj.disabled = false;
}

function clear_list(list_obj) {
  if (list_obj.hasChildNodes()) {
    while (list_obj.childNodes.length >= 1) {
      list_obj.removeChild(list_obj.lastChild);       
    } 
  }
}

function add_option(list_obj, value, text) {
  list_obj.options[list_obj.options.length] = new Option(text, value);
  // DOM method which doesn't work on IE... sigh. Left here as a warning.
  //var opt = document.createElement("option");
  //opt.value = value;
  //opt.text = text;
  //list_obj.appendChild(opt);
}
