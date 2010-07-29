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
function send_catlist_request(url, list_id) {
  var list = document.getElementById(list_id);
  if (list === undefined) {
    alert("List id \"" + list_id + "\" not found.");
    return;
  }
  //clear the existing nodes
  list.disabled = true;
  if (list.hasChildNodes()) {
    while (list.childNodes.length >= 1) {
      list.removeChild(list.lastChild);       
    } 
  }
  //add an status line
  var emptyOpt = document.createElement("option");
  emptyOpt.value = "";
  emptyOpt.text = "loading";
  list.appendChild(emptyOpt);
  //now send the request
  url += "&mode=categories_xml";
  var req = get_request_obj();
  //setup the response handler
  req.onreadystatechange=function() {
    switch (req.readyState) {
    case 0: //not initilized
    case 1: //setup
    case 2: //sent
    case 3: //in progress
      break;
    case 4: //complete
      parse_catlist_response(list, req);
    }
  }
  //send an asynchronos request for the url
  req.open("GET", url, true);
  req.send(null);
}

function parse_catlist_response(list, req) {
  if (req.status != 200) {
    alert("AJAX failed");
    return;
  }
  var xmlDoc = req.responseXML.documentElement;
  //clear the existing nodes
  if (list.hasChildNodes()) {
    while (list.childNodes.length >= 1) {
      list.removeChild(list.lastChild);       
    } 
  }
  //add an empty line
  var emptyOpt = document.createElement("option");
  emptyOpt.value = "";
  emptyOpt.text = "";
  list.appendChild(emptyOpt);
  //fill with loaded data
  var cats = xmlDoc.getElementsByTagName("category");
  for (var i = 0; i < cats.length; i++) {
    var cat = cats[i];
    var opt = document.createElement("option");
    opt.value = "" + (i+1);
    opt.text = cat.getAttribute("name");
    list.appendChild(opt);
  }
  list.disabled = false;
}

function send_dblist_request(url, catid, list_id) {
  var list = document.getElementById(list_id);
  if (list === undefined) {
    alert("List id \"" + list_id + "\" not found.");
    return;
  }
  //clear the existing nodes
  list.disabled = true;
  if (list.hasChildNodes()) {
    while (list.childNodes.length >= 1) {
      list.removeChild(list.lastChild);       
    } 
  }
  //exit if no category to query
  if (catid == "") {
    //add an empty line
    var emptyOpt = document.createElement("option");
    emptyOpt.value = "";
    emptyOpt.text = "";
    list.appendChild(emptyOpt);
    return;
  }
  //add an status line
  var emptyOpt = document.createElement("option");
  emptyOpt.value = "";
  emptyOpt.text = "loading";
  list.appendChild(emptyOpt);
  //now send the request
  url += "&mode=xml&catid=" + catid;
  var req = get_request_obj();
  //setup the response handler
  req.onreadystatechange=function() {
    switch (req.readyState) {
    case 0: //not initilized
    case 1: //setup
    case 2: //sent
    case 3: //in progress
      break;
    case 4: //complete
      parse_dblist_response(list, req);
    }
  }
  //send an asynchronos request for the url
  req.open("GET", url, true);
  req.send(null);
}

function parse_dblist_response(list, req) {
  if (req.status != 200) {
    alert("AJAX failed");
    return;
  }
  var xmlDoc = req.responseXML.documentElement;
  //clear the existing nodes
  if (list.hasChildNodes()) {
    while (list.childNodes.length >= 1) {
      list.removeChild(list.lastChild);       
    } 
  }
  //add an empty line
  var emptyOpt = document.createElement("option");
  emptyOpt.value = "";
  emptyOpt.text = "";
  list.appendChild(emptyOpt);
  //fill with loaded data
  var dbs = xmlDoc.getElementsByTagName("db");
  for (var i = 0; i < dbs.length; i++) {
    var db = dbs[i];
    var is_prot = (db.getAttribute("prot") == "y");
    var is_nucl = (db.getAttribute("nucl") == "y");
    var is_short = (db.getAttribute("short") == "y");
    var opttext = db.getAttribute("menu");
    if (is_prot && is_nucl) {
      opttext += " (peptide and nucleotide)";
    } else if (is_prot) {
      opttext += " (peptide only)";
    } else if (is_nucl) {
      opttext += " (nucleotide only)";
    } else {
      opttext += " (unknown database type)";
    }
    var optvalue = db.getAttribute("base") + "," + 
        (is_prot ? "yes" : "no") + "," + 
        (is_nucl ? "yes" : "no") + "," + 
        (is_short ? "yes" : "no") + "," +
        db.getAttribute("menu") + "," +
        db.getAttribute("desc") + "," +
        (db.hasAttribute("start") ? db.getAttribute("start") : "") + "," +
        (db.hasAttribute("end") ? db.getAttribute("end") : "");
    var cmpdbs = db.getElementsByTagName("cmp_db");
    for (var j = 0; j < cmpdbs.length; j++) {
      var cmpdb = cmpdbs[j];
      optvalue += "," + cmpdb.getAttribute("base");
    }
    var opt = document.createElement("option");
    opt.value = optvalue;
    opt.text = opttext;
    list.appendChild(opt);
  }
  list.disabled = false;
}
