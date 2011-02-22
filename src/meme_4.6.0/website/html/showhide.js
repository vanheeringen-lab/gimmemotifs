function showhide(blockid,textid) {
    if (document.getElementById(blockid).style.display == "none") {
        document.getElementById(blockid).style.display = "block";
        document.getElementById(textid).innerHTML = "-";
    } else {
        document.getElementById(blockid).style.display = "none";
	document.getElementById(textid).innerHTML = "+";
    }
}
