// FILE: validate.js
// AUTHOR: Phan Lu
// CREATE DATE: 04-10-2000
// PROJECT: MHMM
// DESCRIPTION: Validate the Meta-MEME inputs at the client site to filter out
//              any bad input before it is passed to the server.  If the client
//              uses a browser that does not support JavaScript, or the
//              JavaScript is disable on that browser, CGI.pm will bypass the
//              JavaScript and pass the input to the server and let
//              submit-verify.cgi does its usual verification.


function validateForm()
{
  if (document.form1.seq.value == "")
  {
    alert('No sequence file name was entered.');
    document.form1.seq.focus();
    document.form1.seq.select();
  }
  else if (document.form1.meme.value == "")
  {
    alert('No MEME file name was entered.');
    document.form1.meme.focus();
    document.form1.meme.select();
  }
  return validateFilePath(document.form1.seq, document.form1.meme);
}


function validateFilePath(seqfileEntry, memefileEntry)
{
  var seqpath = false;
  var mempath = false;
  var seqname = seqfileEntry.value;
  var memname = memefileEntry.value;

  if ((seqname != "") && (memname != ""))
  {
    for (i = 0; i < seqname.length; i++)
    {
      if ((seqname.charAt(i) == '\\') || (seqname.charAt(i) == '/'))
      {
        seqpath = true;
        break;
      }
    }

    for (i = 0; i < memname.length; i++)
    {
      if ((memname.charAt(i) == '\\') || (memname.charAt(i) == '/'))
      {
        mempath = true;
        break;
      }
    }

    if (!seqpath)
    {
      if (!confirm("The sequence file name requires a directory path.\n" +
                   "Do you want to continue?"))
      {
        seqfileEntry.focus();
        seqfileEntry.select();
        return false;
      }
    }

    if (!mempath)
    {
      if (!confirm("The MEME file name requires a directory path.\n" +
                   "Do you want to continue?"))
      {
        memefileEntry.focus();
        memefileEntry.select();
        return false;
      }
    }
  }

  return true;
}


function changeValue(thisCheckBox)
{
  if (thisCheckBox.value == "on")
    thisCheckBox.value = "off";
  else
    thisCheckBox.value = "on";
}
