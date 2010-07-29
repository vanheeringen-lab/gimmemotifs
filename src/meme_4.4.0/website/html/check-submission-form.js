/*
 * $Id: $
 *
 * $Log: $
 *
 */
function check(form) {
        var email = form.address.value;
        var email_verify = form.address_verify.value;
        var description = form.description.value;

        // check address field
        if(!email) {
                alert("Please fill in the 'e-mail address' field.");
                return false;
        }
        if(!email_verify) {
                alert("Please fill in the 'Re-enter email address' field.");
                return false;
        }
        if(!email.match(email_verify)) {
                alert("Email addresses do not match.");
                return false;
        }
        if(!email_verify.match(email)) {
                alert("Email addresses do not match.");
                return false;
        }

        // check description field; reg exp must match one in Utils.pm.in
        if(description.match(/[^\w _\-\(\)]/)) {
                alert("You may only use letters, numbers, underlines, dashes, blanks and  parentheses in the description field.");
                return false;
        }
}
