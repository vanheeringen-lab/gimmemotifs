      //
      // Functions to measure the position of page elements relative to the size of the page
      //

      // gets the offset of the top of the page due to scrolling
      // from: http://www.howtocreate.co.uk/tutorials/javascript/browserwindow
      function get_scroll_xy() {
        var scrOfX = 0, scrOfY = 0;
        if( typeof( window.pageYOffset ) == 'number' ) {
          //Netscape compliant
          scrOfY = window.pageYOffset;
          scrOfX = window.pageXOffset;
        } else if( document.body && ( document.body.scrollLeft || document.body.scrollTop ) ) {
          //DOM compliant
          scrOfY = document.body.scrollTop;
          scrOfX = document.body.scrollLeft;
        } else if( document.documentElement && ( document.documentElement.scrollLeft || document.documentElement.scrollTop ) ) {
          //IE6 standards compliant mode
          scrOfY = document.documentElement.scrollTop;
          scrOfX = document.documentElement.scrollLeft;
        }
        return [ scrOfX, scrOfY ];
      }

      // gets the width and height of the visible page
      // from: http://www.howtocreate.co.uk/tutorials/javascript/browserwindow
      function get_page_size() {
        var myWidth = 0, myHeight = 0;
        if( typeof( window.innerWidth ) == 'number' ) {
          //Non-IE
          myWidth = window.innerWidth;
          myHeight = window.innerHeight;
        } else if( document.documentElement && ( document.documentElement.clientWidth || document.documentElement.clientHeight ) ) {
          //IE 6+ in 'standards compliant mode'
          myWidth = document.documentElement.clientWidth;
          myHeight = document.documentElement.clientHeight;
        } else if( document.body && ( document.body.clientWidth || document.body.clientHeight ) ) {
          //IE 4 compatible
          myWidth = document.body.clientWidth;
          myHeight = document.body.clientHeight;
        }
        return [myWidth, myHeight];
      }

      // gets the x and y offset of an element
      // from http://www.quirksmode.org/js/findpos.html
      function get_elem_xy(elem) {
        var myX = myY = 0;
        if (elem.offsetParent) {
          do {
            myX += elem.offsetLeft;
            myY += elem.offsetTop;
          } while (elem = elem.offsetParent);
        }
        return [myX, myY];
      }

      //
      // Functions to delay a drawing task until it is required or it would not lag the display to do so
      //

      // a list of items still to be drawn
      var drawable_list = [];
      // the delay between drawing objects that are not currently visible
      var draw_delay = 1;
      // the delay after a user interaction
      var user_delay = 300;
      // the delay after a user has stopped scrolling and is viewing the stuff drawn on the current page
      var stop_delay = 2000;
      // the timer handle; allows resetting of the timer after user interactions
      var draw_timer = null;

      //
      // Drawable
      //
      // elem - a page element which defines the position on the page that drawing is to be done
      // task - an object with the method run which takes care of painting the object
      //
      function Drawable(elem, task) {
        this.elem = elem;
        this.task = task;
        this.is_visible = Drawable_is_visible;
      }

      //
      // Drawable_is_visible
      //
      // Determines if the drawable object is in the visible part of the page.
      //
      // page_top - the distance to the top of the page for the visible region
      // page_height - the height of the visible region
      function Drawable_is_visible(page_top, page_height) {
        var elem_top = get_elem_xy(this.elem)[1];
        var elem_height = this.elem.height;
        if (typeof (elem_height) != 'number') elem_height = 1;
        return ((elem_top + elem_height) >= page_top && elem_top <= (page_top + page_height));
      }

      //
      // draw_on_screen
      //
      // Checks each drawable object and draws those on screen.
      //
      function draw_on_screen() {
        var found = false;
        var page_top = get_scroll_xy()[1];
        var page_height = get_page_size()[1];
        for (var i = 0; i < drawable_list.length; i++) {
          var drawable = drawable_list[i];
          if (drawable.is_visible(page_top, page_height)) {
            drawable.task.run();
            drawable_list.splice(i--, 1);
            found = true;
          }
        }
        return found;
      }

      //
      // process_draw_tasks
      //
      // Called on a delay to process the next avaliable
      // draw task.
      //
      function process_draw_tasks() {
        var delay = draw_delay;
        draw_timer = null;
        if (drawable_list.length == 0) return; //no more tasks
        if (draw_on_screen()) {
          delay = stop_delay; //give the user a chance to scroll
        } else {
          //get next task
          var drawable = drawable_list.shift();
          drawable.task.run();
        }
        //allow UI updates between tasks
        draw_timer = window.setTimeout("process_draw_tasks()", draw_delay);
      }

      //
      // delayed_process_draw_tasks
      //
      // Call process_draw_tasks after a short delay.
      // The delay serves to group multiple redundant events.       
      // Should be set as event handler for onscroll and onresize.
      //
      function delayed_process_draw_tasks() {
        //reset the timer
        if (drawable_list.length > 0) { 
          if (draw_timer != null) clearTimeout(draw_timer);
          draw_timer = window.setTimeout("process_draw_tasks()", user_delay);
        }
      }

      //
      // add_draw_task
      //
      // Add a drawing task to be called immediately if it is
      // visible, or to be called on a delay to reduce stuttering
      // effect on the web browser.
      function add_draw_task(elem, task) {
        var page_top = get_scroll_xy()[1];
        var page_height = get_page_size()[1];
        drawable = new Drawable(elem, task);
        if (drawable.is_visible(page_top, page_height)) {
          task.run();
        } else {
          drawable_list.push(drawable);
          //reset timer
          if (draw_timer != null) clearTimeout(draw_timer);
          draw_timer = window.setTimeout("process_draw_tasks()", user_delay);
        }
      }

