
      var alphabet = new Alphabet("ACGT");
      var motif_lookup = [];
      var graph_lookup = [];
      var group_lookup = [];

      /*
       * store_motif
       *
       * Put the motif data in a lookup table so it can be
       * retrieved when required.
       */
      function store_motif(motif_id, name, ltrim, rtrim, pspm_text) {
        var pspm = new Pspm(pspm_text, name, ltrim, rtrim);
        motif_lookup[motif_id] = pspm;
      }

      /*
       * ready_motif
       *
       * Add the information about where a motif needs to be drawn and any
       * modifiers that change the way it is drawn (reverse_complement & 
       * indent) to the datastructures responsible for drawing them.
       *
       * Cache is used as a hint to the drawing code to reduce
       * draw calls by caching the output and using that when future
       * draw calls are made. A duplicate cache entry will be made
       * for each different indent and reverse_complement combination.
       *
       * If group_id is specified then the motifs will be drawn when the
       * group is activated as opposed to when they scroll onto screen
       * which is the default. Hence the group_id is only used for hidden
       * motifs.
       */
      function ready_motif(replace_id, motif_id, cache, indent, reverse_complement, title_text, group_id) {
        var elem = document.getElementById(replace_id);
        var pspm = motif_lookup[motif_id];
        // check we have enough information to make the task possible
        if (elem === null || elem === undefined) {
          if (console) console.error("Can not find element from id %s", replace_id);
          return;
        }
        if (pspm === undefined) {
          if (console) console.error("Missing pspm for %s: unable to replace %s", motif_id, replace_id);
          return;
        }
        // create the task
        var task = new MotifTask(replace_id, motif_id, cache, indent, reverse_complement, title_text);

        // dispatch the task
        if (group_id === undefined) {
          // this is defined in delay_draw.js
          add_draw_task(document.getElementById(replace_id), task);
        } else { 
          add_group_task(group_id, task);
        }
      }

      /*
       * store_graph
       *
       * Put the graph data in a lookup table so it can be 
       * retrieved when required.
       */
      function store_graph(id, graph_text) {
        var graph = new SpamoHistogram(graph_text);
        graph_lookup[id] = graph;
      }

      /*
       * ready_graph
       *
       *
       */
      function ready_graph(replace_id, graph_id, cache, alternate_image, title_text, group_id) {
        var elem = document.getElementById(replace_id);
        var graph = graph_lookup[graph_id];
        // check we have enough information to make the task possible
        if (elem === null || elem === undefined) {
          if (console) console.error("Can not find element from id %s", replace_id);
          return;
        }
        if (graph === undefined) {
          if (console) console.error("Missing graph for %s: unable to replace %s", graph_id, replace_id);
          return;
        }
        // create the task
        var task = new GraphTask(replace_id, graph_id, cache, alternate_image, title_text);

        // dispatch the task
        if (group_id === undefined) {
          // this is defined in delay_draw.js
          add_draw_task(document.getElementById(replace_id), task);
        } else {
          add_group_task(group_id, task);
        }
      }

      /*
       *
       */
      function toggle_group(group_id) {
        //on first run render all the tasks in the group
        if (group_lookup.hasOwnProperty(group_id)) {
          var list = group_lookup[group_id];
          for (var i = 0; i < list.length; i++) {
            list[i].run();
          }
          delete group_lookup[group_id];
        }
        //toggle the visibility
        var elem = document.getElementById(group_id);
        if (elem.style.display == 'none') {
          elem.style.display = 'block';
        } else {
          elem.style.display = 'none';
        }
      }

      /*
       * help_popup
       *
       * Moves around help pop-ups so they appear
       * below an activator.
       */
      function help_popup(activator, popup_id) {
        if (help_popup.popup === undefined) {
          help_popup.popup = null;
        }
        if (help_popup.activator === undefined) {
          help_popup.activator = null;
        }

        if (typeof(activator) == 'undefined') { // no activator so hide
          help_popup.popup.style.display = 'none';
          help_popup.popup = null;
          return;
        }
        var pop = document.getElementById(popup_id);
        if (pop == help_popup.popup) {
          if (activator == help_popup.activator) {
            //hide popup (as we've already shown it for the current help button)
            help_popup.popup.style.display = 'none';
            help_popup.popup = null;
            return; // toggling complete!
          }
        } else if (help_popup.popup != null) {
          //activating different popup so hide current one
          help_popup.popup.style.display = 'none';
        }
        help_popup.popup = pop;
        help_popup.activator = activator;

        //must make the popup visible to measure it or it has zero width
        pop.style.display = 'block';
        var xy = get_elem_xy(activator);
        var padding = 10;
        var edge_padding = 15;
        var scroll_padding = 15;

        var pop_left = (xy[0] + (activator.offsetWidth / 2)  - (pop.offsetWidth / 2));
        var pop_top = (xy[1] + activator.offsetHeight + padding);

        // ensure the box is not past the top or left of the page
        if (pop_left < 0) pop_left = edge_padding;
        if (pop_top < 0) pop_top = edge_padding;
        // ensure the box does not cause horizontal scroll bars
        var page_width = null;
        if (window.innerWidth) {
          page_width = window.innerWidth;
        } else if (document.body) {
          page_width = document.body.clientWidth;
        }
        if (page_width) {
          if (pop_left + pop.offsetWidth > page_width) {
            pop_left = page_width - pop.offsetWidth - edge_padding - scroll_padding; //account for scrollbars
          }
        }

        pop.style.left = pop_left + "px";
        pop.style.top = pop_top + "px";
      }

      /*
       * show_hidden
       *
       * Looks for specially named elements and switches to the shown view
       */
      function show_hidden(prefix) {
        document.getElementById(prefix + '_activator').style.display = 'none';
        document.getElementById(prefix + '_deactivator').style.display = 'block';
        document.getElementById(prefix + '_data').style.display = 'block';
      }
      /*
       * hide_shown
       *
       * Looks for specially named elements and switches to the hidden view
       */
      function hide_shown(prefix) {
        document.getElementById(prefix + '_activator').style.display = 'block';
        document.getElementById(prefix + '_deactivator').style.display = 'none';
        document.getElementById(prefix + '_data').style.display = 'none';
      }
      
  
      /*
       * add_group_task
       *
       *
       */
      function add_group_task(group_id, task) {
        var list;
        if (group_lookup.hasOwnProperty(group_id)) {
          list = group_lookup[group_id]
        } else {
          list = [];
          group_lookup[group_id] = list;
        }
        list.push(task);
      }

      /*
       * MotifTask
       */
      function MotifTask(replace_id, motif_id, cache, indent, reverse_complement, title_text) {
        this.replace_id = replace_id;
        this.motif_id = motif_id;
        this.cache = cache;
        this.indent = indent;
        this.reverse_complement = reverse_complement;
        this.title_text = title_text;
        this.run = MotifTask_run;
      }

      /*
       * MotifTask_run
       */
      function MotifTask_run() {
        var pspm = motif_lookup[this.motif_id];
        if (this.reverse_complement) {
          pspm = pspm.copy().reverse_complement(alphabet);
        }
        var logo = new Logo(alphabet, "MEME");
        logo.add_pspm(pspm, this.indent);
        replace_logo(logo, this.replace_id, 0.5, this.title_text, "block");
      }

      /*
       * GraphTask
       */
      function GraphTask(replace_id, graph_id, cache, alternate_image, title_text) {
        this.replace_id = replace_id;
        this.graph_id = graph_id;
        this.cache = cache;
        this.alternate_image = alternate_image;
        this.title_text = title_text;
        this.run = GraphTask_run;
      }

      /*
       * GraphTask_run
       */
      function GraphTask_run() {
        var element = document.getElementById(this.replace_id);
        var canvas = create_canvas(350, 260, this.replace_id, "", "block");
        if (canvas) {
          var spamo = graph_lookup[this.graph_id];
          draw_spamo_on_canvas(spamo, canvas);
          canvas.title = this.title_text;
          element.parentNode.replaceChild(canvas, element);
        } else {
          var image = document.createElement("img");
          image.width = 350;
          image.height = 260;
          image.src = this.alternate_image;
          image.title = this.title_text;
        }
      }
