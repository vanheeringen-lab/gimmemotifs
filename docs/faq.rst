.. _`FAQ`:

FAQ
===

I have upgraded GimmeMotifs and now it doesn't find my genome
-------------------------------------------------------------

The genome index in GimmeMotifs has changed, see upgradegenome_.


I cannot run gimme index anymore
--------------------------------

The genome index in GimmeMotifs has changed, see upgradegenome_.


I get 'RuntimeError: Invalid DISPLAY variable'
----------------------------------------------

The default matplotlib configuration expects a display. Probably you are running GimmeMotifs on a server without an X server. There are several ways to solve it.

Option 1
~~~~~~~~

Change the matplotlib configuration. Find the path of the matplotlib installation of your current environment (make sure to activate the environment you use for GimmeMotifs first).

::

    $ python -c "import matplotlib; print(matplotlib.__file__)"
    /home/simon/anaconda2/envs/gimme3/lib/python3.5/site-packages/matplotlib/__init__.py

So, matplotlib is in ``/home/simon/anaconda2/envs/gimme3/lib/python3.5/site-packages/matplotlib/``.
Now you can edit ``<MATPLOT_DIR>/mpl-data/matplotlibrc``. In my example case this would be:

``/home/simon/anaconda2/envs/gimme3/lib/python3.5/site-packages/matplotlib/mpl-data/matplotlibrc``

Change the line

::

    backend     : Qt5Agg

to

::

    backend     : Agg


You can also put a matplotlibrc file in ``$HOME/.config/matplotlib``.

Option 2
~~~~~~~~

Run GimmeMotifs via ``xvfb-run``. If this program is installed, you can simply run GimmeMotifs in a virtual X server environment.

For example:

:: 

    $ xvfb-run gimme motifs [args]


I get a KeyError when running gimme maelstrom
---------------------------------------------

You get an error like this:

::

    File "pandas/_libs/index.pyx", line 132, in pandas._libs.index.IndexEngine.get_loc (pandas/_libs/index.c:5280)
    File "pandas/_libs/index.pyx", line 154, in pandas._libs.index.IndexEngine.get_loc (pandas/_libs/index.c:5126)
    File "pandas/_libs/hashtable_class_helper.pxi", line 1210, in pandas._libs.hashtable.PyObjectHashTable.get_item (pandas/_libs/hashtable.c:20523)
    File "pandas/_libs/hashtable_class_helper.pxi", line 1218, in pandas._libs.hashtable.PyObjectHashTable.get_item (pandas/_libs/hashtable.c:20477)
    KeyError: '5'

This a bug in ``gimme maelstrom``. The column headers can't be numbers. Change this to a word, for instance ``cluster5`` or ``col5``.
