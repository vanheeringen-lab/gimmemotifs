.. _`FAQ`:

FAQ
===

SQLite error when running on a cluster
--------------------------------------

The current implementation of the cache that GimmeMotifs uses does not play nice with concurrent access, for instance on
 a cluster. The result is that the cache will get corrupted and that GimmeMotifs will fail. Until this is fixed, there i
s a workaround. In your job submission script, use something like the following:

::

    NEW_CACHE=$TMPDIR/cache
    mkdir -p $NEW_CACHE
    if [ -z $XDG_CACHE_HOME ]; then
        XDG_CACHE_HOME=$HOME/.cache
    fi
    cp -r $XDG_CACHE_HOME/gimmemotifs $NEW_CACHE/
    export XDG_CACHE_HOME=$NEW_CACHE
    echo 'Using $XDG_CACHE_HOME for cache'

This will use a local directory to store the cache.


ImportError: dlopen: cannot load any more object with static TLS
----------------------------------------------------------------

This is `a bug <https://github.com/pytorch/pytorch/issues/2575#issuecomment-523667875>`_
that will occur on Ubuntu 16.04 (Trusty). 
The recommended solution is to upgrade Ubuntu as Trusty is no longer supported. 
However, if you are unable to do so, there is a workaround.
Run the following command:

::

    $ export LD_PRELOAD=$CONDA_PREFIX/lib/libgomp.so

Now you should be able to run `gimme` as normal.

'i' format requires -2147483648 <= number <= -2147483646
-------------------------------------------------------

If you get the following error with `gimme maelstrom`:

```
'i' format requires -2147483648 <= number <= -2147483646
```

this means that you should decrease the amount of input sequences. 
It is a bug that should be solved with a newer version of Python. However, it might be bestg anyway to use a limited
set of sequences (< ~100k) as input.

Sorry, motif prediction tool [X] is not supported
-------------------------------------------------

If `gimme motifs` does not recognize a motif tools that should be installed, you can update the configuration file. This file is likely located at `~/.config/gimmemotifs/gimmemotifs.cfg`.

Edit this file and update the following lines under the `[params]` section:

::

    available_tools = MDmodule,MEME,MEMEW,Weeder,GADEM,MotifSampler,Trawler,Improbizer,BioProspector,Posmo,ChIPMunk,AMD,HMS,Homer
    tools = MDmodule,MEME,Weeder,MotifSampler,trawler,Improbizer,BioProspector,Posmo,ChIPMunk,AMD,Homer

Add the tool that you want to the `available_tools` parameter. Keep in mind the exact upper-/lower-case combination that GimmeMotifs uses. By updating the `tools` parameter you can set the tools that `gimme motifs` uses by default. This can always be changed at the command-line. 

In addition, you might also have to update the binary location. Either update this section if it exists or add it. For instance, to set this information for MEME:

::

    [MEME]
    bin = /home/simon/anaconda2/envs/gimme/gimmemotifs/tools/meme.bin
    dir = /home/simon/anaconda2/envs/gimme/gimmemotifs/tools

The `dir` variable ususally doesn't matter and you can set it the same directory as where the binary is located.

I get motifs that have differential scores in gimme maelstrom, however, the number is not different across clusters
-------------------------------------------------------------------------------------------------------------------

The different methods use different ways to rank the motifs. The hypergeometric test is the only one that uses motif counts. All the other methods use the PWM logodds score of the best match. While the counts may not be different across clusters, the scores most likely are.

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
