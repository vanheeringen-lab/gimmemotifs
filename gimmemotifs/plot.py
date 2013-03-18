# Copyright (c) 2009-2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Various plotting functions """
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from tempfile import NamedTemporaryFile
import numpy
import os

VALID_EXTENSIONS = [".png", ".pdf", ".svg", ".ps"]

def axes_off(ax):
    """Get rid of all axis ticks, lines, etc.
    """
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)

def roc_plot(outfile, plot_x, plot_y, ids=[]):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    fig = plt.figure()
    fig.add_subplot(111, aspect="equal")

    colors = [cm.Paired(256 / 11 * i) for i in range(11)]
    
    if type(plot_x[0]) == type(numpy.array([])):
        for i,(x,y) in enumerate(zip(plot_x, plot_y)):
            plt.plot(x, y, color=colors[(i * 2) % 10 + 1])
    else:
        plt.plot(plot_x,plot_y, color=colors[(0 * 2) % 10 + 1])
    
    plt.axis([0,1,0,1])
    plt.xlabel("1 - Specificity")
    plt.ylabel("Sensitivity")

    if len(ids) > 0:
         plt.legend(ids, loc=(1.03,0.2))

    if not os.path.splitext(outfile)[-1] in VALID_EXTENSIONS:
        outfile += ".png"
      
    plt.savefig(outfile, dpi=300, bbox_inches='tight')

def plot_histogram(values, outfile, xrange=None, breaks=10, title=None, xlabel=None, color=10):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    
    colors = [cm.Paired(256 / 11 * i) for i in range(11)]

    plt.clf()
    try:
        # matplotlib >= 0.99
        plt.hist(values, range=xrange, bins=breaks, color=colors[color], edgecolor="black")
    except:
        plt.hist(values, range=xrange, bins=breaks)
    plt.xlim(xrange)

    if title:
        plt.title(title)

    plt.ylabel("Frequency")
    if xlabel:
        plt.xlabel(xlabel)
    if not outfile.endswith(".svg"):
        outfile += ".svg"
    plt.savefig(outfile, format="svg")

def match_plot(plotdata, outfile):
    """Plot list of motifs with database match and p-value
    "param plotdata: list of (motif, dbmotif, pval)
    """
    fig_h = 2 
    fig_w = 7

    nrows = len(plotdata)
    ncols = 2
    fig = plt.figure(figsize=(fig_w, nrows * fig_h))
    
    for i, (motif, dbmotif, pval) in enumerate(plotdata):
        text = "Motif: %s\nBest match: %s\np-value: %0.2e" % (motif.id, dbmotif.id, pval)
        

        grid = ImageGrid(fig, (nrows, ncols, i * 2 + 1), 
                         nrows_ncols = (2,1),
                         axes_pad=0, 
                         )

        for j in range(2):  
            axes_off(grid[j])

        tmp = NamedTemporaryFile(suffix=".png")
        motif.to_img(tmp.name, format="PNG")
        grid[0].imshow(plt.imread(tmp.name), interpolation="none")
        tmp = NamedTemporaryFile(suffix=".png")
        dbmotif.to_img(tmp.name, format="PNG")
        grid[1].imshow(plt.imread(tmp.name), interpolation="none")

        ax = plt.subplot(nrows, ncols, i * 2 + 2)
        axes_off(ax)

        ax.text(0, 0.5, text,
        horizontalalignment='left',
        verticalalignment='center') 
    
    plt.savefig(outfile, dpi=300, bbox_inches='tight')


        

