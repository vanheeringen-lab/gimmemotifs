# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Various plotting functions """
import os
import sys
from tempfile import NamedTemporaryFile
import numpy as np

# Clustering
from scipy.cluster import hierarchy as hier
from gimmemotifs import mytmpdir

# Matplotlib imports
import matplotlib as mpl
mpl.use("Agg", warn=False)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import ImageGrid

VALID_EXTENSIONS = [".png", ".pdf", ".svg", ".ps"]

def axes_off(ax):
    """Get rid of all axis ticks, lines, etc.
    """
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)

def roc_plot(outfile, plot_x, plot_y, ids=None):
    if ids is None:
        ids = []

   
    fig = plt.figure()
    fig.add_subplot(111, aspect="equal")

    colors = [cm.Paired(256 / 11 * i) for i in range(11)]
    
    if type(plot_x[0]) == type(np.array([])):
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
    plt.close(fig)

def plot_histogram(values, outfile, xrange=None, breaks=10, title=None, xlabel=None, color=10):
    
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

        tmp = NamedTemporaryFile(dir=mytmpdir(), suffix=".png")
        motif.to_img(tmp.name, format="PNG", height=6)
        grid[0].imshow(plt.imread(tmp.name), interpolation="none")
        tmp = NamedTemporaryFile(dir=mytmpdir(), suffix=".png")
        dbmotif.to_img(tmp.name, format="PNG")
        grid[1].imshow(plt.imread(tmp.name), interpolation="none")

        ax = plt.subplot(nrows, ncols, i * 2 + 2)
        axes_off(ax)

        ax.text(0, 0.5, text,
        horizontalalignment='left',
        verticalalignment='center') 
    
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close(fig)

def diff_plot(motifs, pwms, names, freq, counts, bgfreq, bgcounts, outfile, mindiff=0, minenr=3, minfreq=0.01):
    w_ratio = np.array([14, len(names), len(names) + 1])
    plot_order = [0,1,2]
    
    nbar = 5
    
    freq = np.array(freq)
    counts = np.array(counts)
    bgfreq = np.array([[x] for x in bgfreq])
    
    enr = np.log2(np.divide(freq, bgfreq))
   
    filt = np.ones(len(enr), dtype="bool")
    filters = [
                np.sum(enr > minenr, 1) > 0, 
                np.sum(freq > minfreq, 1) > 0,
                (np.max(enr, 1) - np.min(enr, 1)) > mindiff,
                np.sum(counts > 2, 1) > 0 
              ]
    for f in filters:
        filt = np.logical_and(filt, f)
         
        print "Filter: ", sum(filt)
    

    motifs = np.array(motifs)[filt]
    freq = freq[filt]
    bgfreq = bgfreq[filt]
    enr = enr[filt]
    
    for m,f,b,e in zip(motifs,freq,bgfreq,enr):
        sys.stderr.write("{0}\t{1}\t{2}\t{3}\n".format(m,f,b,e))
    
    
    if len(freq) == 0:
        sys.stderr.write("No enriched and/or differential motifs found.\n")
        return
    elif len(freq) >= 3:
        z = hier.linkage(freq, method="complete", metric="correlation")
        ind = hier.leaves_list(z)
    else:
        ind = np.arange(len(freq))
   
    fig = plt.figure(figsize=(
                (5 + 0.75 * len(names)) * 3,
                (0.3 * len(motifs) + 1.5) * 3
            ))
    
    gs = GridSpec(len(motifs) + 3 + nbar, 3,
                  height_ratios=[1] * nbar + [3] * (len(motifs) + 3),
                  width_ratios=w_ratio[plot_order],
                  )
    
    # Colormaps
    c1 = mpl.cm.RdBu
    c2 = mpl.cm.Blues ##create_colormap("white", "blue")
    
    ### Frequency plot ###
    
    # Create axis
    ax = plt.subplot(gs[nbar:-3, plot_order[2]])
    
    # Plot frequencies
    vmin = 0
    vmax = 0.3
    
    pfreq = np.hstack((freq, bgfreq))
    ax.pcolormesh(pfreq[ind], cmap=c2, vmin=vmin, vmax=vmax)
    
    sm = plt.cm.ScalarMappable(cmap=c2, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    
    # Show percentages
    for y,row in enumerate(pfreq[ind]):
        for x,val in enumerate(row):
            v = vmax
            if val >= (vmin + ((vmax - vmin) / 2)):
                v = vmin        
            plt.text(x + 0.5, y + 0.5, "{:.1%}".format(val), ha='center', va='center', color=sm.to_rgba(v))
    
    # Hide most labels
    plt.setp(ax.get_xticklines(),visible=False)
    plt.setp(ax.get_yticklines(),visible=False)
    plt.setp(ax.get_yticklabels(),visible=False)
    
    # Set the X labels
    ticks = np.arange(len(names)+ 1) + 0.5
    plt.xticks(ticks, names + ["background"], rotation=30, ha="right")

    ax.set_ylim(0, len(motifs))

    # Title
    plt.title('Frequency')
    
    # Colorbar
    sm._A = []
    cax = plt.subplot(gs[0,plot_order[2]])
    cb = fig.colorbar(sm, cax=cax, ticks = [0, 0.3], orientation='horizontal')
    cb.ax.set_xticklabels(["0%","30%"])
   

    #### Enrichment plot
    ax = plt.subplot(gs[nbar:-3, plot_order[1]])
    vmin = -10
    vmax = 10
    ax.pcolormesh(enr[ind], cmap=c1, vmin=vmin, vmax=vmax)
    for y,row in enumerate(enr[ind]):
        for x,val in enumerate(row):
            col = "black"
            if val >= (vmin + ((vmax - vmin) / 8.0 * 7)):
                col = "white"
            elif val <= (vmin + ((vmax - vmin) / 8.0)):
                col = "white"
            plt.text(x + 0.5, y + 0.5, "{:.1f}".format(val), ha='center', va='center', color=col)
    
    ticks = np.arange(len(names)) + 0.5
    plt.xticks(ticks, names, rotation=30, ha="right")
    #plt.setp(plt.xticks()[1], rotation=30)
    #for label in labels: 
    #    label.set_rotation(30)
    ticks = np.arange(len(motifs)) + 0.5
    plt.yticks(ticks, motifs[ind])
    plt.setp(ax.get_xticklines(),visible=False)
    plt.setp(ax.get_yticklines(),visible=False)
    
    ax.set_ylim(0, len(motifs))
    
    # Title
    plt.title('Enrichment (log2)')
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=c1, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cax = plt.subplot(gs[0,plot_order[1]])
    cb = fig.colorbar(sm, cax=cax, ticks = [vmin,0, vmax], orientation='horizontal')
    cb.ax.set_xticklabels([vmin, 0, vmax])
   
   
    #### Motif logos
   
    for i,motif in enumerate(motifs[ind][::-1]):
        ax = plt.subplot(gs[i + nbar, plot_order[0]]) 
        axes_off(ax)
        tmp = NamedTemporaryFile(dir=mytmpdir(), suffix=".png")
        pwms[motif].to_img(tmp.name, format="PNG", height=6)
        ax.imshow(plt.imread(tmp.name), interpolation="none")
    
    #plt.show()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close(fig)
