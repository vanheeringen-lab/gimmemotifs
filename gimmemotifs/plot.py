# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Various plotting functions """
from __future__ import print_function
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
from matplotlib.colors import to_hex, Normalize
from mpl_toolkits.axes_grid1 import ImageGrid
import seaborn as sns
sns.set_style('white')
from PIL import Image


#from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle

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

    if isinstance(plot_x[0], np.ndarray):
        for i,(x,y) in enumerate(zip(plot_x, plot_y)):
            plt.plot(x, y)
    else:
        plt.plot(plot_x,plot_y)
    
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
    
    plt.clf()
    try:
        # matplotlib >= 0.99
        plt.hist(values, range=xrange, bins=breaks, edgecolor="black")
    except Exception:
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
        motif.to_img(tmp.name, fmt="PNG", height=6)
        grid[0].imshow(plt.imread(tmp.name), interpolation="none")
        tmp = NamedTemporaryFile(dir=mytmpdir(), suffix=".png")
        dbmotif.to_img(tmp.name, fmt="PNG")
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
         
    motifs = np.array(motifs)[filt]
    freq = freq[filt]
    bgfreq = bgfreq[filt]
    enr = enr[filt]
    
    sys.stderr
    for m,f,b,e in zip(motifs,freq,bgfreq,enr):
        sys.stderr.write("{0}\t{1}\t{2}\t{3}\n".format(
            m, 
            "\t".join(str(x) for x in e), 
            "\t".join(str(x) for x in f),
            b[0]))
    
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
    
    sm = plt.cm.ScalarMappable(cmap=c2, norm=Normalize(vmin=vmin, vmax=vmax))
    
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
    # pylint: disable=protected-access
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
    sm = plt.cm.ScalarMappable(cmap=c1, norm=Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cax = plt.subplot(gs[0,plot_order[1]])
    cb = fig.colorbar(sm, cax=cax, ticks = [vmin,0, vmax], orientation='horizontal')
    cb.ax.set_xticklabels([vmin, 0, vmax])
   
   
    #### Motif logos
   
    for i,motif in enumerate(motifs[ind][::-1]):
        ax = plt.subplot(gs[i + nbar, plot_order[0]]) 
        axes_off(ax)
        tmp = NamedTemporaryFile(dir=mytmpdir(), suffix=".png")
        pwms[motif].to_img(tmp.name, fmt="PNG", height=6)
        ax.imshow(plt.imread(tmp.name), interpolation="none")
    
    #plt.show()
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close(fig)

def _tree_layout(node):
    if node.is_leaf():
        nameFace = AttrFace("name", fsize=24, ftype="Nimbus Sans L")
        faces.add_face_to_node(nameFace, node, 10, position="branch-right")

def _get_motif_tree(tree, data, circle=True, vmin=None, vmax=None):
    print(circle, vmin, vmax)
    t = Tree(tree)
    # Determine cutoff for color scale
    if not(vmin and vmax):
        for i in range(90, 101):
            minmax = np.percentile(data.values, i)
            if minmax > 0:
                break
    if not vmin:
        vmin = -minmax
    if not vmax:
        vmax = minmax
    print(vmin, vmax)
    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap="coolwarm")
    
    m = 25 / data.values.max()
    
    for node in t.traverse("levelorder"):
        val = data[[l.name for l in node.get_leaves()]].values.mean()
        style = NodeStyle()
        style["size"] = 0
        
        style["hz_line_color"] = to_hex(mapper.to_rgba(val))
        style["vt_line_color"] = to_hex(mapper.to_rgba(val))
        
        v = max(np.abs(m * val), 5)
        style["vt_line_width"] = v
        style["hz_line_width"] = v

        node.set_style(style)
    
    ts = TreeStyle()

    ts.layout_fn = _tree_layout
    ts.show_leaf_name= False
    ts.show_scale = False
    ts.branch_vertical_margin = 10

    if circle:
        ts.mode = "c"
        ts.arc_start = 180 # 0 degrees = 3 o'clock
        ts.arc_span = 180
    
    return t, ts

def motif_tree_plot(outfile, tree, data, circle=True, vmin=None, vmax=None, dpi=300):
    """
    Plot a "phylogenetic" tree 
    """
    # Define the tree
    t, ts = _get_motif_tree(tree, data, circle, vmin, vmax)
    
    # Save image
    t.render(outfile, tree_style=ts, w=100, dpi=dpi, units="mm");
    
    # Remove the bottom (empty) half of the figure
    if circle:
        img = Image.open(outfile)
        size = img.size[0]
        spacer = 50
        img.crop((0,0,size,size/2 + spacer)).save(outfile)
