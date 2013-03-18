# Copyright (c) 2009-2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Various plotting functions """
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from tempfile import NamedTemporaryFile

def axes_off(ax):
    """Get rid of all axis ticks, lines, etc.
    """
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)

def roc_plot(outfile, x, y):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	
	fig = plt.figure()
	try:
		# matplotlib >= 0.99
		rect = fig.patch # a rectangle instance
	except:
		# matplotlib 0.98
		rect = fig.figurePatch # a rectangle instance

	colors = [cm.Paired(256 / 11 * i) for i in range(11)]
	plt.plot(x, y, color=colors[(0 * 2) % 10 + 1])
	plt.axis([0,1,0,1])
	plt.xlabel("1 - Specificity")
	plt.ylabel("Sensitivity")
	plt.savefig(outfile, format="png")

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


        

