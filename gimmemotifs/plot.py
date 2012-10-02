# Copyright (c) 2009-2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Various plotting functions """

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


