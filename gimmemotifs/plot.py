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
	plt.savefig(roc_img_file % (id,name), format="png")

