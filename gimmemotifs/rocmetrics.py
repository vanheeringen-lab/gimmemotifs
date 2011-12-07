# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Module to calculate ROC and MNCP scores """


# External imports
from scipy.stats import stats
from numpy import *

def MNCP(fg_vals, bg_vals):
	from scipy.stats import stats
	from numpy import mean
	#from pylab import *
	fg_len = len(fg_vals)
	total_len = len(fg_vals) + len(bg_vals)

	fg_rank = stats.rankdata(fg_vals)
	total_rank = stats.rankdata(fg_vals + bg_vals)

	slopes = []
	for i in range(len(fg_vals)):
		slope = ((fg_len - fg_rank[i] + 1) / fg_len ) / ((total_len - total_rank[i] + 1)/ total_len)
		slopes.append(slope)
	return mean(slopes)

def ROC_AUC(fg_vals, bg_vals):
	#if len(fg_vals) != len(bg_vals):
	#	return None
	
	if len(fg_vals) == 0 or len(bg_vals) == 0:
		return None
	
	fg_len = len(fg_vals)
	total_len = len(fg_vals) + len(bg_vals)

	fg_rank = stats.rankdata(fg_vals) 
	total_rank = stats.rankdata(fg_vals + bg_vals) 
	
	return (sum(total_rank[:fg_len]) - sum(fg_rank))/ (fg_len * (total_len - fg_len))

def ROC_AUC_xlim(x_bla, y_bla, xlim=None):
	x = x_bla[:]
	y = y_bla[:]

	x.sort()
	y.sort()

	u = {}
	for i in x + y:
		u[i] = 1

	vals = u.keys()
	vals.sort()
	
	len_x = float(len(x))
	len_y = float(len(y))
	
	new_x = []
	new_y = []
	
	x_p = 0
	y_p = 0
	for val in vals[::-1]:
		while len(x) > 0 and x[-1] >= val:
			x.pop()
			x_p += 1
		while len(y) > 0 and y[-1] >= val:
			y.pop()
			y_p += 1
		new_y.append((len_x - x_p) / len_x)
		new_x.append((len_y - y_p) / len_y)
	
	#print new_x
	#print new_y
	new_x = 1 - array(new_x)
	new_y = 1 - array(new_y)
	#plot(new_x, new_y)
	#show()

	x = new_x
	y = new_y

	if len(x) != len(y):
		raise "Unequal!"

	if not xlim:
		xlim = 1.0

	auc = 0.0
	bla = zip(stats.rankdata(x), range(len(x)))

	def sortfunc(x,y):
		res = x[0] - y[0]
		if res < 0:
			return -1
		elif res > 0:
			return 1
		elif res == 0:
			return y[1] - x[1]
	
	bla.sort(sortfunc)
	
	prev_x = x[bla[0][1]]
	prev_y = y[bla[0][1]]
	index = 1

	while index < len(bla) and x[bla[index][1]] <= xlim:

		(rank, i) = bla[index]
		
		auc += y[i] * (x[i] - prev_x) - ((x[i] - prev_x) * (y[i] - prev_y) / 2.0)
		prev_x = x[i]
		prev_y = y[i]
		index += 1
	
	if index < len(bla):
		(rank, i) = bla[index]
		auc += prev_y * (xlim - prev_x) + ((y[i] - prev_y)/(x[i] - prev_x) * (xlim -prev_x) * (xlim - prev_x)/2)

	return auc

def ROC_values(x_bla, y_bla):
	if len(x_bla) == 0 or len(y_bla) == 0:
		return [],[]

	x = x_bla[:]
	y = y_bla[:]

	x.sort()
	y.sort()

	u = {}
	for i in x + y:
		u[i] = 1

	vals = u.keys()
	vals.sort()
	
	len_x = float(len(x))
	len_y = float(len(y))
	
	new_x = []
	new_y = []
	
	x_p = 0
	y_p = 0
	for val in vals[::-1]:
		while len(x) > 0 and x[-1] >= val:
			x.pop()
			x_p += 1
		while len(y) > 0 and y[-1] >= val:
			y.pop()
			y_p += 1
		new_y.append((len_x - x_p) / len_x)
		new_x.append((len_y - y_p) / len_y)
	
	#print new_x
	#print new_y
	new_x = 1 - array(new_x)
	new_y = 1 - array(new_y)
	
	return (new_x, new_y)

def max_fmeasure(x,y):
	x = array(x[:])
	y = array(y[:])
	p = y / (y + x)
	f = (2 * p * y) / (p + y)
	if len(f) > 0:
		return nanmax(f), nanmax(y[f == nanmax(f)])
	else:
		return None,None


#print ROC_AUC([0,0,0,1,2,3,3,3,3,3],[2,2,2,2,2])
#print ROC_AUC_xlim([0,0,0,1,2,3,3,3,3,3],[2,2,2,2,2], 1.0)


