# Cian Steenkamp 
import pygmo as pg
import numpy as np
import sys
import os
np.set_printoptions(threshold=sys.maxsize)

# Testing
# points = [	[0, 1], # A
# 			[0, 3], # B
# 			[1, 0], # C
# 			[0, 2]	# D
# 		 ]
# 		 	# < = dominates
# 		 	# A < B
# 		 	# A < D
# 		 	# D < B
# 		 	# NDF = A and C since no solution(s) dominate A and C

# # points = [[0, 1], [0, 3], [0.8, 0], [0, 0.9], [0.7, 0]]
# # points = [[0,1],[-1,3],[2.3,-0.2],[1.1,-0.12],[1.1, 2.12],[-1.1,-1.1]]
# # points = [[0, 1], [1, 0]]

# nadir = pg.nadir(points) # The nadir is that point that has the maximum value of the objective function in the points of the non-dominated front.
# ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(points)
# print nadir
# print ndf[0] # list of indices of nd solutions

fronts = os.listdir('raw/')
for front in fronts:
	print front
	original_pf = np.loadtxt('raw/' + front)
	ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(original_pf)

	if len(ndf[0]) < len(original_pf): # remove dominated and store non-dom new file
		print "Removing dominated and saving to file"
		new_pf = []
		d = [] # dominated sols indices
		for i in range(len(original_pf)):
			if i in ndf[0]:
				new_pf.append(original_pf[i])
			else:
				d.append(i)

		print 'Number of dominated sols => ', len(d) # number of dominated sols
		new_pf = np.array(new_pf).astype(np.float)
		np.savetxt('nonDominated/' + front, new_pf)
	else: # store non-dominated new file (no dominated)
		# print "None dominated; saving to file"
		new_pf = np.array(original_pf).astype(np.float)
		np.savetxt('nonDominated/' + front, new_pf)
