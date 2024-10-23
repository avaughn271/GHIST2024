import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from numpy import genfromtxt
from math import log10

for name in [1, 2]:
	for j in ["all"]:
		empirical = genfromtxt(str(name)+"_"  + str(j) + '.csv', delimiter=' ')
		empirical = empirical / np.sum(empirical)

		import matplotlib.pyplot as plt

		empirical[0,0] = 0
		plt.imshow( np.log(empirical) )
		plt.savefig('Model' + str(name) + "_"  + str(j) , dpi=1000)