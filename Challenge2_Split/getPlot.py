import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from numpy import genfromtxt
from math import log10


mu = 5.7*(10**-9)
L = 100000000 
numberofunits = 500
repitations = 50
empirical = genfromtxt('challenge2_joint.csv', delimiter=' ')
empirical = empirical / np.sum(empirical)

eastids = []
for i in range(0,44):
	eastids.append(i)
westids = []
for i in range(44,44+36):
	westids.append(i)

def kl(E, X):
	test1 = 0
	test2 =  0

	largesum = 0
	for i in range(45):
		for j in range(37):
			if E[i,j] > 0 and X[i,j] > 0:
				largesum = largesum + X[i,j] * np.log(X[i,j]  / E[i,j])
				test1 = test1 +  X[i,j] 
				test2 = test2 +  E[i,j] 
	print("test", test1, test2)
	return np.sum(largesum)

def getTrees(N_east, N_west,  N_ancestral, T_split):
	demography = msprime.Demography()
	demography.add_population(name="east", initial_size=N_east)
	demography.add_population(name="west", initial_size=N_west)
	demography.add_population(name="ancestral", initial_size=N_ancestral)
	demography.add_population_split(time=T_split, derived=["east", "west"], ancestral="ancestral")
	return msprime.sim_ancestry({"east": 22, "west": 18}, demography=demography, 
						    sequence_length = L / numberofunits, recombination_rate = 3.386 * 10**(-9),  
							num_replicates = numberofunits)

def getvalue(N_east, N_west,  N_ancestral, T_split):
	SFSFULL = np.zeros((45,37))
	for _ in range(repitations):
		SFS = np.zeros((45,37))
		tre = getTrees(N_east, N_west,  N_ancestral, T_split)
		totalmut = 0

		for ts in tre:
			mts = msprime.sim_mutations(ts, rate=mu  , model = "binary")
			totalmut = totalmut + np.sum(mts.allele_frequency_spectrum([eastids, westids], mode='site', polarised=True, span_normalise = False))
			SFS = SFS + mts.allele_frequency_spectrum([eastids, westids], mode='branch', polarised=True, span_normalise = False)
		SFS = np.divide(SFS, np.divide(np.sum(SFS)  , totalmut) )

		SFS[0,0]  = L - np.sum(SFS)
		SFS = np.divide(SFS, L)
		SFSFULL = np.add(SFSFULL , SFS)
	return np.divide(SFSFULL , repitations)

N_east = 125992
N_west=20109
N_ancestral= 100000
T_split= 13302.758875
A = getvalue(N_east, N_west,  N_ancestral, T_split)


import matplotlib.pyplot as plt 

A[0,0] = 0
empirical[0,0] = 0
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow( np.log(empirical) )
ax1.set_xlabel('West')
ax1.set_ylabel('East')
ax1.set_title("Empirical")

ax2.imshow( np.log(A) )
ax2.set_xlabel('West')
ax2.set_ylabel('East')
ax2.set_title("Expectation")
fig.tight_layout()
plt.savefig('Model', dpi=1000)