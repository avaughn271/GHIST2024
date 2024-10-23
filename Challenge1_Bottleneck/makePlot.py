import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from numpy import genfromtxt
from math import log10


mu = 1.26*(10**-8)
L = 100000000
numberofunits = 500
repitations = 50

empirical = (pandas.read_csv("sfs.csv", sep = " ", header  = None))
empirical = np.array(empirical.iloc[:,1])
empirical = np.insert(empirical, 0, L - np.sum(empirical))
empirical = empirical / np.sum(empirical)

def getTrees(N_recent, SPLIT,  N_ancestral):
	demography = msprime.Demography()
	demography.add_population(name="wisent", initial_size = N_recent)
	demography.add_population_parameters_change(time = SPLIT, initial_size = N_ancestral, population=0)
	return msprime.sim_ancestry({"wisent": 20}, demography=demography, 
						    sequence_length = L / numberofunits, recombination_rate = 1.007 * 10**(-8),  
							num_replicates = numberofunits)

def getvalue(N_recent, SPLIT,  N_ancestral):
	SFSFULL = np.zeros(41)
	for _ in range(repitations):
		SFS = np.zeros(41)
		tre = getTrees(N_recent, SPLIT,  N_ancestral)
		totalmut = 0

		for ts in tre:
			mts = msprime.sim_mutations(ts, rate=mu  , model = "binary")
			totalmut = totalmut + np.sum(mts.allele_frequency_spectrum(mode='site', polarised=True, span_normalise = False))
			SFS = SFS + mts.allele_frequency_spectrum(mode='branch', polarised=True, span_normalise = False)
		SFS = np.divide(SFS, np.divide(np.sum(SFS) , totalmut))

		SFS[0]  = L - np.sum(SFS)
		SFS = np.divide(SFS, L)

		SFSFULL = np.add(SFSFULL , SFS)
	expectation = np.divide(SFSFULL , repitations)
	return expectation[0:40]

###haploid Ne is: [ 2412.06052557 28159.91681494   302.5125612 ]

expectation = getvalue(1130, 285, 14125)

import matplotlib.pyplot as plt 

# plot lines
plt.plot(np.arange(1,40, step = 1), np.log10(empirical[1:]), 'o')
plt.plot(np.arange(1,40, step = 1), np.log10(expectation[1:]), label = "line 1")

plt.savefig('Model', dpi=1000)