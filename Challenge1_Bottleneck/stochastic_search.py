import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from numpy import genfromtxt
from math import log10


mu = 1.26*(10**-8)
L = 100000000
numberofunits = 500
repitations = 25

empirical = (pandas.read_csv("sfs.csv", sep = " ", header  = None))
empirical = np.array(empirical.iloc[:,1])
empirical = np.insert(empirical, 0, L - np.sum(empirical))
empirical = empirical / np.sum(empirical)
print(empirical)

def kl(E, X):
    return (X * np.log(X / E)).sum()

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
	return kl(expectation[0:40], empirical)

def negloglieklihoodfunction(args):
	recentne = 10**args[0]
	oldne = 10**args[1]
	gensplit = 10**args[2]
	A = getvalue(recentne, gensplit, oldne)
	print(recentne, oldne, gensplit, A)

	return(A)

res = minimize(negloglieklihoodfunction, [3, 4, 2],
            method='Nelder-Mead', options = {"xatol" : 1e-6}).x
print(res)
