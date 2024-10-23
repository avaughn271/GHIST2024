import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from math import log10


mu = 1.29 * (10**-8)
L = 250000000
numberofunits = 500
repitations =  1 ##25 25

empirical = (pandas.read_csv('Marginal/sfs1.csv', sep = " ", header  = None))
empirical = np.array(empirical.iloc[:,1])
empirical = np.insert(empirical, 0, L - np.sum(empirical))
empirical  = empirical/np.sum(empirical)
print(empirical)


islandids = []
for i in range(0,40):
	islandids.append(i)

def kl(E, X):
	largesum = 0
	test1 = 0
	for j in range(40):
		if E[j] > 0 and X[j] > 0:
			largesum = largesum + X[j] * np.log(X[j]  / E[j])
			test1 = test1 +  X[j] 
		if E[j] <= 0 and X[j] > 0:
			print(j, X[j])
			largesum = largesum + X[j] * (np.log(X[j]) + 30.0  )
	print("test", test1)
	return np.sum(largesum)

def getTrees(N_neander, N_human,  N_ancestral, T_split, T_mig, prop):
	demography = msprime.Demography()

	demography.add_population(name="ancestral", initial_size=N_ancestral)
	demography.add_population(name="human", initial_size=N_human)
	demography.add_population(name="ancestral_human", initial_size=N_human)
	demography.add_population(name="neander", initial_size=N_neander)

	demography.add_admixture(time=T_mig, derived="human", ancestral=["ancestral_human", "neander"], proportions=[1-prop, prop])

	demography.add_population_split(time=T_split, derived=["ancestral_human", "neander"], ancestral="ancestral")
	return msprime.sim_ancestry({"human": 20}, demography=demography, 
						    sequence_length = L / numberofunits, recombination_rate = 1.38 * 10**(-8),  
							num_replicates = numberofunits)

def getvalue(N_main, N_island,  N_ancestral, T_split, T_mig, migrate):
	SFSFULL = np.zeros(41)
	for _ in range(repitations):
		SFS = np.zeros(41)
		tre = getTrees(N_main, N_island,  N_ancestral, T_split, T_mig, migrate)
		totalmut = 0

		for ts in tre:
			mts = msprime.sim_mutations(ts, rate=mu  , model = "binary")
			totalmut = totalmut + np.sum(mts.allele_frequency_spectrum([islandids], mode='site', polarised=True, span_normalise = False))
			SFS = SFS + mts.allele_frequency_spectrum([islandids], mode='branch', polarised=True, span_normalise = False)
		SFS = np.divide(SFS, np.divide(np.sum(SFS)  , totalmut) )

		SFS[0]  = L - np.sum(SFS)
		SFS = np.divide(SFS, L)
		SFSFULL = np.add(SFSFULL , SFS)
	return kl(np.divide(SFSFULL , repitations), empirical)

def negloglieklihoodfunction(args):
	N_neander = 10**args[0]
	N_human= 10**args[1]
	N_ancestral= 10**args[2]
	T_split= 10**args[3]
	T_mig= 10**args[4]
	prop=  0.011

	if T_split < T_mig: return(100000000000.00)
	if prop < 0 or prop > 1: return(100000000000.00)

	A = getvalue(N_neander, N_human,  N_ancestral, T_split, T_mig, prop)
	print(N_neander, N_human,  N_ancestral, T_split, T_mig, prop, A)
	print("B")

	return((A*100000)**2)

print("A")
res = (minimize(negloglieklihoodfunction, [
										   log10(20000), #neadnder
											 log10(4500), # human 
											  log10( 30000),  # ancestral
											  	log10(10000),  #split  time
											  log10(100) 
											 ], method='Nelder-Mead'))
print(10**res.x)

