import numpy as np
import msprime
import pandas
from scipy.optimize import minimize
from numpy import genfromtxt
from math import log10


####ALSO, JUST TRY GETTING THE ESTIMATES USING ONLY THE ISLAND POPULATION. THEN, SLOWLY ADD MORE
####MAINLAND POPS.
NUMTRIALS = 7000
mu = 5.4*(10**-9)
L = 100000000
empirical = genfromtxt('small_joint.csv', delimiter=' ')
empirical = empirical / np.sum(empirical)

mainids = [0,1,2,3,4,5,6,7]
islandids = [8,9,10,11,12,13,14,15]
nummainhaps = 8
numislandhaps = 8

def kl(E, X):
	print(E)
	return (X * np.log(X / E)).sum()

def getTrees(N_ancestral, N_main,  N_island, T_split, T_mig, migrate):
	demography = msprime.Demography()
	demography.add_population(name="Main", initial_size=N_main)
	demography.add_population(name="Island", initial_size=N_island)
	demography.add_migration_rate_change(time = 0.00000001,  rate = migrate, source="Island", dest="Main")
	demography.add_migration_rate_change(time = T_mig,  rate = 0.0, source="Island", dest="Main")
	demography.add_population(name="ancestral", initial_size=N_ancestral)
	demography.add_population_split(time=T_split, derived=["Main", "Island"], ancestral="ancestral")
	#print(demography)
	return msprime.sim_ancestry({"Main": 4, "Island": 4}, demography=demography, 
						    sequence_length = 1, recombination_rate = 0.0,  
							num_replicates = NUMTRIALS)

def getvalue(N_ancestral, N_main,  N_island, T_split, T_mig, migrate):
	trees = getTrees(N_ancestral, N_main,  N_island, T_split, T_mig, migrate)
	BranchLengths = 0
	Expected = np.full((nummainhaps + 1, numislandhaps + 1), 1/((nummainhaps + 1) * (numislandhaps + 1)))

	for ts in trees:
		BranchLengths = BranchLengths + ts.first().total_branch_length

		SFS = ts.allele_frequency_spectrum([mainids, islandids], mode='branch', polarised=False)
		Expected = Expected + SFS

	Expected[0, 0] = 0.0
	Expected[0, 0] = 1 - mu * (BranchLengths / NUMTRIALS)
	Expected = Expected.flatten()
	Expected[1:]  = (Expected[1:] / np.sum(Expected[1:]) ) * (1 - Expected[0])

	return(kl(Expected[0:(len(Expected)-1)],  empirical.flatten()[0:(len(Expected)-1)] ))

"""
def negloglieklihoodfunction(args):
	N_ancestral = 10**args[0]
	N_main= 10**args[1]
	N_island= 10**args[2]
	T_split= 10**args[3]
	T_mig= 10**args[4]
	migrate= 10**args[5]
	A = getvalue(N_ancestral, N_main,  N_island, T_split, T_mig, migrate)
	print(N_ancestral, N_main,  N_island, T_split, T_mig, migrate, A)

	return(A)

#[436725.41689147 122360.15211858  16474.75155179] 1.3239340884450319e-06   IN HAPLOID UNITS!
#[218362 61180 16474.75155179] 1.3239340884450319e-06   IN diploid UNITS!

res = minimize(negloglieklihoodfunction, [3, 4, 2.5, log10(1000), log10(100), log10(0.01)],
            method='Nelder-Mead', options = {"xatol" : 1e-6}).x

print(res)
"""




def negloglieklihoodfunctionFIX(args):
	N_ancestral = 61180
	N_main= 218362
	N_island= 10**args[0]
	T_split= 16474.75155179
	T_mig= 10**args[1]
	migrate= 10**args[2]
	A = getvalue(N_ancestral, N_main,  N_island, T_split, T_mig, migrate)
	print(N_ancestral, N_main,  N_island, T_split, T_mig, migrate, A)

	return(A)

#[436725.41689147 122360.15211858  16474.75155179] 1.3239340884450319e-06   IN HAPLOID UNITS!
#[218362 61180 16474.75155179] 1.3239340884450319e-06   IN diploid UNITS!

res = minimize(negloglieklihoodfunctionFIX, [3, log10(100),  log10(0.1)],
            method='Nelder-Mead', options = {"xatol" : 1e-6}).x

print(res)
