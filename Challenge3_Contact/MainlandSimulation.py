import numpy as np
import msprime
import pandas
from scipy.optimize import minimize

NUMTRIALS = 10000
mu = 5.4*(10**-9)
L = 100000000
numhaps = 44
empirical = (pandas.read_csv("mainland.csv", sep = " ", header  = None))
empirical = np.array(empirical.iloc[:,1])
empirical = np.insert(empirical, 0, L - np.sum(empirical))
empirical = empirical / np.sum(empirical)

#Don't bother adding the mutations. Just find the relative proportion of the set of branches that
#subtend a certain number of leaves. this is good for all nonzero vals.
#For vals, just get average total tree length and multiply by the mutation rate.
##check for sfs against total one.
##saves time in doing the mutations maybe??
##check agaisnt fastsimcoal.
##hopefulyl fastsimcoal does this already. If not, then can improve on it.

def kl(E, X):
    return (X * np.log(X / E)).sum()


def getTrees(N_recent, N_old,  SPLIT):
	demography = msprime.Demography()
	demography.add_population(name="A", initial_size = N_recent)
	demography.add_population_parameters_change(time = SPLIT, initial_size = N_old,  
											 population=0)
	return msprime.sim_ancestry({"A": numhaps/2}, demography=demography, 
						    sequence_length = 1, recombination_rate = 0.0,  
							num_replicates = NUMTRIALS)

def getvalue(N_recent, N_old,  SPLIT):
	trees = getTrees(N_recent, N_old,  SPLIT)
	BranchLengths = 0
	Expected = [0] * numhaps

	for ts in trees:
		BranchLengths = BranchLengths + ts.first().total_branch_length
		Expected[1:] = Expected[1:] + ts.allele_frequency_spectrum(mode='branch', polarised=True)[1:numhaps]

	Expected[0] = 1 - mu * (BranchLengths / NUMTRIALS)
	Expected[1:]  = (Expected[1:] / np.sum(Expected[1:]) ) * (1 - Expected[0])
	return(kl(Expected, empirical))

def negloglieklihoodfunction(args):
	recentne = 10**args[0]
	oldne = 10**args[1]
	gensplit = 10**args[2]
	A = getvalue(recentne,oldne ,gensplit)
	print(recentne, oldne, gensplit, A)

	return(A)


#[436725.41689147 122360.15211858  16474.75155179] 1.3239340884450319e-06   IN HAPLOID UNITS!
#[218362 61180 16474.75155179] 1.3239340884450319e-06   IN diploid UNITS!

res = minimize(negloglieklihoodfunction, [4, 3.5, 3.5],
            method='Nelder-Mead', options = {"xatol" : 1e-6}).x

print(res)
