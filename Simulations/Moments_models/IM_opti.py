import sys
import os
import numpy
from numpy import array
import moments
import pylab
from datetime import datetime
import Optimize_Functions
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D

# infile, in this case simulation from Hudson's MS
infile=sys.argv[1]

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************

data = moments.Spectrum.from_ms_file (infile, average=False)

numpy.set_printoptions(precision=3)    
 


#print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
# print "projection", projections
print "sample sizes", data.sample_sizes
sfs_sum = numpy.around(data.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

#================================================================================
# Here is an example of using a custom model within this script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script:
 
 Optimize_Routine(data, outfile, model_name, func, rounds, param_number, fs_folded=True, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")
 
   Mandatory Arguments =
    data:  spectrum object name
    outfile:  prefix for output naming
    model_name: a label to slap on the output files; ex. "no_mig"
    func: access the model function from within 'dadi_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_2D, calling Models_2D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).
   Optional Arguments =
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order	 
	 
	THIS IS A MODIFIED VERSION OF DANIEL PORTIK's PIPELINE, CITE ALSO THIS PAPER:
    Portik, D.M., Leach, A.D., Rivera, D., Blackburn, D.C., Rdel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266
'''



def IM(params, ns):
    """
    nu1= pop size for North Sea
	nu2=pop size for Baltic Sea 
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu1,nu2,T1,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=numpy.array([[0, m12], [m21, 0]]))
    return fs


#	Below labels and upper and lower bounds of parameter values


p_labels = "nu1, nu2, T1, m12, m21"
upper = [20,20,10,200,200]
lower = [1e-3,1e-3,1e-3,1e-5,1e-5]

#	Here details of the optimization routines: 4 rounds of optimization for large dataset and unfolded spectrum

reps = [10,10,20,20]
maxiters = [10,20,20,30]
folds = [3,3,2,1]

# 	Run 10 Independent optimizations routines, at the end use the "Summarize_Output.py" to keep the best run from each optimization routine


for i in range(1,11):
    prefix = infile+"_OPTI_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(data, prefix, "IM", IM, 4, 5, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds)
