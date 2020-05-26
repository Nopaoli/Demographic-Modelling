import sys
import os
import numpy
from numpy import array
import dadi
import pylab
from datetime import datetime
import Optimize_Functions
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration,Demographics1D,Demographics2D

# infile, in this case simulation from Hudson's MS

infile=sys.argv[1]


'''

General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.
 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688
Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. For ddRADseq data where a single SNP is selected per locus, this
 is true, but if SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See the discussion group for more information on this subject. 
Citations:
 If you use these scripts for your work, please cite the following publication:
    Portik, D.M., Leach, A.D., Rivera, D., Blackburn, D.C., Rdel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K.Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 52455263.
    doi: 10.1111/mec.14266
-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-dadi
-------------------------
Daniel Portik
daniel.portik@gmail.com
https://github.com/dportik
Updated May 2018
'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

#**************

data = dadi.Spectrum.from_ms_file (infile, average=False)

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
'''

def SC(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    nu1, nu2, m12, m21, T1, T2 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs




#	Below labels and upper and lower bounds of parameter values


pts = [120,130,140]
p_labels = "nu1, nu2, m12, m21, T1,T2"
upper = [20,20,10,10,5,5]
lower = [1e-3,1e-3,1e-5,1e-5,1e-3,1e-3]

#	Here details of the optimization routines: 3 rounds of optimization for the simple IM and SC models

reps = [10,10,10]
maxiters = [10,20,30]
folds = [3,2,1]

# 	Run 3 Independent optimizations routines, at the end use the "Summarize_Output.py" to keep the best run from each optimization routine

for i in range(1,4):
    prefix = infile+"_OPTI_Number_{}".format(i)
    Optimize_Functions.Optimize_Routine(data, pts, prefix, "SC", SC, 3, 6, data_folded=False, param_labels = p_labels, in_upper=upper, in_lower=lower, reps = reps, maxiters = maxiters, folds = folds, optimizer="log")
