#!/usr/bin/env python

# constant pop size
# n(para): 5

import matplotlib
matplotlib.use('PDF')
import moments
import random
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
pop_ids=[sys.argv[2]]
projections=[int(sys.argv[3])]


# mutation rate per sequenced portion of genome per generation
mu=float(sys.argv[4])
# generation time, in thousand years
gtime=float(sys.argv[5]) 



dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

def Snm(params, ns):
    """
    There are no additional parameter, this is a standard neutral model with no pop changes 
    """

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    
    return fs

func=Snm


# extracting model predictions, likelihood and theta
neutral_params = array([])
neutral_upper_bound = []
model = func(neutral_params, ns)


ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

                                   
# printing parameters 
print "RESULT","Snm",ind,ll_model,sys.argv[1],sys.argv[2],sys.argv[3],theta



# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_1d_comp_multinom(model, data,  residual='Anscombe',plot_masked=False)
plt.savefig("Snm_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+'.pdf')

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func,neutral_params, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="Snm_"+ind+"_"+sys.argv[1]+".png",pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

