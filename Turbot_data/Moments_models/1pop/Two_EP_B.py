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
from numpy import exp
from numpy import log
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys
infile=sys.argv[1]
pop_ids=[sys.argv[2]]
projections=[int(sys.argv[3])]
params=[1,1,1,1,1]

# mutation rate per sequenced portion of genome per generation
mu=float(sys.argv[4])
# generation time, in thousand years
gtime=float(sys.argv[5]) 



dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     

def Two_EP_B(params, ns):
    """
    There are no additional parameter, this is a standard neutral model with no pop changes 
    """
    nu1,nu2,nu3,T1,T2 = params

    nu1_func = lambda t: nu2 * (nu3/nu2)**(t/T2)
    nu_func = lambda t: [nu1_func(t)]


    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    fs.integrate([nu1], T1)
    fs.integrate(nu_func, T2)
    
    return fs

func=Two_EP_B


upper_bound = [100,100,100,10,10]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3]
params = moments.Misc.perturb_params(params, fold=int(sys.argv[6]), upper_bound=upper_bound,
                              lower_bound=lower_bound)

# fitting (poptg = optimal parameters):
poptg = moments.Inference.optimize_log(params, data, func,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=False, maxiter=int(sys.argv[7]))
# extracting model predictions, likelihood and theta
model = func(poptg, ns)
ll_model = moments.Inference.ll_multinom(model, data)
theta = moments.Inference.optimal_sfs_scaling(model, data)

# random index for this replicate
ind=str(random.randint(0,999999))

                                   
# printing parameters 
print "RESULT","Two_EP_B",ind,len(params),ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta



# plotting quad-panel figure witt AFS, model, residuals:
moments.Plotting.plot_1d_comp_multinom(model, data,  residual='Anscombe',plot_masked=False)
plt.savefig("Two_EP_B_"+ind+"_"+sys.argv[1]+"_"+sys.argv[2]+"_"+sys.argv[3]+'.pdf')

# plotting demographic model
plot_mod = moments.ModelPlot.generate_model(func, poptg, ns)
moments.ModelPlot.plot_model(plot_mod, save_file="Two_EP_B_"+ind+"_"+sys.argv[1]+".png",pop_labels=pop_ids, nref=theta/(4*mu), draw_scale=False, gen_time=gtime, gen_time_units="KY", reverse_timeline=True)

