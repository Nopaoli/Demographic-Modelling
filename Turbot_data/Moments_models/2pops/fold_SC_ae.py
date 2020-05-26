#!/usr/bin/env python

# secondary contact with ancestral population expansion: ancestral expansion, split in two arbitrary constant Ne, asymmetric migration following secondary contact
# n(para): 6

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
pop_ids=[sys.argv[2],sys.argv[3]]
projections=[int(sys.argv[4]),int(sys.argv[5])]
params=[1,1,1,1,1,1,1,1]



dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     


#-------------------
# split with growth and asymmetrical migration
def SC_ae(params, ns):
    """
    nu1= pop size following ancestral population expnasion (assumed constant in both NS and Baltic Sea)
	Tae= timing of ancestral population expansion
	T1= time of population split
	T2= time of secondary contact
	m12_0= migration rate from North Sea to Baltic
	m21_0= migration rate from Baltic Sea to North Sea
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
    """
    nu_ae,nu1,nu2,Tae,T1,T2,m12,m21 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, m = np.array([[0, 0], [0, 0]]))
    fs.integrate([nu1, nu2], T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    return fs

func=SC_ae
upper_bound = [100,100,100,10,10,10,200,200]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5]
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
# optimization number
opti=int(sys.argv[8])
# round number
round=(sys.argv[9])



# printing parameters 
print "RESULT","SC_ae",ind,len(params),opti,round,ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta
                                    
