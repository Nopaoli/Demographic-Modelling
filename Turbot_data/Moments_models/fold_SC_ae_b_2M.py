#!/usr/bin/env python

# secondary contact: ancestral pop expansion,split in two arbitrary Ne, asymmetric migration following secondary contact and growth of the Baltic Sea pop , heterogeneous migration rate
# n(para): 11

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
params=[1,1,1,0.1,1,1,1,1,1,0.1,0.1,0.1]




dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     


#-------------------
def SC_ae_b_2M(params, ns):
    """
    nu1= pop size after pop expansion (and pop size of North Sea)
	s=proportion of the North Sea pop which invaded the Baltic (i.e. original bottleneck)
	nu2= final size of Baltic Sea pop
	T1= time of population split
	T2= time of second epoch of migration and population growth
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
    """
    nu_ae,nu1,nu2,s,Tae,T1,T2,m12,m21,i1,i2,P = params
    nu2_0 = nu1*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T2)
    nu_func= lambda t: [nu1,nu2_func(t)]
# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs.integrate([nu_ae], Tae)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu1*s], T1, m = np.array([[0, 0], [0, 0]]))
    fs.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
# calculate teh spectrum for genomic islands	
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi.integrate([nu_ae], Tae)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu1*s], T1, m = np.array([[0, 0], [0, 0]]))
    fsi.integrate(nu_func, T2, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
	
    fs2=P*fsi+(1-P)*fs
    return fs2

func=SC_ae_b_2M
upper_bound = [100,100,100,0.999,10,10,10,200,200,0.999,0.999,0.499]
lower_bound = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3]
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
print "RESULT","SC_ae_b_2M",ind,len(params),opti,round,ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta
                                    
