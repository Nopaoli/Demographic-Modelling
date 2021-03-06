#!/usr/bin/env python

# IM model: split in two arbitrary  Ne,  asymmetric migration , hetero Ne and hetero m
# n(para): 10

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
params=[1,1,1,1,1,0.1,0.1,0.1,0.1,0.1]




dd = Misc.make_data_dict(infile)
data = Spectrum.from_data_dict(dd, pop_ids,projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)     


#-------------------

def IM_2M_2N(params, ns):
    """
    Isoation with migration model,  split and asymmetric migration following secondary contact
  	
	nu1= pop size for North Sea and the ancestral pop 
	nu2= pop sie for Baltic Sea 
	= time of population split
	T1= time of split
	m12= migration rate from North Sea to Baltic
	m21= migration rate from Baltic Sea to North Sea
	i1= reduction in migration rate for Islands from the Baltic to North Sea 
	i2= reduction in migration rate for Islands from the North Sea to Baltic
	P= proportion of the genome made up of "islands"
	hrf= Hill-Robertson factor (i.e. average Ne for regions under selection as a proportion of "neutral" Ne)
	Q= proportion of the genome under selection
    """
    nu1,nu2,T1,m12,m21,i1,i2,P,hrf,Q = params

# calculate the spectrum for neutral portion of the genome
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
# calculate the spectrum from genomic islands
    stsi = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsi = moments.Spectrum(stsi)
    fsi = moments.Manips.split_1D_to_2D(fsi, ns[0], ns[1])
    fsi.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12*i2], [m21*i1, 0]]))
    # calculate the spectrum for neutral portion of the genome
    stsn1 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn1 = moments.Spectrum(stsn1)
    fsn1 = moments.Manips.split_1D_to_2D(fsn1, ns[0], ns[1])
    fsn1.integrate([nu1, nu2], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
# calculate the spectrum for lower Ne regions of the genome
    stsn2 = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fsn2 = moments.Spectrum(stsn2)
    fsn2 = moments.Manips.split_1D_to_2D(fsn2, ns[0], ns[1])
    fsn2.integrate([nu1*hrf, nu2*hrf], T1, dt_fac=0.01, m=np.array([[0, m12], [m21, 0]]))
    fs2=P*fsi+(1-P)*fs+Q*fsn2+(1-Q)*fsn1
	
    return fs2

func=IM_2M_2N
upper_bound = [100,100,10,200,200,0.999,0.999,0.499,0.999,0.999]
lower_bound = [1e-3,1e-3,1e-3,1e-5,1e-5,1e-3,1e-3,1e-3,1e-3,1e-3]
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
print "RESULT","IM_2M_2N",ind,len(params),opti,round,ll_model,sys.argv[1],sys.argv[2],sys.argv[3],poptg,theta
                                    
