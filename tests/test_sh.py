import pytest
import numpy as np

def test_sh(analisi_traj,num_regression):
    import pyanalisi
    startr=0.7
    endr=3.5
    nbin=4
    dr=(endr-startr)/nbin
    test=pyanalisi.SphericalCorrelations(analisi_traj
                                         ,startr #minimum of radial distance
                                         ,endr #maximum of radial distance
                                         ,nbin # number of distance bins
                                         ,200 # maximum length of correlation functions in timesteps
                                         ,4 # number of threads
                                         ,90 # time skip to average over the trajectory
                                         ,10 #buffer size
                                         ,False)
    test.reset(700) # number of timesteps to analyze
    test.calculate(0) #calculate starting at this timestep
    res=np.array(test,copy=True) # get result
    print(res)
    num_regression.check({
            'spherical_harmonics_correlation' : res.flatten()
            })
    
