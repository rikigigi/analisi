import pytest
import pyanalisi
import numpy as np

def test_gofr(analisi_traj,num_regression):
    startr=0.0
    endr=3.8
    nbin=200
    gofr=pyanalisi.Gofrt(analisi_traj,startr,endr,nbin,10,4,10,False)
    gofr.reset(700)
    gofr.calculate(0)
    gr=np.array(gofr,copy=True)
    print(gr)
    num_regression.check({'g_of_r_and_t':gr.flatten()})
