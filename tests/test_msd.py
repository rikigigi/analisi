import pytest
import pyanalisi
import numpy as np

def test_gofr(analisi_traj,num_regression):
    msd=pyanalisi.MeanSquareDisplacement(analisi_traj,10,4000,4,True,False,False)
    msd.reset(3000)
    msd.calculate(0)
    m=np.array(msd,copy=False)
    print(m)
    num_regression.check({'msd':m.flatten()})

