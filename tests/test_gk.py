import pytest
import numpy as np

def test_gk(analisi_log,num_regression):
    import pyanalisi
    gk=pyanalisi.GreenKubo(analisi_log,'',1,['c_flux[1]','c_vcm[1][1]'], False, 2000, 2,False,0,4,False,1,100)
    gk.reset(analisi_log.getNtimesteps()-2000)
    gk.calculate(0)
    m=np.array(gk,copy=False)
    num_regression.check({'gk':m.flatten()})
