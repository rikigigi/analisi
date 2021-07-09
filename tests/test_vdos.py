import pytest
import numpy as np

def test_vdos(analisi_traj,num_regression):
    import pyanalisi
    vdos=pyanalisi.VibrationSpectrum(analisi_traj,False)
    vdos.reset(analisi_traj.getNtimesteps())
    vdos.calculate(0)
    m=np.array(vdos)
    num_regression.check({'vdos':m.flatten()})
