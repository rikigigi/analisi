import pytest
import numpy as np

def test_write_binary_header(analisi_traj, numpy_traj, tmpdir_factory):
    import pyanalisi
    fn = tmpdir_factory.mktemp("data").join("binary_test.lammps")
    tmpfile = str(fn)
    pos,vel,types,box = numpy_traj
    analisi_traj.write_lammps_binary(tmpfile,0,-1)
    ltraj = pyanalisi.Traj(tmpfile)
    ltraj.setAccessWindowSize(3)
    ltraj.setAccessStart(0)
    assert (np.array(ltraj)[:3]==pos[:3]).all()
    ltraj.toggle_pos_vel()
    assert (np.array(ltraj)[:3]==vel[:3]).all()
