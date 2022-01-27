import pytest
import warnings
import numpy as np

def test_write_binary_header(analisi_traj, numpy_traj, tmpdir_factory):
    import pyanalisi
    fn = tmpdir_factory.mktemp("data").join("binary_test.lammps")
    tmpfile = str(fn)
    pos,vel,types,box = numpy_traj
    analisi_traj.write_lammps_binary(tmpfile,0,-1)
    nt = pos.shape[0]
    if pyanalisi.has_mmap():
       ltraj = pyanalisi.Traj(tmpfile)
       ltraj.setAccessWindowSize(nt)
       ltraj.setAccessStart(0)
       assert (np.array(ltraj)[:nt]==pos[:nt]).all()
       ltraj.toggle_pos_vel()
       assert (np.array(ltraj)[:nt]==vel[:nt]).all()
       assert np.allclose(analisi_traj.get_box_copy()[:nt], ltraj.get_box_copy()[:nt], rtol=1e-10, atol=1e-10, equal_nan=False) 
    else:
       warnings.warn(UserWarning('This compiled version of pyanalisi DOES NOT SUPPORT READING LAMMPS FILES'))

def test_write_binary_header_triclinic(triclinic_traj, tmpdir_factory):
    import pyanalisi
    fn = tmpdir_factory.mktemp("data").join("binary_test_tri.lammps")
    tmpfile = str(fn)
    pos = triclinic_traj.get_positions_copy()
    vel = triclinic_traj.get_velocities_copy()
    box = triclinic_traj.get_box_copy()
    triclinic_traj.write_lammps_binary(tmpfile,0,-1)
    nt = pos.shape[0]
    if pyanalisi.has_mmap():
       ltraj = pyanalisi.Traj(tmpfile)
       ltraj.setAccessWindowSize(nt)
       ltraj.setAccessStart(0)
       assert (np.array(ltraj)[:nt]==pos[:nt]).all()
       ltraj.toggle_pos_vel()
       assert (np.array(ltraj)[:nt]==vel[:nt]).all()
       assert np.allclose(box[:nt], ltraj.get_box_copy()[:nt], rtol=1e-10, atol=1e-10, equal_nan=False)
    else:
       warnings.warn(UserWarning('This compiled version of pyanalisi DOES NOT SUPPORT READING LAMMPS FILES'))


