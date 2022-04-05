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

@pytest.fixture
def run_cli(testdir, filepath_tests):
   def do_run(ex,*args):
       args = ['python', filepath_tests + '/../tools/'+ex] + list(args)
       return testdir.run(*args)
   return do_run

def get_atraj_data(binf):
   import pyanalisi as pa
   trajg=pa.Traj(binf)
   trajg.setWrapPbc(False)
   trajg.setAccessWindowSize(trajg.getNtimesteps())
   trajg.setAccessStart(0)
   return trajg.get_box_copy(), trajg.get_positions_copy(), trajg.get_velocities_copy()

#./cp2analisi.py --nstep 4 -p cp/water125_testanalisi -s cp/water125_testanalisi.species --natoms 375
def test_cp2analisi(tmpdir_factory, run_cli,filepath_tests, num_regression):
   import pyanalisi as pa
   inputf=filepath_tests+ '/../tools/cp/water125_testanalisi'
   outf=inputf+'.bin'
   output = run_cli('cp2analisi.py', '-p', inputf, '-s', filepath_tests+ '/../tools/cp/water125_testanalisi.species', '--natoms', '375', '--nstep', '4')
   if pa.has_mmap():
      box,pos,vel = get_atraj_data(outf)
      num_regression.check({'pos':pos.flatten(), 'box':box.flatten(), 'vel':vel.flatten()})
   else:
      warnings.warn(UserWarning('This compiled version of pyanalisi DOES NOT SUPPORT READING LAMMPS FILES'))
      
def test_lammps2analisi(tmpdir_factory, run_cli,filepath_tests, num_regression):
   import pyanalisi as pa
   inputf=filepath_tests+ '/../tools/lammps/dump.lammpstrj'
   outf=inputf+'.bin'
   output = run_cli('lammps2analisi.py', '-p', inputf, '--nstep', '4')
   if pa.has_mmap():
      box,pos,vel = get_atraj_data(outf)
      num_regression.check({'pos':pos.flatten(), 'box':box.flatten(), 'vel':vel.flatten()})
   else:
      warnings.warn(UserWarning('This compiled version of pyanalisi DOES NOT SUPPORT READING LAMMPS FILES'))
      
   


   
