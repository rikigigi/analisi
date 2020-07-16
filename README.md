# analisi: your Swiss Army Knife of molecular dynamics analysis

Fastest possible example usage:

## command line

- mean square displacement
```
./analisi -i tests/lammps.bin -Q > output_file
```

- spherical harmonics correlation functions
```
./analisi -i tests/lammps.bin -Y 4 -F 0.7 3.5 > output_file
```

- g(r), with time lags too
```
./analisi -i tests/lammps.bin -g 200 -F 0.7 3.5 > output_file
```

- green kubo autocorrelation integrals 
```
./analisi -l tests/gk_integral.txt -H -a 'c_flux[1]' 'c_vcm[1][1]' > output_file
```

many others...

## python

```
#read trajectory
import numpy as np
pos = np.load( 'tests/data/positions.npy')
vel = np.load( 'tests/data/velocities.npy')
box = np.load( 'tests/data/cells.npy')
types = np.zeros(pos.shape[1], dtype = np.int32)
types[-16:-8]=1
types[-8:]=2
print('position array shape is {}'.format(pos.shape))
print('first cell matrix is {}'.format(box[0]))

#create trajectory object
import pyanalisi
analisi_traj = pyanalisi.Trajectory(pos, vel, types, box,True, False)

#do the calculation that you want
msd=pyanalisi.MeanSquareDisplacement(analisi_traj,10,4000,4,True,False,False)
msd.reset(3000)
msd.calculate(0)

result=np.array(msd,copy=False)


#other calculations
#...
#...

```


## note

If you use this program and you like it, spread the word around and give it credits! Implementing stuff that is already implemented can be a waste of time. And why don't you try to implement something that you like inside it? You will get for free MPI parallelization and variance calculation, that are already implemented in a very generic and abstracted way.

# Description


Features:

 - python interface (reads numpy array)
 - command line interface (reads binary lammps-like files)
 - multithreaded
 - command line interface has MPI too (for super-heavy calculations)
 - command line calculates variance of every quantity and every function (in python you can do it by yourselves with numpy.var )

Calculations:

- g of r ( and time too!)
- vibrational spectrum (this is nothing special)
- histogram of number of neighbours
- mean square displacement
- green kubo integral of currents
- multicomponent green kubo time domain formula (MPI here can be useful...)
- spherical harmonics number density time correlation analysis
- and more ...
- ...

# Building from source
Dependencies:

- C++17 capable compiler
- linux (mmap)
- FFTW3 (not included)
- Eigen3 (included in the package)
- Boost (included in the package)
- Mpi (optional)
- libxdrfile (for gromacs file conversion -- optional)
- python (optional) 

## MPI build (why not?)
```
mkdir build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DUSE_MPI=ON
make
```
## non-MPI build (shame on you!)
```
mkdir build
cd build
cmake ../
make
```

# Credits
Written by Riccardo Bertossa during his lifetime at SISSA
