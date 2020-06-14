# analisi: your Swiss Army Knife of molecular dynamics analysis

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

- C++17 
- linux (mmap)
- FFTW3
- Eigen3
- Boost (program_options)
- Mpi (optional)
- libxdrfile (for gromacs -- optional)
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
