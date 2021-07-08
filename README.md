# analisi: your Swiss Army Knife of molecular dynamics analysis


[Link to the pdf version of this document](README.pdf)

Binary packages available on linux and macOS [conda-forge](https://conda-forge.org/#about):
```
conda config --add channels conda-forge    #if not already done
conda config --set channel_priority strict #if not already done
conda install analisi
```

| Name | Downloads | Version | Platforms |
| --- | --- | --- | --- |
| [![Conda Recipe](https://img.shields.io/badge/recipe-analisi-green.svg)](https://anaconda.org/conda-forge/analisi) | [![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/analisi.svg)](https://anaconda.org/conda-forge/analisi) | [![Conda Version](https://img.shields.io/conda/vn/conda-forge/analisi.svg)](https://anaconda.org/conda-forge/analisi) | [![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/analisi.svg)](https://anaconda.org/conda-forge/analisi) |

or by compiling from source following this [section](#building-from-source)

Fastest possible example usage:

## command line example

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


and many others...

## python example

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

Heat transport coefficient calculation: correlation functions and gk integral for a multicomponent fluid example
```
import numpy as np
with open(filepath_tests + '/data/gk_integral.dat', 'r') as flog:
    headers = flog.readline().split()
log = np.loadtxt(filepath_tests + '/data/gk_integral.dat', skiprows=1)

import pyanalisi
traj = pyanalisi.ReadLog(log, headers)
gk = pyanalisi.GreenKubo(analisi_log,'',
       1,['c_flux[1]','c_vcm[1][1]'],
       False, 2000, 2,False,0,4,False,1,100)
gk.reset(analisi_log.getNtimesteps()-2000)
gk.calculate(0)
result = np.array(gk,copy=False)

```

## note

If you use this program and you like it, spread the word around and give it credits! Implementing stuff that is already implemented can be a waste of time. And why don't you try to implement something that you like inside it? You will get for free MPI parallelization and variance calculation, that are already implemented in a very generic and abstracted way.

# Description

This is a framework for computing averages on molecular dynamics trajectories and block averages.
An mpi parallelization over the blocks is implemented in a very abstract and general way, so it is easy to implement a new calculation that can be parallelized with no effort using MPI. The code has two interfaces: a command line interface that is parallelized with MPI and threads, and a python interface that is parallelized on threads only. With such many parallelization levels more accurate averages can be performed, as described in details in the next sections.

Features:

 - python interface (reads numpy array)
 - command line interface (reads binary lammps-like files or time series in column format)
 - multithreaded ( defaults to number of threads specifies by the shell variable OMP_NUM_THREADS )
 - command line interface has MPI too (for super-heavy calculations)
 - command line calculates variance of every quantity and every function (in python you can do it by yourselves with numpy.var )
 - jupyter notebook example of the python interface

Calculations:

- g of r ( and time too!)
- vibrational spectrum (this is nothing special)
- histogram of number of neighbours
- mean square displacement
- green kubo integral of currents
- multicomponent green kubo time domain formula (MPI here can be useful...)
- spherical harmonics number density time correlation analysis (MPI here can be useful...)
- atomic position histogram
- and more ...
- ...

# Building from source
Dependencies:

- C++17 compatible compiler (GCC 7+, clang 6+, intel 19.0.1+, [source](https://en.cppreference.com/w/cpp/compiler_support) )
- cmake
- linux or macOS (mmap - maybe this dependency can be removed with some additional work)
- FFTW3 (included in the package)
- Eigen3 (included in the package)
- Boost (included in the package)
- Mpi (optional)
- libxdrfile (for gromacs file conversion -- included in the package)
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

After compiling the code you will find the executable `analisi` and the shared library `pyanalisi.*.so` (with a  part of the name that depends on your particular python installation) in the building folder. To install the python library, simply copy it in your site-library folder, for example by executing the command
```
cp pyanalisi.*.so "`python3 -m site --user-site`/pyanalisi.so"
```
Be careful to choose the correct python executable, the same that you used to compile the library.
This task can be performed as well by the command `make install/local` inside a python virtual environment, for example a miniconda installation.

To test that the library was compiled correctly and that there are no regressions, you can run the (small) test suite with the command
```
make test
```

## additional cmake options

 - `-DBUILD_TESTS=OFF` skips the building of the C++ tests
 - `-DSYSTEM_FFTW3=ON` use system's FFTW3 library
 - `-DSYSTEM_EIGEN3=ON` use installed system's eigen3 library
 - `-DSYSTEM_BOOST=ON` use detected system's boost library
 - `-DSYSTEM_XDRFILE=ON` use detected system's xdrfile library
 - `-DPYTHON_EXECUTABLE=/path/to/python` use this python installation to build the python interface
 - `-DPYTHON_INTERFACE=OFF` don't build any python interface (the building process will not depend on python)

# Documentation
This document is better rendered in the pdf version. [Link to the pdf version](README.pdf). Note that is generated by the script `build_pdf.sh`, and the original file is `README_.md`.

## Command line interface

The command line utility is able to read only binary trajectory files in the LAMMPS format, specified with the command line option <code>-i <i>input_file</i></code>, or a time series in a column formatted text file with a header, specified with the command line option <code>-l <i>input_file</i></code>. The trajectory file can be generated in many ways:
 - by LAMMPS :-)
 - by using the command line utility with the command line options <code>-i <i>input_file</i> -binary-convert <i>output_file</i></code>, where <code><i>input_file</i></code> is the name of a plain text trajectory in the format:
 
   ```
   [natoms]
   [xlo] [xhi]
   [ylo] [yhi]
   [zlo] [zhi]
   [id_1] [type_1] [x_1] [y_1] [z_1] [vx_1] [vy_1] [vz_1]
    .      .       .    .    .    .     .     .
    .      .       .    .    .    .     .     .
    .      .       .    .    .    .     .     .
   [id_natoms] [type_natoms] [x_natoms] [y_natoms] [z_natoms] [vx_natoms] [vy_natoms] [vz_natoms]
   .
   .
   .
   ```
   
   That is: for every step you have to provide the number of atoms, low and high coordinates of the orthorombic cell, and then for every atom its id, type id, positions and velocities.
 - by using the command line utility with the command line options <code>-i <i>input_file</i> -binary-convert-gromacs <i>output_file</i> <i>typefile</i></code> and a gromacs trajectory (you have to provide the xdr library in the building process)
 - by using the python interface: see section [Buffer protocol interface](#Buffer-protocol-interface)
 
### Command line interface common arguments

Where meaningful, the following arguments are shared by all calculation performed by the command line utility
 - <code>-i <i>input_file</i></code> specify input trajectory binary file. See [lammps documentation](https://lammps.sandia.gov/doc/dump.html) for the documentation of how to produce a binary dump from lammps. The order of the dumped atomic quantities **must be** `id type xu yu zu vx vy vz`
 - <code>-l <i>input_file</i></code> specify input time series in text column format with headers
 - <code>-N <i>thread_number</i></code> specify the number of threads to use where parallelization is implemented with threads
 - <code>-B <i>number_of_block</i></code> when the code calculates variances of the quantity, it splits the trajectory in many blocks. With this option you can specify the number of blocks. If MPI parallelization is used, since different blocks are calculated in parallel, you may want to spacify a number of blocks that is a multiple of the number of MPI processes. In this way, at the final iteration over the blocks, all processes will be busy.
 - <code>-s <i>timestep_number_skip</i></code> usually neighbour configurations in the trajectory are strongly correlated. With this option you specify how far must be two consecutive configuration used to calculate averages.
 - <code>-S <i>timestep_number_length</i></code> specify how long must be the function that the code calculates, measured in number of timesteps. (for example the length of the MSD plot)

## Python interface

### Creating a trajectory object

You can create a trajectory python object to be used in this library in two ways: 
 - start from python arrays generated by other code
 - use a LAMMPS binary file that you have on the filesystem. This is the same file that is used by the command line interface.
 

#### using Python arrays: the buffer protocol interface


You must have 4 arrays. In general the interface supports any object that supports python's buffer protocol. This include, among others, numpy arrays. Suppose you have a trajectory of `t` timesteps and `n` atoms. You need the following arrays:
 - position array, of shape `(t,n,3)`
 - velocity array, of shape `(t,n,3)`
 - cell array, of shape `(t,3,3)` only diagonal matrices (orthorombic cells) are supported at the moment, or if a lammps formated cell is provided `(t,6)` The lammps format simply list the low and high coordinates of the orthorombic cell, as in the header of each timestep that you find in the lammps trajetory format (shown in [command line interface](#command-line-interface) section )
 - types array, of shape `(n)` (integer array)
 
In general no particular units of measure are required, and the output will reflect the units that you used in the input. The calculations that the program does are reported exactly in the following sections. From those formulas you can deduce the units of the output given the units that you used in the input.

Then you must decide if you want that the coordinates are rewrapped inside the cell or not. At the moment only orthorombic cells are supported in all calculations but those that need only unwrapped coordinates, like MSD.

The lammps format for the cell is `[x_lo, x_hi, y_lo, y_hi, z_lo, z_hi]`: you have to provide only the coordinates boundaries.

The syntax for creating the object is 
```
import pyanalisi
analisi_traj = pyanalisi.Trajectory(positions_array,
                                    velocity_array,
                                    types_array,
                                    box_array,
                                    use_matrix_or_lammps_cell_format,
                                    wrap_atomic_coordinates)
```
where `use_matrix_or_lammps_cell_format` is `True` if usual matrix format for the cell is given and `False` if a lammps formatted cell is provided and `wrap_atomic_coordinates` is `True` if you want to wrap all the atomic coordinates inside the cell.

You can write a LAMMPS bynary trajectory (that can be used by the command line interface with mpi, for example) by calling
```
analisi_traj.write_lammps_binary('output_path', start_timestep, end_timestep)
```
where `start_timestep` is the first timestep index to dump (indexes start from 0) and  `end_timestep` is the first timestep that will not be written. If `end_timestep == -1`, it will dump everything till the end of the trajectory. This is a very convenient way of moving heavy computations on a cluster where MPI can be used, or more in general to convert a generic trajectory format in the format used by the command line tool. For example

   ```
   #read trajectory. It can come from everywhere
    import numpy as np
    pos = np.load( 'tests/data/positions.npy') #shape (N_timesteps, N_atoms, 3)
    vel = np.load( 'tests/data/velocities.npy') #shape (N_timesteps, N_atoms, 3)
    box = np.load( 'tests/data/cells.npy') #shape (N_timesteps, 3, 3)
    types = np.load( 'tests/data/types.npy') #shape (N_atoms), dtype=np.int32
    
    #create trajectory object and dump to file
    import pyanalisi
    analisi_traj = pyanalisi.Trajectory(pos, vel, types, box,
                                        True, # matrix format for the box array
                                        False # don't wrap the coordinates
                                        )
    analisi_traj.write_lammps_binary("output_filename.bin"
                               , 0, # starting timestep
                               -1   # last timestep:
                               )    #  -1 dumps full trajectory
   ```


#### LAMMPS binary trajectory interface


This interface is a little more complicated, since it was designed for computing block averages of very big files. It can read the same files that the command line program reads. The object is created with
```
analisi_traj = pyanalisi.Traj('path_of_binary_file')
```

Then you have to call some more functions, BEFORE calling the compute routines :
```
analisi_traj.setWrapPbc(True) #optional: if you want to wrap positions inside the cell
analisi_traj.setAccessWindowSize(size_in_number_of_steps_of_the_reading_block)
analisi_traj.setAccessStart(first_timestep_to_read)
```
The first line is to set the pbc wrapping (don't use this if you have to compute the MSD!). Since the access to the trajectory is done in blocks, the second line sets the size of the reading block. The bigger the block the bigger the memory allocated to store the positions and the velocities. The last call sets the index of the first timestep that is going to be read, and then read it. In this moment the selected chunk of the trajectory is loaded into your memory. At this point you are able to call the needed compute routine, making sure that it is not going to read past the last timestep of the block. Later you can call again `setAccessStart` to load a different part of the trajectory and compute the quantity again. Only the data of the current block is stored in the memory, and if evenutally there is some overlap with the previous block, the data is copied from memory and the overlapped part is not read again from the filesystem.


### Creating a time series object

Before analyzing a time series, for example to calculate the integral of the autocorrelation function, it is necessary to create a time series object. This object will hold the data that a different function can analyze. It is created with the following:
```
time_series = pyanalisi.ReadLog(data_array, headers_array)
```
where `data_array` and `header_array` are python objects that support the buffer protocol interface, for example numpy arrays, or a plain python list.
The `data_array` is a two dimensional array where the first index is the timestep index and the second index is the column index.
The `header_array` object must describe the columns that are stored in the bidimensional floating point array `data_array` using Python Strings.


## MSD

Given a trajectory <img src="svgs/995aad7488ee6146d6b411759a5db052.svg?invert_in_darkmode" align=middle width=20.789844899999988pt height=27.91243950000002pt/> where <img src="svgs/0e45672d84ee21e1464939301cc46b43.svg?invert_in_darkmode" align=middle width=137.42503499999998pt height=24.65753399999998pt/> is the atomic index and <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/> is the timestep index, defining the center of mass position of the atomic species <img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> at the timestep <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/> as
<p align="center"><img src="svgs/1069af11f835cf7996677118283e621d.svg?invert_in_darkmode" align=middle width=165.4146549pt height=51.449925449999995pt/></p>
where <img src="svgs/bf032e59612f5135c5425a8a0c07f2f7.svg?invert_in_darkmode" align=middle width=19.312276499999992pt height=22.465723500000017pt/> is the number of atoms of the specie <img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/>,
the code computes the following
<p align="center"><img src="svgs/fb6fd8b28a21787d6b119a3f059a9bd3.svg?invert_in_darkmode" align=middle width=429.54891045pt height=110.10997139999999pt/></p>
If the option `--mean-square-displacement-self` is provided in the command line or in the python interface the documented argument is set to `True`, the atomic mean square displacement for each atomic species is calculated in the reference system of the center of mass of that particular atomic specie. That is, in this case the following is computed:
<p align="center"><img src="svgs/89ebf178c4461aec4cd4651cafb867ab.svg?invert_in_darkmode" align=middle width=627.5146679999999pt height=110.10997139999999pt/></p>
In the output you have many columns, one for each of the <img src="svgs/1a60f096b735c3e3bae490620b576fc8.svg?invert_in_darkmode" align=middle width=44.47093199999999pt height=22.465723500000017pt/> atomic species, first the block of the atomic MSD and then eventually the block of the center of mass MSD if asked to compute. The center of mass MSD is computed only if the command line option `-Q` is provided or the documented argument is set to `True` in the python interface constructor. The output is the following:
<p align="center"><img src="svgs/cd779299413189940523acc73b2c31da.svg?invert_in_darkmode" align=middle width=406.26059265pt height=20.87668275pt/></p>
In the command line output after each column printed as described in the line above you will find the variance calculated with a block average over the specified number of blocks.

The options that you can use for this calculation, in summary, are:

<code>
./analisi -i <i>lammp_binary</i> [-N <i>number of threads</i>] [-S <i>stop timestep</i>] [-s <i>skip every</i>] [-B <i>number of blocks</i>] -q | -Q [--mean-square-displacement-self]
</code>

where any option of `-q`, `-Q`, `--mean-square-displacement-self` activates the calculation of the MSD.

To use this calculation in the python interface, you have to first create a trajectory object (that can be both from a lammps binary or from python arrays objects - described in the section [creating a trajectory object](#creating-a-trajectory-object)), then you have to create the following instance if the numnpy array trajectory is used:
```
msd_calculation = pyanalisi.MeanSquareDisplacement(
    trajectory_instance, # Trajectory instance
    skip, # in calculating averages skip this amount of timestep,
          # as in -s option
    tmax, # size of the MSD(t) plot in number of timesteps
    nthreads, # number of parallel threads to use in the computation
    center_of_mass, # if True calculates center of mass displacement too
    use_cm_reference, # use the reference system of the center of mass
                      # of all atoms of the same type
    debug_flag # if True produces some dump files...
   )

```
if the lammps binary trajectory instance is used, you must use the `pyanalisi.MeanSquareDisplacement_lammps`, with exactly the same arguments. The calculation will be exactly the same. After initializing the object, the calculation can be initialized with
```
msd_calculation.reset(block_size)
```
that sets the number of steps that will be used to calculate the average over the trajectory. If  `tmax` is greater than `block_size`, then `tmax` will be setted to `block_size`. The calculation is started with:
```
msd_calculation.calculate(first_timestep)
```
where `first_timestep` is the index of the first timestep that will be used to calculate the averages. After the calculation over this block is finished, the result can be collected in a numpy array with:
```
result = np.array(msd_calculation, copy=False)
```
In general the object msd_calculation supports the buffer protocol interface, so to get the result you can use anything that can interface with the buffer.

## Green-Kubo

 Given <img src="svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.73973739999999pt height=22.465723500000017pt/> vector time series of length <img src="svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.99998994999999pt height=22.465723500000017pt/> <img src="svgs/6ea8ff937f14aa2a3e0622a395dee889.svg?invert_in_darkmode" align=middle width=27.22419314999999pt height=22.55708729999998pt/>, <img src="svgs/e17d8a1279f0bfa5b92c35f05c242243.svg?invert_in_darkmode" align=middle width=101.57889059999998pt height=24.65753399999998pt/>, <img src="svgs/3ef98fe3db393644cede24e594d556d5.svg?invert_in_darkmode" align=middle width=90.34214144999999pt height=24.65753399999998pt/>,
 implements an expression equivalent to the following formula:
 <p align="center"><img src="svgs/168f1aff83bb74b9c5f70ae822348485.svg?invert_in_darkmode" align=middle width=251.32803135pt height=254.2893243pt/></p>
 but with the trapezoidal rule in place of the sums marked with <img src="svgs/cdcac8939f3840cd8cddf40059a4cf58.svg?invert_in_darkmode" align=middle width=6.735194399999992pt height=22.63846199999998pt/>. Note that <img src="svgs/bcd81920696e26093495d90eb4cd7b1b.svg?invert_in_darkmode" align=middle width=16.153034699999992pt height=22.465723500000017pt/> is a matrix. To get the correct units of measure, you have still to multiply all the quantities but the <img src="svgs/9442d407594839d38f91ad3d2deb134e.svg?invert_in_darkmode" align=middle width=16.71464189999999pt height=22.465723500000017pt/>s by the integration timestep. <img src="svgs/8379d82223552fdbdfd98783823c755b.svg?invert_in_darkmode" align=middle width=33.56332814999999pt height=22.465723500000017pt/> is the number of timesteps on which the code runs the average.
 Every quantity is written in the output in the following order:
 <p align="center"><img src="svgs/9c46ade8e39f77cf4db90ff1e0f72236.svg?invert_in_darkmode" align=middle width=528.55261305pt height=18.7598829pt/></p>
 If the command line tool is used, the variance of the block average is automatically computed, and after each column you find its variance. Moreover you find in the output a useful description of the columns with their indexes.


## g(r,t)

Given a central atom of type <img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> at timestep <img src="svgs/2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/>, calculates the histogram of the minimum image's radial distance of atoms of type <img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> at timestep <img src="svgs/e2c65f38f063936d351de3aa6b418190.svg?invert_in_darkmode" align=middle width=31.25563319999999pt height=22.831056599999986pt/>. The histogram can be of the same atom (you now, it spreads arond while time is passing) or of all atoms different from the original one. Everything is averaged over <img src="svgs/2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> and atoms of the same type. In particular, defining the index of a pair of atoms <img src="svgs/77a3b857d53fb44e33b53e4c8b68351a.svg?invert_in_darkmode" align=middle width=5.663225699999989pt height=21.68300969999999pt/> and <img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/> and a timelag <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/> at a timestep <img src="svgs/2f2322dff5bde89c37bcae4116fe20a8.svg?invert_in_darkmode" align=middle width=5.2283516999999895pt height=22.831056599999986pt/> as
<p align="center"><img src="svgs/2b241d102e574b95a855c128779eafc0.svg?invert_in_darkmode" align=middle width=312.05818545pt height=19.104592649999997pt/></p>
the histogram <img src="svgs/77a68d5c8eea72cbfc44afccff1db2d1.svg?invert_in_darkmode" align=middle width=41.41748159999999pt height=24.65753399999998pt/> between specie <img src="svgs/21fd4e8eecd6bdf1a4d3d6bd1fb8d733.svg?invert_in_darkmode" align=middle width=8.515988249999989pt height=22.465723500000017pt/> and <img src="svgs/8eb543f68dac24748e65e2e4c5fc968c.svg?invert_in_darkmode" align=middle width=10.69635434999999pt height=22.465723500000017pt/> is defined as
<p align="center"><img src="svgs/f34d038041b002c0270cb552072d9188.svg?invert_in_darkmode" align=middle width=261.4717941pt height=40.548151049999994pt/></p>
where <img src="svgs/5f8096da90fe7f742b4614ecdcaf33f6.svg?invert_in_darkmode" align=middle width=37.15183559999999pt height=24.65753399999998pt/> is the Kronecker delta. In practice the program, for each atoms, it adds 1 to the corresponding bin of the histogram. 

For each <img src="svgs/d26f972e4392a750678405391a6a77c7.svg?invert_in_darkmode" align=middle width=39.76018364999999pt height=22.465723500000017pt/>, with <img src="svgs/7e5a4eb9cf494cd7ab77c9102b9e5450.svg?invert_in_darkmode" align=middle width=41.12996909999999pt height=22.465723500000017pt/>, the position of the histograms in the memory is <img src="svgs/5748d094c0cc1c4f6cab2f1906cccb2c.svg?invert_in_darkmode" align=middle width=316.84093485pt height=24.65753399999998pt/> (0 is the first). Given this order of the histograms, the layout in the memory is the following:
<p align="center"><img src="svgs/979acf5a19e3998a4d36c700f47a8a93.svg?invert_in_darkmode" align=middle width=519.2253676500001pt height=19.526994300000002pt/></p>



The layout of the command line output is a gnuplot-friendly one, where the output is organized in blocks, one for each <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/>, separated by two blank lines, and every line corresponds to a histogram bin for every combination of <img src="svgs/36c1aeba203f707a20a8a19715dff5bc.svg?invert_in_darkmode" align=middle width=26.518202699999986pt height=22.465723500000017pt/>:

<p align="center"><img src="svgs/a3a094dfc4263aa79479f7a66c9c8a1a.svg?invert_in_darkmode" align=middle width=306.69572174999996pt height=216.56478195pt/></p>

where we collapsed the indexes <img src="svgs/36c1aeba203f707a20a8a19715dff5bc.svg?invert_in_darkmode" align=middle width=26.518202699999986pt height=22.465723500000017pt/> into a single number, the position order of the histogram. As usual in the command line output every column is an average over all the blocks and it is followed by the variance of the mean.


## Spherical harmonics correlations

### Calculation procedure:

The implemented formula for the real spherical harmonics is the following:
<p align="center"><img src="svgs/0c0d2559a0357c70d0d2e3f54ec979e4.svg?invert_in_darkmode" align=middle width=511.03006845pt height=139.3985439pt/></p>
Where <img src="svgs/16000f6ee652afc9c9c2e4300018fc97.svg?invert_in_darkmode" align=middle width=32.92233779999999pt height=22.831056599999986pt/>, <img src="svgs/2e22eb5d0b0553029f9b63efc7b27688.svg?invert_in_darkmode" align=middle width=33.67577894999999pt height=21.95701200000001pt/>, <img src="svgs/4cc129af69c969f6e93ce3ec90587d2e.svg?invert_in_darkmode" align=middle width=35.50224479999999pt height=14.15524440000002pt/> are calculated using cartesian components:
<p align="center"><img src="svgs/dabaca28551d4a5f720bd5103eb4a5da.svg?invert_in_darkmode" align=middle width=164.5483884pt height=123.47163014999998pt/></p>
and <img src="svgs/b09723ed7a995d2e7a5f91d4b69e91a3.svg?invert_in_darkmode" align=middle width=57.24132809999999pt height=24.65753399999998pt/>, <img src="svgs/948f531991fcec35915fa196d15d6c93.svg?invert_in_darkmode" align=middle width=59.06779229999999pt height=24.65753399999998pt/> are evaluated using Chebyshev polynomials with a recursive definition:
<p align="center"><img src="svgs/16fd8b4eda8830c33d24fbf02c914e8f.svg?invert_in_darkmode" align=middle width=342.84232124999994pt height=49.315569599999996pt/></p>
and <img src="svgs/ef6d9536520a15dcbc2f9d92fbe52e9c.svg?invert_in_darkmode" align=middle width=24.50162219999999pt height=22.465723500000017pt/> are the associated Legendre polynomials, calculated with the following set of recursive definition:
<p align="center"><img src="svgs/8520326519dfeaca6a867a16e6c56fd7.svg?invert_in_darkmode" align=middle width=436.23321884999996pt height=23.1312312pt/></p>
that allows us to calculate every <img src="svgs/4b650efe36e9b5341bf3ca60aff67dcf.svg?invert_in_darkmode" align=middle width=33.79005134999999pt height=24.65753399999998pt/> and every <img src="svgs/7ea990c3d2fab85fabaa40d4f36644d8.svg?invert_in_darkmode" align=middle width=62.10045104999998pt height=24.65753399999998pt/> element of the <img src="svgs/341122b6a0b79a1002ddaadb4be1f5dc.svg?invert_in_darkmode" align=middle width=41.37378464999999pt height=24.65753399999998pt/> values. Then we have the recursion to go up in <img src="svgs/d30a65b936d8007addc9c789d5a7ae49.svg?invert_in_darkmode" align=middle width=6.849367799999992pt height=22.831056599999986pt/>, for any value of it:
<p align="center"><img src="svgs/58038ea2478f4b7921bfee4e82035dbd.svg?invert_in_darkmode" align=middle width=284.74263014999997pt height=37.410748649999995pt/></p>

The program, given a number <img src="svgs/575bfed20f98976cb4eb301ca9ff8312.svg?invert_in_darkmode" align=middle width=33.09897194999999pt height=22.831056599999986pt/> and a triplet <img src="svgs/a35d9ea85439dede6d90c9f53db8be8c.svg?invert_in_darkmode" align=middle width=53.80901294999998pt height=24.65753399999998pt/>, is able to calculate every value of <img src="svgs/3dabeef25aa17a14cbe93b2a1a3372f7.svg?invert_in_darkmode" align=middle width=175.9928742pt height=24.65753399999998pt/> and for all allowed values of <img src="svgs/0e51a2dede42189d77627c4d742822c3.svg?invert_in_darkmode" align=middle width=14.433101099999991pt height=14.15524440000002pt/> with a single recursion. Let
<p align="center"><img src="svgs/ea9309e975a6f8aa8582bb0e7e8f93f9.svg?invert_in_darkmode" align=middle width=296.447415pt height=39.0049506pt/></p>
for some timestep <img src="svgs/4f4f4e395762a3af4575de74c019ebb5.svg?invert_in_darkmode" align=middle width=5.936097749999991pt height=20.221802699999984pt/> in some volume <img src="svgs/82e473be44804d79b6439eb860647d42.svg?invert_in_darkmode" align=middle width=16.30942004999999pt height=22.465723500000017pt/> taken as a the difference of two concentric spheres centered on the atom <img src="svgs/21fd4e8eecd6bdf1a4d3d6bd1fb8d733.svg?invert_in_darkmode" align=middle width=8.515988249999989pt height=22.465723500000017pt/> of radius <img src="svgs/9dd4c0a065698518e075c4a3f2428ba9.svg?invert_in_darkmode" align=middle width=41.01363254999999pt height=14.15524440000002pt/> and <img src="svgs/a6b2ca117874a8471ddab2b11814ebde.svg?invert_in_darkmode" align=middle width=39.33723254999999pt height=14.15524440000002pt/>, and where <img src="svgs/0042e6a739d97acc385edba86cca9bf8.svg?invert_in_darkmode" align=middle width=14.60339264999999pt height=14.15524440000002pt/> is the atomic density of the species <img src="svgs/36b5afebdba34564d884d347484ac0c7.svg?invert_in_darkmode" align=middle width=7.710416999999989pt height=21.68300969999999pt/>. Since the densities are taken as sums of dirac delta functions, it is sufficient to evaluate the spherical harmonics functions at the position of the atoms.
Then we calculate the following:
<p align="center"><img src="svgs/179637dd5a919e5e8305e3b17e5a3f98.svg?invert_in_darkmode" align=middle width=317.23280159999996pt height=50.5480899pt/></p>
where <img src="svgs/1154b2711344812ed220c9c2cb9fdc49.svg?invert_in_darkmode" align=middle width=17.35165739999999pt height=24.65753399999998pt/> is an average operator, and we do an additional average over all the <img src="svgs/40a509fd1a7316439f1607d512f059a5.svg?invert_in_darkmode" align=middle width=21.56622599999999pt height=22.465723500000017pt/> central atoms of the type <img src="svgs/8eb543f68dac24748e65e2e4c5fc968c.svg?invert_in_darkmode" align=middle width=10.69635434999999pt height=22.465723500000017pt/>. The <img src="svgs/1154b2711344812ed220c9c2cb9fdc49.svg?invert_in_darkmode" align=middle width=17.35165739999999pt height=24.65753399999998pt/> average is implemented as an average over the starting timestep.


## Vibrational Spectrum

### Calculation procedure:
The implemented equation is:

<p align="center"><img src="svgs/2bd6ddca99e95c087db7cb6d817b33d9.svg?invert_in_darkmode" align=middle width=309.3120327pt height=47.49793455pt/></p>

where <img src="svgs/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/> is the type of atom.
The diffusivity can be computed as half of the zero value of D.

### command line version
The options that you can use for this calculation are simply:

<code>
./analisi -i <i>lammp_binary</i> [-N <i>number of threads</i>] [-V] [-B <i>number of blocks</i>] 
</code>

the code will generate a file with a number of line equal to the number of frequencies, and with `ntypes_atoms*2+1` columns where:
- the first column rappresent the index of the frequncies 
- then there are the block average and variance of the spectrum for each atomic type

The code is trasparent ot units of measuares of the quantities. If a user wants the diffusivity in the correct units ( e.g. metal) must porcede in the following way:
- the first column can be multiplied by `1/(nstep*dt)` to obtain the frequencies in multiples of Hz
- the other columns can be multiplied by `dt/nstep`;
where `nstep` is the total number of step of the block used to compute  VDOS, `dt` is the time difference of two consecutive molecular dynamis steps.

### python interface

There are two optimal methods of computing the VDOS with the python interface depending on which trajectory are you using. 
The two function can be found in the `common.py` file

- with the lammps `Traj`: `analyze_vdos_lammps(traj,nstep=None,start=0,resetAccess=True,print=print)`.
   - `traj` is the lammps trajectory file;
   - `nstep` is the number of step to use in the computation of the VDOS; if `None` all steps are included;
   - `start` is the starting step to use as starting point of the trajectory;
   - `resetAccess` if true resets  `traj.setAccessWindowSize(traj.getNtimesteps())` and `traj.setAccessStart(0)`

- with numpy `Trajectory` interface: `analyze_vdos_numpy(pos,vel,types,box,nstep=ltot,start=start)` .
   - `pos`   : positions matrix (N_timesteps, N_atoms, 3)
   - `vel`   : velocities matrix (N_timesteps, N_atoms, 3)
   - `types` : types vector (N_atoms)
   - `box `  : box matrix  
   - `matrix_format` :bool
          matrix format for the box array
   - `wrap `: bool 
    `     `   wrap coordinate
   - `nstep` : int
    `     `   nstep to use for the computation of vdos if None it uses all
   - `start` : int
   `     `   initial step
  The first 6 input are the same of the numpy `Trajectory` interface

Both functions returns a numpy array with dimensions: `(ntypes, freq, Cartesian_coordinate)`.
The units are the same of the `command line version`, thus the matrix must be multiplied by `dt/nstep`.


# Abstract multithreaded calculation C++ class implementation

The `CalcolaMultiThread` class tries to abstract away common features of most of the calculations over the trajectory, allowing to easily write parallel calculations that take advantage of the modern multicore architectures of every computer.

# Abstract block averages calculation for transparent use of MPI

All the logic about MPI and block averages is self contained in `mediablocchi.h`, a small file of only ~200 lines, that take care of doing block averages of both the MPI and the non MPI case. Calculations that needs (or take advantage of) block averages and want to use this tool must implement few functions, as found in all the classes, and must use the Traiettoria class to read the trajectory. In this way the user can concentrate on implementing the actual calculation without thinking of the boilerplate code needed to compute block averages and MPI-parallelize it. The calculation class with the data that must be averaged over the blocks must be a class derived from `OperazioniSuLista`


# Credits

Written by Riccardo Bertossa during his lifetime at SISSA
