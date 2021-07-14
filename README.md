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

## installation and test suite

After compiling the code you will find the executable `analisi` and the shared library `pyanalisi.*.so` (with a  part of the name that depends on your particular python installation) in the build folder. To install the python library, simply copy the pyanalisi folder it in your site-library folder, and then copy inside it the `pyanalisi.*.so library`. This is done by the script `install/install_python.sh` after exporting the variables `SP_DIR`, setted to your python's package directory, `BUILD_DIR`, setted to the build directory, and `SOURCE_DIR`, setted to the main directory of the repository. The python's packages folders can be found with the command `python -m site`
```
export SP_DIR="/path/to/python/package/folder"
export SOURCE_DIR="/path/to/repository/dir"
export BUILD_DIR="/path/to/build/directory"
"<img src="svgs/10c587ecb43d57f3321271e69c3b79c4.svg?invert_in_darkmode" align=middle width=7571.819375999999pt height=3359.8173575999995pt/>\bf ^ix_t<img src="svgs/1da9fcafafb5c088ab6ddc5e21359621.svg?invert_in_darkmode" align=middle width=44.86319144999999pt height=22.831056599999986pt/>i\in\{1,\dots,N_{atoms}\}<img src="svgs/254d34f827bc7a00266c9098033ff247.svg?invert_in_darkmode" align=middle width=154.4806758pt height=22.831056599999986pt/>t<img src="svgs/b3197524a9c8b2f48ef4a513a59b8152.svg?invert_in_darkmode" align=middle width=540.578742pt height=22.831056599999986pt/>j<img src="svgs/233aad70e73be16a015a9513dbe4fd6f.svg?invert_in_darkmode" align=middle width=100.93944465pt height=22.831056599999986pt/>t<img src="svgs/83cd70716a0def23c5ab8c8def525d27.svg?invert_in_darkmode" align=middle width=16.39463429999999pt height=14.15524440000002pt/><img src="svgs/abe6f7eaf772e7b2a192800144f61446.svg?invert_in_darkmode" align=middle width=180.00735884999997pt height=31.75825949999999pt/><img src="svgs/ff43df63604b0f704e319703268f961e.svg?invert_in_darkmode" align=middle width=44.86319144999999pt height=22.831056599999986pt/>N_j<img src="svgs/2d4a50d58cfef22e61830dc02a030c65.svg?invert_in_darkmode" align=middle width=240.14770724999997pt height=22.831056599999986pt/>j<img src="svgs/4fc92a961b0e29e1d9399eda1ec14ee8.svg?invert_in_darkmode" align=middle width=225.59357865pt height=22.831056599999986pt/><img src="svgs/95cddcd0b2c0ed8ee77da853b79cde13.svg?invert_in_darkmode" align=middle width=429.54891045pt height=118.3291758pt/><img src="svgs/b032451ff09baef4a5888bacf04411e8.svg?invert_in_darkmode" align=middle width=2405.12971545pt height=22.831056599999986pt/><img src="svgs/700f5c207bfe4aaa8f55c1d26357e553.svg?invert_in_darkmode" align=middle width=627.5146679999999pt height=118.3291758pt/><img src="svgs/569977133408508004cf9a32e3c99070.svg?invert_in_darkmode" align=middle width=384.36798674999994pt height=22.831056599999986pt/>N_{types}<img src="svgs/96a9516fa39f195ae1bd07c2fb7f1255.svg?invert_in_darkmode" align=middle width=1380.3905659499999pt height=45.84475500000001pt/><img src="svgs/d3f4412da6e66155b5b544663cf49de3.svg?invert_in_darkmode" align=middle width=406.26059265pt height=33.5341875pt/><img src="svgs/7842a825ae29b1147def5c3814c749e7.svg?invert_in_darkmode" align=middle width=3367.1345515499997pt height=631.232877pt/>M<img src="svgs/778d31d086b49c1a51eaf66402afdfc6.svg?invert_in_darkmode" align=middle width=187.41725805pt height=22.831056599999986pt/>N<img src="svgs/2e3724020dfd7378bb41e0ce660c93bf.svg?invert_in_darkmode" align=middle width=8.21920935pt height=14.15524440000002pt/>^m {\bf J}_{t}<img src="svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>m\in\{1\dots M\}<img src="svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>t\in\{1\dots N\}<img src="svgs/5a467be78997744678eefa3053da06c0.svg?invert_in_darkmode" align=middle width=450.94967775pt height=22.831056599999986pt/><img src="svgs/c545038b6efb5b787f79d608651eeb6a.svg?invert_in_darkmode" align=middle width=251.32803135pt height=260.8646931pt/><img src="svgs/62017fec7c48d421dc0a3b94339bc94d.svg?invert_in_darkmode" align=middle width=413.95399485pt height=22.831056599999986pt/>^*<img src="svgs/7d99e7f90db2a5dfc0e9f13487a50f1c.svg?invert_in_darkmode" align=middle width=71.15694299999998pt height=22.831056599999986pt/>L_t<img src="svgs/64589a6bf6e9e0a2709393b46a737d59.svg?invert_in_darkmode" align=middle width=651.1841952pt height=22.831056599999986pt/>C_t<img src="svgs/26dd67c80c7a29aec70c005b66f4c806.svg?invert_in_darkmode" align=middle width=197.83695525pt height=22.831056599999986pt/>N_{ave}<img src="svgs/f203d3fb00c64faf5cc15287442b8bfd.svg?invert_in_darkmode" align=middle width=839.5716548999999pt height=22.831056599999986pt/><img src="svgs/659f76ddc18c6f4c18c13b2af74d9ecf.svg?invert_in_darkmode" align=middle width=528.55261305pt height=27.6567522pt/><img src="svgs/6cec578c9296b0ba259bb2f49b2df01e.svg?invert_in_darkmode" align=middle width=1530.0016644pt height=85.29681270000002pt/>j<img src="svgs/d892d9794be3cd13d416f3ac42a1c812.svg?invert_in_darkmode" align=middle width=77.87809589999999pt height=21.68300969999999pt/>l<img src="svgs/decc45ff72b6d9a733efd6f1ce6f547c.svg?invert_in_darkmode" align=middle width=567.9387648pt height=24.7161288pt/>i<img src="svgs/d892d9794be3cd13d416f3ac42a1c812.svg?invert_in_darkmode" align=middle width=77.87809589999999pt height=21.68300969999999pt/>l+t<img src="svgs/40ab8aaea90b361f9d18e014f04ae07a.svg?invert_in_darkmode" align=middle width=1116.9110864999998pt height=24.65753399999998pt/>l<img src="svgs/1cd59ee5498727f7c0f690d566ba2b30.svg?invert_in_darkmode" align=middle width=534.9732910500001pt height=22.831056599999986pt/>i<img src="svgs/fd92a53167b3c6ae9574071613d555dc.svg?invert_in_darkmode" align=middle width=27.11199479999999pt height=22.831056599999986pt/>j<img src="svgs/5d145b27b89dacc3da288e0ebe70a7cd.svg?invert_in_darkmode" align=middle width=91.83558615pt height=22.831056599999986pt/>t<img src="svgs/3a99e396c7225236dde4144a54df1901.svg?invert_in_darkmode" align=middle width=86.56724834999999pt height=21.68300969999999pt/>l<img src="svgs/83cd70716a0def23c5ab8c8def525d27.svg?invert_in_darkmode" align=middle width=16.39463429999999pt height=14.15524440000002pt/><img src="svgs/b4e6418025806933b8dc043c79bb5104.svg?invert_in_darkmode" align=middle width=312.05818545pt height=27.15900329999998pt/><img src="svgs/22f314221c639d3cf0da37e9d6139fce.svg?invert_in_darkmode" align=middle width=99.23088614999998pt height=22.831056599999986pt/>g(r,t)<img src="svgs/ccffde2bd7540699b7edb11e4ebaeb6c.svg?invert_in_darkmode" align=middle width=102.09238214999998pt height=22.831056599999986pt/>I<img src="svgs/fd92a53167b3c6ae9574071613d555dc.svg?invert_in_darkmode" align=middle width=27.11199479999999pt height=22.831056599999986pt/>J<img src="svgs/12709a6ed08b5ca9ceaa8b86d16c4e23.svg?invert_in_darkmode" align=middle width=87.53105789999998pt height=22.831056599999986pt/><img src="svgs/1777ac22904e1a70c8e6697739b29800.svg?invert_in_darkmode" align=middle width=297.81893130000003pt height=27.6567522pt/><img src="svgs/ff43df63604b0f704e319703268f961e.svg?invert_in_darkmode" align=middle width=44.86319144999999pt height=22.831056599999986pt/>\delta(\cdot,\cdot)<img src="svgs/fcd8cc4d70cc76ac5c05d9953b6bcc5e.svg?invert_in_darkmode" align=middle width=777.76797945pt height=39.45205440000001pt/>t,I,J<img src="svgs/f3cb00a03b225cd75b1c04333da2817c.svg?invert_in_darkmode" align=middle width=40.58716694999999pt height=22.831056599999986pt/>J\geq I<img src="svgs/8f1f056462c34995c995cb763769eb79.svg?invert_in_darkmode" align=middle width=327.101346pt height=22.831056599999986pt/>N_{types}(N_{types}+1)/2 - (J+1)(J+2)/2 + I<img src="svgs/7f134f5699d2338aa540d5ff80f24f63.svg?invert_in_darkmode" align=middle width=630.87888105pt height=24.65753399999998pt/><img src="svgs/cde03ef67b50f6aec46c86b58cfa0edb.svg?invert_in_darkmode" align=middle width=519.2253676500001pt height=29.190975000000005pt/><img src="svgs/7f380fbbecc44a9a64dd759fcb70bec9.svg?invert_in_darkmode" align=middle width=700.2744969pt height=85.29680939999997pt/>t<img src="svgs/b1dd19c7daae74ab9c82a38737b56cb6.svg?invert_in_darkmode" align=middle width=683.3613484499999pt height=22.831056599999986pt/>I,J<img src="svgs/1ed5a3d62c832e552de9c6a97fd94368.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/><img src="svgs/a1d518671ff6c434098a596e8f58630a.svg?invert_in_darkmode" align=middle width=306.69572174999996pt height=225.01244039999997pt/><img src="svgs/87d77752cde360e96223ac312edd4254.svg?invert_in_darkmode" align=middle width=211.18775204999997pt height=45.84475499999998pt/>I,J<img src="svgs/b4e39ff12eabae483b7b2892f8015d36.svg?invert_in_darkmode" align=middle width=1242.3666172499998pt height=124.74886710000001pt/><img src="svgs/b46a0311acaf69834fe5cd97a5856b8c.svg?invert_in_darkmode" align=middle width=511.03006845pt height=146.9602365pt/><img src="svgs/ffa563ab114359f9487c40623a4caada.svg?invert_in_darkmode" align=middle width=50.46060359999999pt height=22.831056599999986pt/>\cos \theta<img src="svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>\sin \varphi<img src="svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>\cos \varphi<img src="svgs/44e623a16ab126c6a1c0a2def6e5e1bf.svg?invert_in_darkmode" align=middle width=304.01778494999996pt height=22.831056599999986pt/><img src="svgs/ffb647b24c297c53c99aa0a37e54b403.svg?invert_in_darkmode" align=middle width=164.54838839999996pt height=131.69084279999998pt/><img src="svgs/d247292d83cafa8d119cd35ec0614746.svg?invert_in_darkmode" align=middle width=27.11199479999999pt height=22.831056599999986pt/>\sin |m|\varphi<img src="svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662248499999905pt height=14.15524440000002pt/>\cos |m|\varphi<img src="svgs/093d90e46009c7d352cb0516fb06ea8f.svg?invert_in_darkmode" align=middle width=501.74309789999995pt height=22.831056599999986pt/><img src="svgs/9471f8e188845996e06585e0b6c0188c.svg?invert_in_darkmode" align=middle width=342.84232125pt height=57.53473439999999pt/><img src="svgs/d247292d83cafa8d119cd35ec0614746.svg?invert_in_darkmode" align=middle width=27.11199479999999pt height=22.831056599999986pt/>P_\ell^m<img src="svgs/fe73f7fe03b9453274bfc928c3b5b8c8.svg?invert_in_darkmode" align=middle width=685.5268645499999pt height=22.831056599999986pt/><img src="svgs/ea67b84cd4a6f8b68d2fe39e68eb2e02.svg?invert_in_darkmode" align=middle width=281.84952675pt height=55.908118200000004pt/><img src="svgs/0aec1319a01afe663368b8eef1b473ab.svg?invert_in_darkmode" align=middle width=213.53399594999996pt height=22.831056599999986pt/>(\ell,\ell)<img src="svgs/0c5aa989c7d63a9734d206c1540b65d2.svg?invert_in_darkmode" align=middle width=67.50029219999999pt height=22.831056599999986pt/>(\ell,\ell+1)<img src="svgs/189d93ad150ea9e6b2155186aa87c0a3.svg?invert_in_darkmode" align=middle width=99.27364589999998pt height=22.831056599999986pt/>(\ell,m)<img src="svgs/a118f9264a997cde37060745a97a29c2.svg?invert_in_darkmode" align=middle width=302.63273369999996pt height=22.831056599999986pt/>\ell<img src="svgs/1c1be01d57aecbdef3c673c55063c42a.svg?invert_in_darkmode" align=middle width=138.2264004pt height=22.831056599999986pt/><img src="svgs/9e77fde9f89717512c649736d6db9369.svg?invert_in_darkmode" align=middle width=217.65235965000002pt height=36.76707539999999pt/><img src="svgs/7dc20fa718ebdd1bb8204b54ddb30c36.svg?invert_in_darkmode" align=middle width=203.33389229999997pt height=45.84475499999998pt/>\ell_{max}<img src="svgs/a768aa9cd6efbfe12e44fe78accbae11.svg?invert_in_darkmode" align=middle width=82.36257314999999pt height=22.831056599999986pt/>(x,y,z)<img src="svgs/80c7a6129a1f2c46cbeaf4fcd9978dbf.svg?invert_in_darkmode" align=middle width=225.98178075pt height=22.831056599999986pt/>Y_{\ell m}(x,y,z) \forall \ell \in [0,\ell_{max}]<img src="svgs/90a8dd61b967e18624439c25775ff3fc.svg?invert_in_darkmode" align=middle width=189.7420602pt height=22.831056599999986pt/>m<img src="svgs/09e82c2a86c7519bfc35384c943b60a7.svg?invert_in_darkmode" align=middle width=186.99032054999998pt height=22.831056599999986pt/><img src="svgs/73410d17ee5ac3553527b63aec14e3fd.svg?invert_in_darkmode" align=middle width=295.07753549999995pt height=31.525041899999984pt/><img src="svgs/1eea7d6700c6bffe8917d1e5daf5113d.svg?invert_in_darkmode" align=middle width=126.67203119999999pt height=22.831056599999986pt/>t<img src="svgs/46934cde3bc8ef0ea9fd875cd2cec979.svg?invert_in_darkmode" align=middle width=106.54262024999997pt height=22.831056599999986pt/>V_I<img src="svgs/50a34f6f0d7de86fdc1014ec37d669d9.svg?invert_in_darkmode" align=middle width=487.67158274999997pt height=22.831056599999986pt/>I<img src="svgs/001e05e2e485069e1be9124f247ad9b5.svg?invert_in_darkmode" align=middle width=65.68251525pt height=22.831056599999986pt/>r_{inner}<img src="svgs/fd92a53167b3c6ae9574071613d555dc.svg?invert_in_darkmode" align=middle width=27.11199479999999pt height=22.831056599999986pt/>r_{outer}<img src="svgs/46ae9b87c4b40c3f335224127b7611a4.svg?invert_in_darkmode" align=middle width=79.28106944999999pt height=22.831056599999986pt/>\rho_j<img src="svgs/26c7280cf01adbd7c1b3525e0404442c.svg?invert_in_darkmode" align=middle width=232.8781224pt height=22.831056599999986pt/>j<img src="svgs/e04236f4a12e5f0d1138cc0029dda413.svg?invert_in_darkmode" align=middle width=1266.3017845499999pt height=22.831056599999986pt/><img src="svgs/dd25ab8ccfaa41417443f76fd7df5d02.svg?invert_in_darkmode" align=middle width=349.56064935pt height=32.51169900000002pt/><img src="svgs/ff43df63604b0f704e319703268f961e.svg?invert_in_darkmode" align=middle width=44.86319144999999pt height=22.831056599999986pt/>\langle \cdot \rangle<img src="svgs/0d8bab37edb5a64ecb206f842cbcffc1.svg?invert_in_darkmode" align=middle width=446.77807679999995pt height=22.831056599999986pt/>N_J<img src="svgs/0658d0846fdd4a92c0eafc3e5f4244f4.svg?invert_in_darkmode" align=middle width=168.4500741pt height=22.831056599999986pt/>J<img src="svgs/2ae9d822b4563e09bfb4c9da84ba1959.svg?invert_in_darkmode" align=middle width=33.58078679999999pt height=22.831056599999986pt/>\langle \cdot \rangle<img src="svgs/2d3b40d68cdb2012e3f121f938f7d431.svg?invert_in_darkmode" align=middle width=441.7981821pt height=85.29681270000002pt/><img src="svgs/f17a3962dae2b14e7b818d9e36ff1e3c.svg?invert_in_darkmode" align=middle width=315.61073444999994pt height=32.256008400000006pt/><img src="svgs/d58337d48610f0bf4202bfdfca7f94fe.svg?invert_in_darkmode" align=middle width=30.182742149999992pt height=39.45205439999997pt/>\alpha$ is the type of atom.
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
