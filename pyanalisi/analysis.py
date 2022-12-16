import numpy as np
import pyanalisi as pa

from .plotters import *

class Analysis:

    @staticmethod
    def pyanalisi_wrapper(Class, traj, *args):
        if isinstance(traj, pa.Trajectory):
            return getattr(pa, Class)(traj, *args)
        elif isinstance(traj, pa.Traj):
            return getattr(pa, Class + '_lammps')(traj, *args)
        else:
            raise RuntimeError(f"Wrapper for trajectory class '{traj.__name__}' not implemented")

    def __init__(self, trajectory,nthreads=4):
        self.traj = trajectory
        self.nthreads=nthreads

    def compute_steinhardt(self, ranges=None, nthreads=None, skip=10,
                           neigh=None, histogram=True, ls=[6, 4], averaged=False, n_segments=1, nbins=100, tskip=1):
        """
        ranges : list of (min_r, max_r) radial distances where the code will compute the steinhardt parameters
        for each possible pair of types (total n_types**2 elements)
        neigh : list of (n_max, r_max**2, 0.0 ) neighbor list specifications for each possible pair of types. Set to use the SANN algorithm.
        """
        if nthreads is None:
            nthreads = self.nthreads
        if ranges is None:
            ranges = [ (0.0,2.0) for i in range(len(self.traj.get_atomic_types())**2) ]
        atraj = self.traj.get_analisi_traj(tskip=tskip,wrapped=True)
        if neigh is None:
            neig=[]
        stein = None
        l = sorted(ls)
    
        if l[0] < 0:
            raise IndexError(f'minimum l cannot be negative {l}')
        elif l[-1] > 10:
            raise IndexError(
                f'spherical harmonics with l>10 are not supported {l}. Please modify and recompile source code.')
    
        stein = Analysis.pyanalisi_wrapper(f'SteinhardtOrderParameterHistogram_{l[-1]}', atraj,
                                  ranges,
                                  1, nbins,
                                  l,
                                  nthreads, skip, histogram, neigh, averaged
                                  )
    
        if n_segments == 1:
            stein.reset(atraj.get_nloaded_timesteps())
            stein.calculate(0)
            stein_res = np.array(stein)
            return stein_res
        elif n_segments > 1:
            segment_size = max(1, atraj.get_nloaded_timesteps() // n_segments)
            stein.reset(segment_size)
            #print(segment_size)
            res = []
            for i in range(0, segment_size * n_segments, segment_size):
                #print(f'calculating {i}...')
                stein.calculate(i)
                res.append(np.array(stein, copy=True))
            return np.array(res)
        else:
            raise IndexError('n_segments must be >= 1')

    def steinhardt_movie(self, skip=1,averaged=True,n_segments=60,plt_steinhardt_kw=None, compute_steinhardt_kw=None):
        if compute_steinhardt_kw is None:
            compute_steinhardt_kw={}
        compute_steinhardt_kw.setdefault('n_segments',n_segments)
        tstein = self.compute_steinhardt(**compute_steinhardt_kw)
        return SteinPlot.animation(tstein,plt_steinhardt_kw=plt_steinhardt_kw)

    @staticmethod
    def max_l(start, stop, tmax=0):
        if tmax <= 0:
            tmax = (stop - start) // 2
        n_ave = stop - start - tmax
        if start >= stop:
            raise RuntimeError("start index must be less than the end index")
        if n_ave <= 0:
            tmax = stop - start - 1
            n_ave = 1
        return tmax, n_ave


    def compute_msd(self,start=0,stop=-1,tmax=0, tskip_msd=10,center_of_mass_MSD=True, center_of_mass_frame=False):
        if stop <0:
            stop = self.traj.numsteps
        tmax, n_ave = Analysis.max_l(start, stop, tmax)
        # mean square displacement
        with Analysis.pyanalisi_wrapper('MeanSquareDisplacement', self.traj.get_analisi_traj(wrapped=False), tskip_msd, tmax, self.nthreads, center_of_mass_MSD, center_of_mass_frame, False) as msd:
           msd.reset(n_ave)
           msd.calculate(start)
           res = np.array(msd,copy=True)
        return res

    @staticmethod
    def block_average_std(data,n_b=6):
        ndp=data.shape[0]//n_b
        means=[]
        for i in range(n_b):
            means.append(np.mean(data[i*ndp:(i+1)*ndp],axis=0))
        means=np.array(means)
        std=(np.var(means,axis=0)/(n_b-1))**.5
        mean=np.mean(means,axis=0)
        del means
        return mean,std

    def block_average(self, key, nblocks=6):
        arr=self.traj.get_array(key)
        return Analysis.block_average_std(arr,n_b=nblocks)
        




