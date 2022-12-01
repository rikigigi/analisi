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

    def __init__(self, trajectory):
        self.traj = trajectory

    def compute_steinhardt(self, ranges=None, nthreads=4, skip=10,
                           neigh=None, histogram=True, ls=[6, 4], averaged=False, n_segments=1, nbins=100, tskip=1):
        """
        ranges : list of (min_r, max_r) radial distances where the code will compute the steinhardt parameters
        for each possible pair of types (total n_types**2 elements)
        neigh : list of (n_max, r_max**2, 0.0 ) neighbor list specifications for each possible pair of types. Set to use the SANN algorithm.
        """
        if ranges is None:
            ranges = [ (0.0,2.0) for i in range(len(self.traj.get_atomic_types())) ]
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



