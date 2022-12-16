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

    def steinhardt_movie(self, skip=1,averaged=True,n_segments=60,plt_steinhardt_kw=None, compute_steinhardt_kw=None, plt_steinhardt_init_add_kw=None):
        if compute_steinhardt_kw is None:
            compute_steinhardt_kw={}
        compute_steinhardt_kw.setdefault('n_segments',n_segments)
        tstein = self.compute_steinhardt(**compute_steinhardt_kw)
        return SteinPlot.animation(tstein,plt_steinhardt_kw=plt_steinhardt_kw,plt_steinhardt_init_add_kw=plt_steinhardt_init_add_kw)

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
    def hist2gofr(gr_N, gr_dr, gr_0, gofr):
        rs_m = np.arange(gr_N) * gr_dr + gr_0
        rs_p = (np.arange(gr_N) + 1) * gr_dr + gr_0
        vols = 4 * np.pi / 3 * (rs_p ** 3 - rs_m ** 3)
        return gofr / vols
 

    def compute_gofr(self, startr, endr, nbin,start=0,stop=None, tmax=1, tskip=10, n_segments=1, nthreads=None,return_histogram=False):
        """
           startr -> first r of the g(r) histogram
           endr -> maximum r of the g(r)
           nbin -> number of bins on the radial g(r) histogram
           tmax -> number of time lags ( if > 1, computes van hove function g(r,t) )
        """
        #TODO: uniform usage of tskip and skip
        traj = self.traj.get_analisi_traj(tskip=1,wrapped=True)
        if stop is None:
            stop=traj.get_nloaded_timesteps()
        if nthreads is None:
            nthreads = self.nthreads
        tmax, n_ave = Analysis.max_l(start, stop, tmax)
        gofr = Analysis.pyanalisi_wrapper('Gofrt', traj, startr, endr, nbin, tmax, nthreads, tskip, False, 1)
        if n_segments == 1:
            gofr.reset(n_ave)
            #print('calculating g(r)...', flush=True)
            gofr.calculate(start)
            if return_histogram:
                return np.array(gofr ,copy=True)
            else:
                return Analysis.hist2gofr(nbin,(endr-startr)/nbin, startr, np.array(gofr, copy=True))
        elif n_segments > 1:
            res = []
            segment_size = max(1, n_ave // n_segments)
            gofr.reset(segment_size)
            #print(segment_size)
            for i in range(0, min(segment_size * n_segments, n_ave), segment_size):
                gofr.calculate(start + i)
                if return_histogram:
                    res.append(np.array(gofr, copy=True))
                else:
                    res.append( Analysis.hist2gofr(nbin,(endr-startr)/nbin, startr, np.array(gofr, copy=True)  ) )
            return res
        else:
            raise IndexError(f'n_segments must be > 0 ({n_segments})')

    #### stuff to find peaks of g(r)

    @staticmethod
    def compute_peaks_gofr(gr, r0, gr_0, gr_dr, peak_fractional_height):
        '''
        find the peak around r0 in the g(r) array and get some info
        gr -> g(r) 
        gr_0 -> starting value of r in g(r) array
        gr_dr -> size of the bin of the histogram
        peak_fractional_height -> threshold to trigger the peak width algorithm, in fractions of peak height
        '''
        def peak_width(gr, NH, NH_thre, gr_r2i, min_spread=0.1):
            '''
            get the width of the first peak in the gr at a given height
            '''
            idx_23 = 0
            while gr[0, NH, idx_23] < NH_thre:
                idx_23 += 1
            idx_spread_low = idx_23
            while gr[0, NH, idx_23] > NH_thre * 0.9:
                if idx_23 == gr.shape[2] - 1 or (idx_23 - idx_spread_low > gr_r2i(min_spread) and gr[0, NH, idx_23] < NH_thre):
                    break
                idx_23 += 1
            idx_spread_hi = idx_23
            return idx_spread_low, idx_spread_hi
        def get_conv_functs(gr_0, gr_dr):
            def gr_r2i(r):
                return int((r - gr_0) / gr_dr)
        
            def gr_i2r(i):
                return gr_0 + i * gr_dr
        
            return np.vectorize(gr_r2i), np.vectorize(gr_i2r)
        gr_r2i, gr_i2r = get_conv_functs(gr_0, gr_dr)
        maxs = np.argmax(gr, axis=2)
        peaks = gr_i2r(maxs)

        # find nearest peak to r0
        NH = np.argmin((peaks - r0) ** 2)
        # get spread at 'peak_fractional_height' of peak height
        NH_peak_idx = maxs[0, NH]
        NH_peak_val = gr[0, NH, NH_peak_idx]
        NH_thre = NH_peak_val * peak_fractional_height
        # get width of peak
        idx_spread_low, idx_spread_hi = peak_width(gr, NH, NH_thre, gr_r2i)
        spread_low = gr_i2r(idx_spread_low)
        spread_hi = gr_i2r(idx_spread_hi)
        w23h = spread_hi - spread_low
        return idx_spread_low, idx_spread_hi, spread_low, spread_hi, w23h, NH_thre, NH_peak_val, NH_peak_idx, peaks, NH


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
        




