import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Circle
import matplotlib.colors as colors
from matplotlib import cm
import numpy as np
import matplotlib
from IPython.core.display import display, HTML

class SteinPlot:
    @staticmethod
    def plot(stein_res, vmin=0.01, figsize=(6., 6.), show=True, transpose=True, xmax=0.2, ymax=0.6,
                       inverted_type_index=False, axs=None, fig=None, single=None, plt_points=None):
        if len(stein_res.shape) != 5:
            raise RuntimeError('implemented only for 2d histograms!')
        cmap = cm.get_cmap('inferno').copy()
        cmap.set_bad('black')
        nt = stein_res.shape[1]
        mask = np.zeros(stein_res.shape)
        # mask[:,:,:,-1,-1]=1
        masked = np.ma.masked_array(stein_res, mask=mask)
        if axs is None:
            if fig is None:
                fig = plt.figure(figsize=figsize)
            axs = ImageGrid(fig, 111, nrows_ncols=(nt, nt) if single is None else (1, 1), axes_pad=0.3)
        idx = 0
        for itype in range(nt):
            for jtype in range(nt):
                if single is not None:
                    if (itype, jtype) != single:
                        continue
                if inverted_type_index:
                    itype, jtype = jtype, itype
                try:
                    axs[idx].imshow(stein_res[0, itype, jtype].transpose() if transpose else stein_res[0, itype, jtype],
                                    norm=colors.LogNorm(vmin=vmin, vmax=masked[0, itype, jtype].max()),
                                    cmap=cmap,
                                    origin='lower', extent=[0.0, 1.0, 0.0, 1.0],
                                    aspect=xmax / ymax)
                    axs[idx].set_xlim(0, xmax)
                    axs[idx].set_ylim(0, ymax)
                    if plt_points:
                        for x, y, r, c_kw in plt_points:
                            circ = Circle((x, y), radius=r, **c_kw)
                            axs[idx].add_patch(circ)
                except Exception as e:
                    print(e)
                idx += 1
        if show:
            fig.show()
        return fig, axs

    @staticmethod
    def animation(tstein, plt_steinhardt_kw=None):
        if plt_steinhardt_kw is None:
            plt_steinhardt_kw = {}
        n_segments = tstein.shape[0]
        class SteinAni:
            def __init__(self, tstein, plt_steinhardt_kw):
                self.tstein = tstein
                self.plt_steinhardt_kw = plt_steinhardt_kw
                self.fig, self.axs = SteinPlot.plot(self.tstein[0], **self.plt_steinhardt_kw, show=False)
    
            def __call__(self, i):
                self.fig, self.axs = SteinPlot.plot(self.tstein[i], **self.plt_steinhardt_kw, axs=self.axs, fig=self.fig,
                                                    show=False)
                return self.axs
    
        stani = SteinAni(tstein, plt_steinhardt_kw)
        ani = matplotlib.animation.FuncAnimation(
            stani.fig, stani, interval=200, blit=True, save_count=n_segments)
        return HTML(ani.to_jshtml())

