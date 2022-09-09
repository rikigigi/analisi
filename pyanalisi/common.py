import re
import pickle
import math

try:
    import aiida
    HAS_AIIDA=True
except:
    HAS_AIIDA=False
if HAS_AIIDA:
    try:
       from aiida import load_profile
       from aiida.orm import Code, load_node
       from aiida.orm import *
       from aiida.engine import submit
       load_profile()
    except:
       HAS_AIIDA=False

try:
    import aiida_QECpWorkChain
    HAS_QECPWORKCHAIN=True
except:
    HAS_QECPWORKCHAIN=False
        
if HAS_QECPWORKCHAIN:
    from aiida_QECpWorkChain.workflow import *
    from aiida_QECpWorkChain import write_cp_traj

try:
    import numpy as np
except:
    print('WARNING: cannot import numpy')
try:
   import matplotlib as mpl
   mpl.rcParams['agg.path.chunksize'] = 1000000
   import matplotlib.pyplot as plt
   import matplotlib.colors as pltc
   from mpl_toolkits.axes_grid1.inset_locator import inset_axes
except:
   print('WARNING: cannot import matplotlib')


import pyanalisi.pyanalisi as pa

print(pa.info())

from matplotlib import collections  as mc
from IPython.core.display import display, HTML
import matplotlib.animation

FIGURE_PATH = '.'
def set_figure_path(new_path):
    global FIGURE_PATH
    FIGURE_PATH = new_path

DEFAULT_PLT_STEINHARDT_KW={'transpose':True,'xmax':.30,'ymax':.60}
DEFAULT_NEIGH=[(57,3.5**2,0.0),(45,3.5**2,0.0)]
#aiida
def plt_key(traj,key,conv=1.0,title='',ylabel=''):
    if not key in traj.get_arraynames():
        return
    if title=='': title=key
    fig,ax =plt.subplots(figsize=(8,6),dpi=300)
    ax=fig.add_axes([0,0,1,1])
    axins=inset_axes(ax,width='20%',height='20%',loc='upper right')
    ax.set_title(title)
    ax.set_xlabel('time (ps)')
    ax.set_ylabel(ylabel)
    axins.set_title('100 points of "{}"'.format(title))
    t=traj.get_array('times')
    q=traj.get_array(key)
    ax.plot(t,q*conv)
    axins.plot(t,q*conv)
    x1=len(t)//2
    x2=len(t)//2+100
    if x2==x1: x2=len(t)-1
    if x2>=len(t): x2=len(t)-1
    axins.set_xlim(t[x1],t[x2])
    return ax

    

def get_types_id_array(types_array):
    types={}
    typeid=0
    res=[]
    for t in types_array:
        if not t in types:
            types[t]=typeid
            typeid=typeid+1
        res.append(types[t])
    return np.array(res,dtype='int32')

#aiida
def get_analisi_traj_from_aiida(traj):
    pos=traj.get_array('positions')
    if 'velocities' in traj.get_arraynames():
        vel=traj.get_array('velocities')
    else:
        vel=np.zeros(pos.shape)
    cel=traj.get_array('cells').transpose((0,2,1)).copy(order='C') 
    types=get_types_id_array(traj.get_attribute('symbols'))
    params= [pos, vel, types,  cel]
    atraj=pa.Trajectory(*params,pa.BoxFormat.CellVectors, True,True)
    atraj_unw=pa.Trajectory(*params,pa.BoxFormat.CellVectors, False,True)
    return atraj, atraj_unw

#only pos and cell
def get_analisi_traj(pos,types,cell,wrap=False):
    return pa.Trajectory(pos,np.zeros(pos.shape),types,cell,pa.BoxFormat.CellVectors, wrap,True)
try:
    import k3d
    HAS_K3D=True
except:
    HAS_K3D=False
try:
    from ipywidgets import interact, interactive, fixed
    import ipywidgets as widgets
    HAS_IPYWIDGETS=True
except:
    HAS_IPYWIDGETS=False    
    
    
def pyanalisi_wrapper(Class,traj,*args):
    if isinstance(traj, pa.Trajectory):
        return getattr(pa,Class)(traj,*args)
    elif isinstance(traj, pa.Traj):
        return getattr(pa,Class+'_lammps')(traj,*args)
    else:
        raise RuntimeError(f"Wrapper for trajectory class '{traj.__name__}' not implemented")
        
def analisi_to_cell_vectors(b):
    boxes=[]
    for i in range(b.shape[0]):
        bb=b[i]
        boxes.append([[bb[3]*2,bb[6],bb[7]],[0,bb[4]*2,bb[8]],[0,0,bb[5]*2]])
    return np.array(boxes)

def fft_density(res,coeff=None):
    if coeff is None:
        coeff=np.ones(res.shape[0])
    resg=np.einsum('ijkl,i',np.fft.fftn(res,axes=(1,2,3)),coeff)
    resg[0,0,0]=0
    resg=np.fft.fftshift(resg)
    resgn=np.real(resg*np.conj(resg))
    resgn/=resgn.max()
    return resgn

def max_l(start,stop,tmax=0):
    if tmax<=0:
        tmax=(stop-start)//2
    n_ave=stop-start-tmax
    if start>=stop:
        raise RuntimeError("start index must be less than the end index")
    if n_ave<=0:
        tmax=stop-start -1
        n_ave=1
    return tmax,n_ave



def analyze_msd(traj_unw,start,stop,tmax=0,nthreads=4,tskip_msd=10,print=print):
    tmax,n_ave = max_l(start,stop,tmax)
    #mean square displacement
    msd=pyanalisi_wrapper('MeanSquareDisplacement',traj_unw,tskip_msd,tmax,nthreads,True,False,False)
    msd.reset(n_ave)
    print('calculating msd...',flush=True)
    msd.calculate(start)
    return np.array(msd)#,copy=True)

def analyze_gofr(traj,start,stop,startr,endr,nbin,tmax=1,nthreads=4,tskip=10,print=print,n_segments=1):
    tmax,n_ave = max_l(start,stop,tmax)
    gofr=pyanalisi_wrapper('Gofrt',traj,startr,endr,nbin,tmax,nthreads,tskip,False,1)
    if n_segments==1:
        gofr.reset(n_ave)
        print('calculating g(r)...',flush=True)
        gofr.calculate(start)
        return np.array(gofr)#,copy=True)
    elif n_segments>1:
        res = []
        segment_size=max(1,n_ave//n_segments)
        gofr.reset(segment_size)
        print(segment_size)
        for i in range(0,min(segment_size*n_segments,n_ave),segment_size):
            gofr.calculate(start+i)
            res.append(np.array(gofr,copy=True))
        return res
    else:
        raise IndexError(f'n_segments must be > 0 ({n_segments})')

def analyze_vdos_single(traj,nstep=None,print=print):
    '''
      Computes vdos 
      there are no start,stop options because spettrovibrazionale.calcola  do not  support the features
      Parameters
      ----------
        traj : Trajectory instances
        nstep : int
                nstep to use for the computation of vdos if None it uses all
      return:
        vdos: numpy array (ntypes , nstep/2+1,3) 
    '''
    if nstep is None or nstep>=traj.getNtimesteps():
       ntep = traj.getNtimesteps()

    vdos=pyanalisi_wrapper('VibrationSpectrum',traj,False)
    vdos.reset(nstep)
    print('calculating vdos...',flush=True)
    vdos.calculate(0)

    return np.array(vdos)#,copy=True)

def analyze_vdos_lammps(traj,nstep=None,start=0,resetAccess=True,print=print):
    '''
      Computes vdos from a lammps trajectory
      Parameters
      ----------
        traj : Trajectory instances
        nstep : int
                nstep to use for the computation of vdos if None it uses all
        start : int
                initial step
        resetAccess : bool
                at the end reset the access to nstep and first step
      return:
        vdos: numpy array (ntypes , nstep/2+1,3) 
    '''
    if nstep is None or nstep>=traj.getNtimesteps():
       nstep = traj.getNtimesteps()
       if start!=0:
          start=0
          print('start setted to zero because all the timesteps are have been requested')
    
    #set the trajectory to what we want
    traj.setAccessWindowSize(nstep)
    traj.setAccessStart(start)
    vdos = analyze_vdos_single(traj,nstep=nstep,print=print)

    #reset Access 
    if (resetAccess):
       traj.setAccessWindowSize(traj.getNtimesteps())
       traj.setAccessStart(0)
       print('reset Access')
       print('traj.setAccessWindowSize(traj.getNtimesteps())')
       print('traj.setAccessStart(0)')
       
    return vdos #,copy=True)

def analyze_vdos_numpy(pos, vel, types, box,
                                    matrix_format=True, # matrix format for the box array
                                    wrap=False # don't wrap the coordinates
                                   ,nstep=None,start=0,
                                    print=print):
    '''
      Computes vdos from numpy arrays
      Parameters
      ----------
        pos : positions matrix (N_timesteps, N_atoms, 3)
        vel : velocities matrix (N_timesteps, N_atoms, 3)
        types : types vector (N_atoms)
        box : box matrix  
        matrix_format :bool
                matrix format for the box array
        wrap : bool 
                wrap coordinate
        nstep : int
                nstep to use for the computation of vdos if None it uses all
        start : int
                initial step
      return:
        vdos: numpy array (ntypes , nstep/2+1,3) 
    '''
    totnstep = len(pos)
    if nstep is None or nstep>=totnstep:
       nstep = totnstep 
       if start!=0:
          start=0
          print('start setted to zero because all the timesteps are have been requested')
    
    traj = pa.Trajectory(pos[start:start+nstep,:,:], vel[start:nstep+start,:,:], types, box[start:start+nstep,:,:],
                                    matrix_format, # matrix format for the box array
                                    wrap # don't wrap the coordinates
                                    )
    #set the trajectory to what we want
    vdos = analyze_vdos_single(traj,nstep=nstep,print=print)

       
    return vdos #,copy=True)

def analyze_sh(traj,start,stop,startr,endr, nbin,ntypes=2, tmax=0, nthreads=4,tskip=10,print=print,
           sann=[], #list to tell how to build the neighbours list: [ (max_number_of_neigh,rcut**2,not_used) ]
           buffer_size=20
      ):
    tmax,n_ave = max_l(start,stop,tmax)
    sh=pyanalisi_wrapper('SphericalCorrelations',traj
                                     ,[(startr #minimum of radial distance
                                     ,endr)]*ntypes**2 #maximum of radial distance
                                     ,nbin # number of distance bins
                                     ,tmax # maximum length of correlation functions in timesteps
                                     ,nthreads # number of threads
                                     ,tskip # time skip to average over the trajectory
                                     ,buffer_size # buffer for sh
                                     ,True,sann) #flag that at the moment does nothing, sann if it is not empty activates the sann algorithm.
    sh.reset(n_ave) # number of timesteps to average
    print('calculating spherical harmonics correlation functions... (go away and have a coffee)',flush=True)
    sh.calculate(start) #calculate starting at this timestep    
    return np.array(sh)#,copy=True)

def analyze(traj,traj_unw,start,stop,nthreads=4,
            tmax_msd=0,
            tskip_msd=10,
            startr_gofr=0.5,
            endr_gofr=3.8,
            nbin_gofr=100,
            tskip_gofr=10,
            tmax_gofr=1,
            startr_sh=0.7,
            endr_sh=1.4,
            nbin_sh=1,
            tskip_sh=20,
            tmax_sh=0):
    
    msd=analyze_msd(traj_unw,start,stop,
                    nthreads=nthreads,
                    tskip_msd=tskip_msd,
                    tmax=tmax_msd)
    
    #gofr (with time)
    gofr=analyze_gofr(traj,start,stop,startr_gofr,endr_gofr,nbin_gofr,
                      tmax=tmax_gofr,
                      nthreads=nthreads,
                      tskip=tskip_gofr)

    
    #spherical harmonics correlations
    sh = analyze_sh(traj,start,stop,startr_sh,endr_sh, nbin_sh,
                    tmax=tmax_sh,
                    nthreads=nthreads,
                    tskip=tskip_sh)

    print('all done',flush=True)
    return msd, gofr, sh

try:
    import scipy as sp
    import scipy.optimize 
except:
    print('WARNING: cannot import scipy')



def _int_gaussian(a,b):
    return a*(np.pi/b)**0.5
def _int_exp(a,b):
    return a/b

def int_fit2_(p):
    freqs =[p[1]] + [ p[(i+1)*2+1] for i in range(0,(p.shape[0]-4)//2,2)] + [p[-2]]
    magns = [(1-p[-1])*_int_gaussian(*p[0:2])]+[ 
         (1-p[-1])*_int_exp(*p[(i+1)*2:(i+2)*2]) 
        for i in range(0,(p.shape[0]-4)//2,2)] + [
       (1-p[-1])*_int_exp(double_exp_aux(p),p[-2])  ]
    return [freqs,magns]

def _meta_generate_initial(n):
    s=''
    for i in range(n):
        s+=f'{1.0/2**(i+1)},10*{10**-i},'
    s=f'[{s}{10**-n},0]'
    #print(s)
    return s

def _meta_generate_n_exp(n):
    names={1:'single', 2:'double', 3:'triple',4:'quadruple'}
    b_min='[0,0'
    b_max='[np.inf,np.inf'
    s=f'''
def {names.get(n,"_"+str(n))}_exp(x, *a):
    return a[{n*2+1}]+(1.0-a[{n*2+1}])*('''
    m='a[0]+'
    for i in range(1,n):
        b_min+=',0,0'
        b_max+=',np.inf,1.0'
        s+=f'a[{i*2}]*np.exp(-a[{i*2+1}]*x)+'
        m+=f'a[{i*2}]+'
    last_mult=f'(1-({m}0))'
    s+=f'a[0]*np.exp(-(a[1]*x)**2)+{last_mult}*np.exp(-a[{n*2}]*x) )'
    b_min+=',0,0]'
    b_max+=',np.inf,np.inf]'
    s_aux=f'def {names.get(n,"_"+str(n))}_exp_aux(a): return {last_mult}'
    #print(s)
    #print(s_aux)
    exec(s, globals())
    exec(s_aux,globals())
    s=f'''
def fit_{names.get(n,"_"+str(n))}_exp(datax,datay,p0={_meta_generate_initial(n)},bounds=(0,np.inf)):
    #filter nan
    res=p0, None
    idxs=~np.isnan(datay)
    if datay[idxs].size==0:
        return res
    try:
        p0[-1] = np.mean(datay[idxs][-5:])
        if p0[-1]<0.0:
            p0[-1] = 0
        if p0[-1] > 1.0 :
            p0[-1] = 1.0
        res= sp.optimize.curve_fit({names.get(n,"_"+str(n))}_exp,datax[idxs],datay[idxs],p0=p0,bounds=bounds)
        return res
    except:
        return res          
'''
    #print(s)
    exec(s,globals())

for i in range(1,4):
    _meta_generate_n_exp(i)
    
def get_exp_spaced_idxs(t,scale=0.1):
    idxs=[]
    skips=np.array(np.exp((t-t[0])*scale),dtype=int)
    idx=0
    while idx<t.shape[0]:
        idxs.append(idx)
        idx+=skips[idx]
    return idxs

def plt_err(ax,x,v,var,*args,**kwargs):
    if var is not None:
        ax.fill_between(x,v-var**.5,v+var**.5,*args,**kwargs)
    else:
        ax.plot(x,v,*args,**kwargs)

def plot_msd(times,res,cm,title='',res_var=None, fig_ax=None):
    if fig_ax is not None:
       fig, ax = fig_ax
    else:
       fig,ax =plt.subplots(figsize=(10,8),dpi=300)
    ax=fig.add_axes([0,0,1,1])
    for i in range(res.shape[-1]):
        plt_err(ax,times[:res.shape[0]]-times[0],res[:,cm,i],res_var[:,cm,i] if res_var is not None else res_var,label='type='+str(i))
    ax.legend()
    ax.set_title('{}MSD'.format(title))
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('$\AA^2/ps$')
    return fig,ax
    
def norm_fit(p):
    """ put to zero the last exponential and move it to the constant part of the fit,
    assuming that the exponential is approximable with a constant in the relevant scale
    """
    out=np.array(p,copy=True)
    n = (out.shape[0]-2)//2
    f = np.sum(out[0:2*n:2])
    out[:2*n:2]=out[:2*n:2]/f
    out[-1]=1+(out[-1]-1)*f
    out[-2]=0
    return out

def fit_sh(times,res,type1,type2,ibin,lmin=0,lmax=11,
           maxt=-1.0,fitmax=-1.0,fitmin=-1.0,pre_fit=-1.0,
           rescale_array=True,scale_exp_coeff=0.1):
    
    if rescale_array:
        idxs=get_exp_spaced_idxs(times[:res.shape[0]],scale_exp_coeff)
        t = times[idxs]-times[0]
        r = res[idxs]
    else:
        t=times[:res.shape[0]]-times[0]
        r = res
    def get_idx(times,maxt):
        idx_end=times.shape[0]
        if maxt<0.0 or maxt>times[-1]:
            maxt=times[-1]
        else:
            while idx_end>0 and times[idx_end-1]>maxt:
                idx_end=idx_end-1
            if idx_end==0:
                pass
        return maxt,idx_end
    maxt,idx_end=get_idx(t,maxt)
    if fitmax<0:
        fitmax=maxt
    fitmax,idx_fit=get_idx(t,fitmax)
    if fitmin<0:
        fitmin=times[0]
    fitmin,idx_fit_beg=get_idx(t,fitmin)
    if idx_fit_beg >= idx_fit:
        raise IndexError('The provided upper time limit is too low')
    print ('fit from {} ({}) to {} ({})'.format(fitmin,idx_fit_beg,fitmax,idx_fit))
    fits=[]
    for l in range(lmin,lmax):
        i=10-l
        print('l={}: {}'.format(l, r[0,type1,type2,ibin,i]))
        kwargs={}
        if pre_fit>0:
            idx_fit_0_t, idx_fit_0 = get_idx(t,pre_fit)
            print ('pre-fit from {} ({}) to {} ({})'.format(fitmin,idx_fit_beg,idx_fit_0_t,idx_fit_0))
            if idx_fit_0 <=idx_fit_beg or idx_fit_0 >= idx_fit:
                raise IndexError('The provided upper time limit for the pre-fit is in the wrong range')
            p0,cv=fit_double_exp(t[idx_fit_beg:idx_fit_0], r[idx_fit_beg:idx_fit_0,type1, type2, ibin, i]/r[0,type1,type2,ibin,i])
            kwargs['p0'] = p0
        p,cv=fit_double_exp(t[idx_fit_beg:idx_fit],r[idx_fit_beg:idx_fit,type1,type2,ibin,i]/r[0,type1,type2,ibin,i],**kwargs)
        if abs(1.0-math.exp(-p[-2]*t[idx_fit-1]) )<1e-5: #safely assume that the last exponential is approximable by a constant
            p = norm_fit(p)
        fits.append(p)
    return fits


def plot_sh(startr,endr,times,res,type1,type2,ibin,lmin=0,lmax=11,title='',res_var=None,
            maxt=-1.0,fitmax=-1.0,fitmin=-1.0,pre_fit=-1.0,log=True,rescale_array=True,scale_exp_coeff=0.1,
            fig_ax=None):
    t=times[:res.shape[0]]-times[0]
    r = res
    r_var = res_var
    def get_idx(times,maxt):
        idx_end=times.shape[0]
        if maxt<0.0 or maxt>times[-1]:
            maxt=times[-1]
        else:
            while idx_end>0 and times[idx_end-1]>maxt:
                idx_end=idx_end-1
            if idx_end==0:
                pass
        return maxt,idx_end
    maxt,idx_end=get_idx(t,maxt)
    dr=(endr-startr)/r.shape[3]
    if fig_ax is not None:
       fig, ax = fig_ax
    else:
       fig,ax =plt.subplots(figsize=(10,8),dpi=300)
    ax=fig.add_axes([0,0,1,1])
    axins=inset_axes(ax,width='20%',height='20%',loc='upper right')
    fits = fit_sh(t,res,type1,type2,ibin,lmin,lmax,maxt,fitmax,fitmin,pre_fit,
                 rescale_array,scale_exp_coeff)
    for l in range(lmin,lmax):
        i=10-l
        p=fits[l-lmin]
        def plot(ax):
            with np.printoptions(precision=3, suppress=True):
                ax.plot(t,double_exp(t,*p),label='fit l={} p={}'.format(l,p))
            plt_err(ax,t,r[:,type1,type2,ibin,i]/r[0,type1,type2,ibin,i],
                    (r_var[:,type1,type2,ibin,i]/r[0,type1,type2,ibin,i]**2) if r_var is not None else r_var,label='l='+str(l))
        plot(ax)
        plot(axins)
    ax.legend()
    axins.set_xlim(0,maxt/50.0)
    #ax.set_xscale('log')
    if log:
        ax.set_yscale('log')
        axins.set_yscale('log')
    ax.set_title('{}$r \in [{:.3f} ,{:.3f}[$'.format(title,startr+dr*ibin,startr+dr*(ibin+1)))
    ax.set_xlim(0,maxt)
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('$c_{{\ell}}^{{{0}\,{1}}}(t)/c_{{\ell}}^{{{0}\,{1}}}(0)$'.format(type1,type2))
    return fig,ax,axins, fits


def plot_gofr(startr,endr,res,title='',res_var=None,fig_ax=None):
    if fig_ax is not None:
       fig, ax = fig_ax
    else:
       fig,ax =plt.subplots(figsize=(10,8),dpi=300)
    rm=np.arange(res.shape[-1])*(endr-startr)/res.shape[-1]+startr
    rp=(np.arange(res.shape[-1])+1)*(endr-startr)/res.shape[-1]+startr
    vols=4*np.pi/3*(rp**3-rm**3)
    for i in range(0,res.shape[1]//2):
        plt_err(ax,(rm+rp)/2,res[0,i,:]/vols,(res_var[0,i,:]/(vols**4)) if res_var is not None else res_var)
    ax.set_xlabel('$\AA$')
    ax.set_ylabel('number of atoms in the shell / volume of the shell')
    ax.set_title('{}$g(r)$'.format(title))
    return fig, ax

def sph_tangent(r,phi,theta,x,y,z):
    a=math.cos(theta)*math.sin(phi)
    b=math.sin(theta)*math.sin(phi)
    c=math.cos(phi)
    return [a,b,c,r-a*x-b*y-c*z]
def sph_tangents_cover(r,n,x,y,z):
    res=[]
    for i in range(n+1):
        phi=math.pi*(i-n/2.0)/n
        n_theta=max(1,int(n*math.cos(phi)))
        for j in range(n_theta):
            theta=2*math.pi*j/n_theta
            res.append(sph_tangent(r,phi+math.pi/2.0,theta,x,y,z))
    return res

def get_max_in_mask(a,mask):
    b=a.view(np.ma.MaskedArray)
    b.mask=mask
    return np.unravel_index(np.argmax(b, axis=None), b.shape)

def idx_to_coord(idx,shape,l=(1.0,1.0,1.0)):
    return (np.array(idx,dtype=float)/np.array(shape,dtype=float)-np.array((0.5,0.5,0.5),dtype=float))*np.array(l)

def is_inside(a,idx,planes):
    idx_c=idx_to_coord(idx,a.shape)
    for plane in planes:
        #print(idx_c)
        if np.sum(idx_c*np.array(plane[:-1]))+plane[3]<0:
            return 1
    return 0

def is_inside_sphere(a,idx,r2,xyz):
    return np.sum((idx_to_coord(idx,a.shape)-xyz)**2)<r2

def get_max_in_spherical_domain(a,r,z,y,x):
    #planes=sph_tangents_cover(r,n,x,y,z)
    mask=np.zeros(a.shape,dtype=int)
    r2=r**2
    xyz=np.array([z,y,x])
    for iz in range(a.shape[0]):
        for iy in range(a.shape[1]):
            for ix in range(a.shape[2]):
                idx=(iz,iy,ix)
                mask[idx]=not is_inside_sphere(a,idx,r2,xyz)
    idx=get_max_in_mask(a,mask)
    return idx, idx_to_coord(idx,a.shape)
            

def mask_around(mask,idx,size):
    idxs_low_high=[[],[],[]]
    for i in range(3):
        low=idx[i]-size[i]
        high=idx[i]+size[i]
        if low < 0:
            if high <= mask.shape[i]:
                idxs_low_high[i].append((0,high))
                idxs_low_high[i].append((mask.shape[i]+low,mask.shape[i]))
            else:
                idxs_low_high[i].append((0,mask.shape[i]))
        elif high > mask.shape[i]:
            if low > 0:
                idxs_low_high[i].append((low,mask.shape[i]))
                idxs_low_high[i].append((0,mask.shape[i]-high))
            else:
                idxs_low_high[i].append((0,mask.shape[i]))
        else:
            idxs_low_high[i].append((low,high))
    for zlo,zhi in idxs_low_high[0]:
        for ylo,yhi in idxs_low_high[1]:
            for xlo,xhi in idxs_low_high[2]:
                mask[zlo:zhi,ylo:yhi,xlo:xhi]=1
                
def d2_minimag(a,b,l):
    d2=0.0
    for i in range(3):
        if abs(a[i]-b[i]) > l[i]/2:
            d2+=(abs(a[i]-b[i])-l[i]/2)**2
        else:
            d2+=(a[i]-b[i])**2
    return d2
def mask_around_spherical(mask,idx,r):
    r2=r**2
    mask_box=np.zeros_like(mask)
    mask_around(mask_box,idx,(int(r),int(r),int(r)))
    for index, m in np.ndenumerate(mask_box):
        if m:
            #check if it is inside the sphere, then mask it
            if d2_minimag(index,idx,mask.shape) <= r2:
                mask[index]=1
    
    

def get_local_maximum_idx(a,threshold=0.0,minimum_distance=(0.1,0.1,0.1),l=(1.0,1.0,1.0)):
    max_list=[]
    max_list_coord=[]
    mask_size=[0,0,0]
    for i in range(3):
        mask_size[i]=max(int(minimum_distance[i]/l[i]*a.shape[i]/2),1)
    b=a.view(np.ma.MaskedArray)
    mask=np.zeros(a.shape)
    while True:
        max_last=np.unravel_index(np.argmax(b, axis=None), b.shape)
        if b[max_last]<threshold:
            break
        max_list.append(max_last)
        idx_=(max_last[2],max_last[1],max_last[0])
        shape_=(b.shape[2],b.shape[1],b.shape[0])
        max_list_coord.append(idx_to_coord(idx_,shape_,l))
        mask_around(mask,max_last,mask_size)
        b.mask=mask
    return max_list,max_list_coord


def get_local_maximum_idx_s(a,threshold=0.0,minimum_distance=0.1,l=1.0):
    max_list=[]
    max_list_coord=[]
    mask_r=minimum_distance/l*max(a.shape)
    b=a.view(np.ma.MaskedArray)
    mask=np.zeros(a.shape)
    while True:
        max_last=np.unravel_index(np.argmax(b, axis=None), b.shape)
        if b[max_last]<threshold:
            break
        max_list.append(max_last)
        idx_=(max_last[2],max_last[1],max_last[0])
        shape_=(b.shape[2],b.shape[1],b.shape[0])
        max_list_coord.append(idx_to_coord(idx_,shape_,l))
        mask_around_spherical(mask,max_last,mask_r)
        b.mask=mask
    return max_list,max_list_coord

if not 'pickled_files' in globals():
    pickled_files={}
if not 'analysis_results' in globals():
    analysis_results={}

def pickle_or_unpickle_reset():
    global pickled_files
    global analysis_results
    pickled_files={}
    analysis_results={}

def pickle_or_unpickle(pickle_dump_name, analisi=None):
    """
    globals that it modifies:
     - pickled (initialized to False)
       - if pickled is True, it does nothing but loading analisi if not defined
     - analisi (initialized with pickle.load)
    """
    global pickled_files
    global analysis_results
    pickled = pickle_dump_name in pickled_files
    analyzed = pickle_dump_name in analysis_results
    if analyzed:
        if not pickled:
            print ('dumping analisi in pickle file "{}"'.format(pickle_dump_name))
            with open(pickle_dump_name, 'wb') as handle:
                pickle.dump(analysis_results[pickle_dump_name], handle, protocol=pickle.HIGHEST_PROTOCOL)
            pickled_files[pickle_dump_name] = True
        else:
            print('analisi was already dumped to the file, or it was already loaded by this function')
    else:
        print('analisi not present: trying to load it from pickle file "{}"'.format(pickle_dump_name))
        try:
            if (pickled):
                raise RuntimeError("logic error")
            with open(pickle_dump_name, 'rb') as handle:
                analysis_results[pickle_dump_name] = pickle.load(handle)
            pickled_files[pickle_dump_name] = True
        except Exception as e:
            print('Error: {}'.format(e))
            if analisi is not None:
                analysis_results[pickle_dump_name] = analisi
                return pickle_or_unpickle(pickle_dump_name, analisi)
            return None
    return analysis_results[pickle_dump_name]


def atomic_density(atraj,dr=0.1,start=-1):
    if start < 0:
        start = atraj.get_current_timestep()
    box=atraj.get_box_copy().mean(axis=0)
    print(box)
    hist=pyanalisi_wrapper('PositionHistogram',atraj,np.array(2*np.rint(box[3:6]/(2*dr)),dtype=int).tolist(),1,1)
    hist.reset(atraj.get_nloaded_timesteps())
    hist.calculate(start)
    res=np.array(hist)
    return res, box


def plot_simulation_box(box,**kwargs):
    kwargs.setdefault('shader','simple')
    a=[2*box[3],0,0]
    if len(box)==9:
       b=[box[6],2*box[4],0]
       c=[box[7],box[8],2*box[5]]
    else:
       b=[0,2*box[4],0]
       c=[0,0,2*box[5]]
    origin=np.array([box[0],box[1],box[2]])
    a=np.array(a)
    b=np.array(b)
    c=np.array(c)
    f1=[origin,origin+a,origin+a+b,origin+b,origin]
    f2=[origin+c,origin+a+c,origin+a+b+c,origin+b+c,origin+c]
    c1=[origin,origin+c]
    c2=[origin+a,origin+a+c]
    c3=[origin+a+b,origin+a+b+c]
    c4=[origin+b,origin+b+c]
    line = k3d.line(f1,**kwargs)
    for l in f2,c1,c2,c3,c4:
        line += k3d.line(l,**kwargs)
    return line

def rotation(a,b,c,order='zxy'):
   mats={}
   mats['x']=np.array([[1, 0,0],
               [0, math.cos(a), -math.sin(a)],
               [0,math.sin(a),math.cos(a)]])
   mats['y']=np.array([[math.cos(b),0,math.sin(b)],
              [0,1,0],
              [-math.sin(b),0,math.cos(b)]])
   mats['z']=np.array([[math.cos(c),-math.sin(c),0],
               [math.sin(c),math.cos(c),0],
               [0,0,1]])
   
   print(order[0])
   xyz=np.copy(mats[order[0]])
   for ax in order[1:]:
      print (ax)
      xyz=mats[ax].dot(xyz)
   r=np.eye(4)
   r[:3,:3]=xyz
   return r

def density_field2(res,box,box_kw={},plot=None,ns=[[1,1,1]]):
    bounds=[ 
                                                        box[0],box[0]+2*box[3],
                                                        box[1],box[1]+2*box[4],
                                                        box[2],box[2]+2*box[5]
                                                                           ]
    print(bounds)
    r0=np.array((box[0]+box[3],box[1]+box[4],box[2]+box[5]))
    if plot is None:
       plot = k3d.plot()
    objs=[]
    if len(res.shape) >3:
        for i in range(res.shape[0]):
            objs.append(k3d.volume(np.array(res[i],dtype='float16'),bounds=bounds))
            plot+=objs[-1]
    else:
        objs.append(k3d.volume(np.array(res,dtype='float16'),bounds=bounds))
        plot+=objs[-1]
        
    #points = k3d.points(g0_sites,point_size=0.02)
    #points1 = k3d.points(g1_sites,point_size=0.005)
    #plot += points
    #plot += points1
    plot += plot_simulation_box(box,**box_kw)
    plot.display()
    global _w_x
    global _x
    global _n_idx
    _w_x=0.05
    _x  =0.14
    _n_idx = 0
    global _rot_a
    global _rot_b
    global _rot_c
    global _rot_Q
    _rot_a = 0.0 
    _rot_b = 0.0
    _rot_c = 0.0
    _rot_Q = np.eye(4)

    nsn=np.array(ns)
    nsn=np.einsum('ij,i->ij',nsn,((nsn**2).sum(axis=1))**-.5) 
    nlabel=widgets.Label(value=f'n={nsn[_n_idx]}')
    display(nlabel)
    def rotate_plane(Q,ps):
       prs=[]
       for p in ps:
          p=np.array(p)
          prs.append(Q.dot(p).tolist())
       return prs
    def get_planes():
        n=nsn[_n_idx]
        c1=-(r0).dot(n) + _w_x +_x
        c2=-(r0).dot(n) - _w_x +_x
        return [
                  [ n[0], n[1], n[2], c1],
                  [-n[0],-n[1],-n[2],-c2]
               ]
    def global_clipping_planes():
            return rotate_plane(_rot_Q, get_planes())
    @interact(idx=widgets.IntSlider(value=0,min=0,max=nsn.shape[0]-1,step=1,description='select index of plane'))
    def g(idx):
        global _n_idx
        _n_idx=idx
        plot.clipping_planes=global_clipping_planes()
        nlabel.value=f'n={nsn[_n_idx]}'
    @interact(x=widgets.FloatSlider(value=box[0]+box[3],min=box[0],max=box[0]+2*box[3],step=box[3]/100))
    def g(x):
        global _x
        _x=x
        plot.clipping_planes=global_clipping_planes()
    @interact(w_x=widgets.FloatSlider(value=box[3]/5,min=0.0,max=2*box[3],step=box[3]/100))
    def g(w_x):
        global _w_x
        _w_x=w_x
        plot.clipping_planes=global_clipping_planes()
    @interact(a=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(a):
        global _rot_a
        global _rot_Q
        _rot_a=a
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()
    @interact(b=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(b):
        global _rot_b
        global _rot_Q
        _rot_b=b
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()
    @interact(c=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(c):
        global _rot_c
        global _rot_Q
        _rot_c=c
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()

    return plot

def density_field(res,box,box_kw={},plot=None):
    bounds=[ 
                                                        box[0],box[0]+2*box[3],
                                                        box[1],box[1]+2*box[4],
                                                        box[2],box[2]+2*box[5]
                                                                           ]
    print(bounds)
    if plot is None:
       plot = k3d.plot()
    objs=[]
    if len(res.shape) >3:
        for i in range(res.shape[0]):
            objs.append(k3d.volume(np.array(res[i],dtype='float16'),bounds=bounds))
            plot+=objs[-1]
    else:
        objs.append(k3d.volume(np.array(res,dtype='float16'),bounds=bounds))
        plot+=objs[-1]
        
    #points = k3d.points(g0_sites,point_size=0.02)
    #points1 = k3d.points(g1_sites,point_size=0.005)
    #plot += points
    #plot += points1
    plot += plot_simulation_box(box,**box_kw)
    plot.display()
    global _w_x
    global _x
    global _w_y
    global _y
    global _w_z
    global _z 
    global _sphere
    _w_x=0.05
    _x  =0.14
    _w_y=0.05
    _y  =0.14
    _w_z=0.05
    _z  =0.14
    _sphere=False
    global _rot_a
    global _rot_b
    global _rot_c
    global _rot_Q
    _rot_a = 0.0 
    _rot_b = 0.0
    _rot_c = 0.0
    _rot_Q = np.eye(4)
 
    def rotate_plane(Q,ps):
       prs=[]
       for p in ps:
          p=np.array(p)
          prs.append(Q.dot(p).tolist())
       return prs
    def global_clipping_planes():
        if _sphere:
            return sph_tangents_cover(_w_x,7,_x,_y,_z)
        else:
            return rotate_plane(_rot_Q,[[1,0,0,_w_x-_x],[-1,0,0,_w_x+_x],
                [0,1,0,_w_y-_y],[0,-1,0,_w_y+_y],
                [0,0,1,_w_z-_z],[0,0,-1,_w_z+_z]
               ])
    button=widgets.Button(description='Center on 1')
    @button.on_click
    def button1_on_click(b):
        global _x
        global _y
        global _z
        print('finding maximum in spherical domain with r={}, center=({},{},{})'.format(_w_x,_x,_y,_z))
        idx, (_z,_y,_x)=get_max_in_spherical_domain(ress[0],_w_x,_z,_y,_x)
        #_x,_y,_z=-_x,-_y,-_z
        print('center=({},{},{}), idx={} shape={}'.format(_x,_y,_z,idx,ress[0].shape))
        print('value={}'.format(ress[0][idx]))
        plot.clipping_planes=global_clipping_planes()
    display(button)#,button2,button3)
    @interact(spherical=widgets.Checkbox(value=False))
    def g(spherical):
        global _sphere
        _sphere=spherical
        plot.clipping_planes=global_clipping_planes()
    @interact(a=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(a):
        global _rot_a
        global _rot_Q
        _rot_a=a
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()
    @interact(b=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(b):
        global _rot_b
        global _rot_Q
        _rot_b=b
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()
    @interact(c=widgets.FloatSlider(value=0,min=0,max=2*np.pi,step=0.01))
    def g(c):
        global _rot_c
        global _rot_Q
        _rot_c=c
        _rot_Q=rotation(_rot_a,_rot_b,_rot_c)
        plot.clipping_planes=global_clipping_planes()

    @interact(x=widgets.FloatSlider(value=box[0]+box[3],min=box[0],max=box[0]+2*box[3],step=box[3]/100))
    def g(x):
        global _x
        _x=x
        plot.clipping_planes=global_clipping_planes()
    @interact(w_x=widgets.FloatSlider(value=box[3]/5,min=0.0,max=2*box[3],step=box[3]/100))
    def g(w_x):
        global _w_x
        _w_x=w_x
        plot.clipping_planes=global_clipping_planes()
    @interact(y=widgets.FloatSlider(value=box[1]+box[4],min=box[1],max=box[1]+2*box[4],step=box[4]/100))
    def g(y):
        global _y
        _y=y
        plot.clipping_planes=global_clipping_planes()
    @interact(w_y=widgets.FloatSlider(value=2*box[4],min=0,max=2*box[4],step=box[4]/100))
    def g(w_y):
        global _w_y
        _w_y=w_y
        plot.clipping_planes=global_clipping_planes()
    @interact(z=widgets.FloatSlider(value=box[2]+box[5],min=box[2],max=box[2]+2*box[5],step=box[5]/100))
    def g(z):
        global _z
        _z=z
        plot.clipping_planes=global_clipping_planes()
    @interact(w_z=widgets.FloatSlider(value=2*box[5],min=0.0,max=2*box[5],step=box[5]/100))
    def g(w_z):
        global _w_z
        _w_z=w_z
        plot.clipping_planes=global_clipping_planes()
    return plot
        
def force_ratio_histogram(wf,print=print,ax=[],create_fig=lambda : plt.subplots(dpi=300)):
    pwcalcjobs=[]
    for c in [x for x in wf.called if str(x.process_class) == "<class 'aiida_quantumespresso.calculations.pw.PwCalculation'>"]:
        pwcalcjobs.append(c)
        print(c.pk)
    res,axs,figs,n_ax=analyze_forces_ratio(pwcalcjobs,minpk=wf.pk,ax_=ax,create_fig = create_fig)
    return res,axs,figs,n_ax

def plot_force_ratio(res,fig=None,ax=None,hheight=10):
    def update_x_low_high(em,dt,d,x,x_size=0.0):
        e=d.setdefault(dt,{}).setdefault(em,{})
        x_lo=e.setdefault('x_lo',x)
        if x-x_size<x_lo: e['x_lo'] = x-x_size
        x_hi=e.setdefault('x_hi',x)
        if x+x_size>x_hi: e['x_hi'] = x+x_size


    #c_k=list(mcolors.TABLEAU_COLORS.keys())
    #idx=0
    if ax==None:
        fig,ax=plt.subplots(dpi=300)
    x_low_high={}
    for e in dict_keys(res,level=2):
        for dt in dict_keys(res,level=1):
            y=[]
            x=[]
            y_err=[]
            for em in res.keys():
                try:
                    update_x_low_high(em,dt,x_low_high,res[em][dt][e]['fratios_mean'],res[em][dt][e]['fratios_std'])
                    x.append(em)
                    y.append(res[em][dt][e]['fratios_mean'])
                    y_err.append(res[em][dt][e]['fratios_std'])
                except KeyError:
                    pass
            if len(x)>0:
                ax.errorbar(x,y,yerr=y_err,label=f'{e} dt={dt:.2f}',alpha=0.7)
    axins=[]
    emass_min=list(res.keys())[0]
    emass_max=list(res.keys())[0]
    for dt in dict_keys(res,level=1):
        for em in res.keys():
            try:
                axi=inset_axes(
                    ax,
                    bbox_transform=ax.transData,
                    bbox_to_anchor=(
                                    em,
                                    x_low_high[dt][em]['x_lo'],
                                    hheight,
                                    x_low_high[dt][em]['x_hi']-x_low_high[dt][em]['x_lo']
                                    ),
                    width="100%", height="100%",borderpad=0
                    )
                for e in dict_keys(res,level=2):
                    x=res[em][dt][e]['hist'][0]
                    y=res[em][dt][e]['hist'][1]
                    y=(y[1:]+y[:-1])/2.0
                    axi.plot(x,y,alpha=0.5)
                    axi.set_ylim([x_low_high[dt][em]['x_lo'],x_low_high[dt][em]['x_hi']])
                    axi.set_axis_off()
                axi.patch.set_alpha(0.5)
                axins.append(axi)
                if em<emass_min: emass_min=em
                if em>emass_max: emass_max=em
            except KeyError:
                pass
    ax.legend(loc=3)
    ax.grid()
    ax.set_xlim(emass_min-hheight,emass_max+hheight)
    ax.set_xlabel('electron fictitious mass')
    ax.set_ylabel('CP/PW force ratio')
    return fig,ax,axins

def length_traj(traj):
    t=traj.get_array('times')
    return t[-1]-t[0]
def get_cp_with_traj(wf,min_t=0.0):
    l=[]
    for c in [x for x in wf.called if str(x.process_class) == "<class 'aiida_quantumespresso.calculations.cp.CpCalculation'>"]:
        if 'output_trajectory' in c.outputs:
            t = length_traj(c.outputs.output_trajectory)
            if t >= min_t:
               l.append(c)
    return sorted(l,key=lambda x: x.pk)

def print_cp_with_traj(wf,print=print,min_t=0.0):
    l = get_cp_with_traj(wf,min_t=min_t)
    for c in l:
        print ("===========")
        print(c.pk,'traj pk={}: {:.3}ps'.format(c.outputs.output_trajectory.pk,length_traj(c.outputs.output_trajectory)) if 'output_trajectory' in c.outputs else 'No output_trajectory')
        print(c.inputs.parameters.get_dict()['IONS'])
        print(c.inputs.parameters.get_dict()['CELL'] if 'CELL' in c.inputs.parameters.get_dict() else '')
    return l

#tools for converting an analisi trajectory to a aiida trajectory, or an object that behaves in a similar way

def analisi_cell2box(box_lammps):
    '''returns cell vectors arranged in ROWS.
       Analisi code always rotate the system to have a triangular cell matrix.'''
    cells=np.zeros((box_lammps.shape[0],3,3))
    for i in range(box_lammps.shape[0]):
        lc=box_lammps[i]
        if lc.shape[0] == 9: #note that this matrix is transposed: the cell vectors are arranged in ROWS
            #                    a_x 0   0
            #                    b_x b_y 0
            #                    c_x c_y c_z
            cells[i]=np.array([ [lc[3]*2, 0,       0      ],
                                [lc[6],   lc[4]*2, 0      ],
                                [lc[7],   lc[8],   lc[5]*2]
                              ])
        elif lc.shape[0] == 6: #for orthorombic cells there is no difference between row arranged and column arranged cell vectors
            cells[i]=np.array([ [lc[3]*2,0,0],
                                 [0,lc[4]*2,0],
                                 [0,0,lc[5]*2]])
        else:
            raise IndexError('wrong shape of cell array')
    return cells

class FakeAiidaT:
    def __init__(self,t,symbols_={},dt=1.0,pk=None):
        self.data={}
        key_l=['positions','cells']
        key_opt=['velocities','energy_constant_motion','ionic_temperature','pressure','electronic_kinetic_energy']
        if isinstance(t,(pa.Traj,pa.Trajectory)):
            self.t=t
            self.symbols=[ symbols_[x] if x in symbols_ else str(x) for x in t.get_lammps_type().tolist()]
            self.data['positions']=t.get_positions_copy()
            self.data['velocities']=t.get_velocities_copy()
            self.box_lammps=t.get_box_copy()
            self.numsteps=self.data['positions'].shape[0]
            self.data['cells']=analisi_cell2box(self.box_lammps)
        elif isinstance(t,dict):
            self.symbols=t['symbols']
            for key in key_l:
                self.data[key]=np.array(t[key])
            for key in key_opt:
                if key in t:
                    self.data[key]=np.array(t[key])
        elif isinstance(t,list):
            self.symbols=t[0]['symbols']
            for key in key_l:
                self.data[key] = np.concatenate( [t[i][key] for i in range(len(t))], axis=0)
            for key in key_opt:
                if key in t[0]:
                    self.data[key]=np.concatenate( [t[i][key] for i in range(len(t))], axis=0)
        else:
            raise RuntimeError(f'first argument cannot be {str(t)}')
        self.data['steps']=np.arange(0,self.data['positions'].shape[0])
        self.data['times']=self.data['steps']*dt
        self.dt=dt
        self.numsites=self.data['positions'].shape[1]
        self.pk=pk
        self.numsteps=self.data['positions'].shape[0]
        
    def get_array(self,name):
        if name in self.data:
            return self.data[name]
        else:
            raise KeyError(f'cannot find key {name}')
    def get_attribute(self,k):
        if k=='symbols':
            return self.symbols
        else:
            raise KeyError(f'attribute {k} not present')
    def get_arraynames(self):
        return self.data.keys()
def analisi2aiida_traj(t,symbols):
    ft = FakeAiidaT(t,symbols)
    res=aiida.orm.nodes.data.array.trajectory.TrajectoryData()
    res.set_attribute('symbols',ft.symbols)
    for name in ['steps','positions','cells','times','velocities']:
        res.set_array(name,ft.get_array(name))
    return res

def read_lammps_bin(f,symbols={},wrap=False,nsteps=0,start=0):
    t=pa.Traj(f)
    t.setWrapPbc(wrap)
    t.setAccessWindowSize(t.getNtimesteps() if nsteps <= 0 else nsteps)
    t.setAccessStart(start)
    return FakeAiidaT(t,symbols_=symbols)

def get_type_mask(at):
    return get_type_mask_from_s(at.symbols)

def get_type_mask_from_s(symbols):
    typesa=set(symbols)
    masks={}
    types=[]
    types_array=np.zeros(len(symbols),dtype=int)
    itype=0
    for t in typesa:
        masks[t]=np.array(symbols)==t
        types.append(t)
        types_array[masks[t]]=itype
        itype+=1
    return types,types_array,masks

def gen_poscar(at,every=1, start=0,desc=''):
    poscars=[]
    types,types_array,masks=get_type_mask(at)
    pos=at.get_array("positions")
    vel=at.get_array("velocities")
    cel=at.get_array("cells")
    nt={}
    for it in types:
        nt[it]=np.sum(masks[it])
    for i in range(start,at.numsteps,every):
        c=cel[i]
        p=f'{desc}\n1.0\n'
        for idim in range(3):
            p+=f'{c[idim,0]} {c[idim,1]} {c[idim,2]}\n'
        for it in types:
            p+=it+' '
        p+='\n'
        for it in types:
            p+=f'{nt[it]} '
        p+='\nCartesian\n'
        for it in types:
            ps=pos[i,masks[it]]
            for iatom in range(nt[it]):
                p+= f'{ps[iatom,0]} {ps[iatom,1]} {ps[iatom,2]} \n'
        p+='Cartesian\n'
        for it in types:
            v=vel[i,masks[it]]
            for iatom in range(nt[it]):
                p+= f'{v[iatom,0]} {v[iatom,1]} {v[iatom,2]} \n'
        
        poscars.append(p)
    return poscars

def wrap_dataset(tt):
    '''This rotates the cell to get a triangular cell matrix, then wraps atomic positions around the cell center in a cubic region of space. Not sure about stress transformation: to check'''
    nat=tt['positions'].shape[1]//3
    nts=tt['positions'].shape[0]
    wrapped=pa.Trajectory(tt['positions'].reshape(nts,nat,3),
                  tt['forces'].reshape(nts,nat,3),
                  np.zeros(nat,dtype='i'),
                  np.array(tt['cells'].reshape(nts,3,3).transpose((0,2,1)),order='C'),
                  pa.BoxFormat.CellVectors,
                  True,
                  True
                  )
    Q=wrapped.get_rotation_matrix()
    rotated_forces=np.einsum('taj,tij -> tai',tt['forces'].reshape(nts,nat,3),Q).reshape((nts,nat*3))
    rotated_stress=np.einsum('tlm,tli,tmj -> tij',tt['stress'].reshape(nts,3,3),Q,Q).reshape((nts,9))
    rotated_positions=np.einsum('taj,tij -> tai',tt['positions'].reshape(nts,nat,3),Q).reshape((nts,nat*3))
    return {'cells':analisi_cell2box(wrapped.get_box_copy()).reshape((nts,9)),
      'cells_unrotated':tt['cells'],
      'Q':Q.reshape((nts,9)),
      'positions':wrapped.get_positions_copy().reshape((nts,nat*3)),
      'forces':wrapped.get_velocities_copy().reshape((nts,nat*3)),
      'energy':tt['energy'],
      'stress':rotated_stress
     }

def show_traj(tr,wrap=True,fast=1.0):
    if wrap:
        atraj,_=get_analisi_traj_from_aiida(tr)
    else:
        atraj=None
    plot=k3d.plot()
    plot=show_atraj(tr,atraj,wrap=wrap,plot=plot,fast=fast)
    plot.display()
    return plot

def show_atraj(tr,atraj,wrap=True,plot=None,fast=1.0):
    atomic_species=list(set(tr.symbols))
    masks=[]
    for sp in atomic_species:
        masks.append(np.array(tr.symbols)==sp)
    if wrap:
        pos_m=atraj.get_positions_copy()
    else:
        pos_m=tr.get_array('positions')
    t=tr.get_array('steps')
    #pos=t.get_positions_copy()
    for sp,mask in zip(atomic_species,masks):
        plot += k3d.points(positions={str(t/30.0/fast):pos_m[t,mask,:] for t in range(pos_m.shape[0])},size=1.0,name=sp)
        print (sp)
    return plot

import matplotlib.colors as colors
from matplotlib import cm


def compute_steinhardt(aiida_traj,ranges=[(2.5,3.5),(0.8,1.2),(2.2,3.0),(1.5,1.8)],nthreads=4,skip=10,neigh=[], histogram=True,ls=[6,4],averaged=False,n_segments=1,nbins=100):
    atraj=aiida_traj
    if isinstance(atraj, FakeAiidaT) and hasattr(atraj,'t'):
        atraj=atraj.t
    elif not isinstance(atraj,pa.Trajectory):
        atraj,atraj_unw=get_analisi_traj_from_aiida(aiida_traj)
    stein = None
    l=sorted(ls)

    if l[0]<0:
       raise IndexError(f'minimum l cannot be negative {l}')
    elif l[-1]>10:
       raise IndexError(f'spherical harmonics with l>10 are not supported {l}. Please modify and recompile source code.')  

    stein=pyanalisi_wrapper(f'SteinhardtOrderParameterHistogram_{l[-1]}',atraj,
                                          ranges,
                                          1,nbins,
                                          l,
                                          nthreads,skip,histogram,neigh,averaged
                                         )

    if n_segments==1:
        stein.reset(atraj.getNtimesteps())
        stein.calculate(0)
        stein_res=np.array(stein)
        return stein_res
    elif n_segments>1:
        segment_size=max(1,atraj.getNtimesteps()//n_segments)
        stein.reset(segment_size)
        print(segment_size)
        res=[]
        for i in range(0,segment_size*n_segments,segment_size):
            print(f'calculating {i}...')
            stein.calculate(i)
            res.append(np.array(stein,copy=True))
        return np.array(res)
    else:
        raise IndexError('n_segments must be >= 1')



from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Circle


def plt_steinhardt(stein_res,vmin=0.01,figsize=(6.,6.),show=True,transpose=True,xmax=0.2,ymax=0.6,inverted_type_index=False,axs=None,fig=None,single=None,plt_points=None):
    if len(stein_res.shape) != 5:
        raise RuntimeError('implemented only for 2d histograms!')
    cmap=cm.get_cmap('inferno').copy()
    cmap.set_bad('black')
    nt=stein_res.shape[1]
    mask=np.zeros(stein_res.shape)
    #mask[:,:,:,-1,-1]=1
    masked=np.ma.masked_array(stein_res,mask=mask)
    if axs is None:
        if fig is None:
            fig = plt.figure(figsize=figsize)
        axs = ImageGrid(fig,111,nrows_ncols=(nt,nt) if single is None else (1,1),axes_pad=0.3)
    idx=0
    for itype in range(nt):
        for jtype in range(nt):
            if single is not None:
                if (itype,jtype) != single:
                    continue
            if inverted_type_index:
                itype, jtype = jtype, itype
            try:
                axs[idx].imshow(stein_res[0,itype,jtype].transpose() if transpose else stein_res[0,itype,jtype],
                                norm=colors.LogNorm(vmin=vmin, vmax=masked[0,itype,jtype].max()),
                                cmap=cmap,
                                origin='lower',extent=[0.0,1.0,0.0,1.0],
                                aspect=xmax/ymax)
                axs[idx].set_xlim(0,xmax)
                axs[idx].set_ylim(0,ymax)
                if plt_points:
                    for x,y,r,c_kw in plt_points:
                        circ = Circle((x,y),radius=r,**c_kw)
                        axs[idx].add_patch(circ)
            except Exception as e:
                print(e)
            idx += 1
    if show: 
        fig.show()
    return fig,axs

def steinhardt_movie(traj,skip=5,neigh=DEFAULT_NEIGH,averaged=True,n_segments=10,plt_steinhardt_kw=DEFAULT_PLT_STEINHARDT_KW,
                     compute_steinhardt_kw={'nthreads':4}):
    tstein=compute_steinhardt(traj,skip=skip,neigh=neigh,averaged=averaged,n_segments=n_segments,**compute_steinhardt_kw)
    class SteinAni:
        def __init__(self,tstein,plt_steinhardt_kw):
            self.tstein=tstein
            self.plt_steinhardt_kw=plt_steinhardt_kw
            self.fig,self.axs=plt_steinhardt(self.tstein[0],**self.plt_steinhardt_kw,show=False)
        def __call__(self,i):
            self.fig,self.axs=plt_steinhardt(self.tstein[i],**self.plt_steinhardt_kw,axs=self.axs,fig=self.fig,show=False)
            return self.axs
        
    stani=SteinAni(tstein,plt_steinhardt_kw)
    ani = matplotlib.animation.FuncAnimation(
        stani.fig, stani, interval=200, blit=True, save_count=n_segments)
    return HTML(ani.to_jshtml())

def gofr_movie(atraj,skip=5,n_segments=10,gofr_kw={},plt_gofr_kw={},callax=lambda x:None):
    gofr_kw.setdefault('gr_kw',{})['tskip']=skip
    if isinstance(atraj, FakeAiidaT) and hasattr(atraj,'t'):
        atraj=atraj.t
    elif not isinstance(atraj,pa.Trajectory):
        atraj,atraj_unw=get_analisi_traj_from_aiida(atraj)
    traj=atraj
    gofr=do_compute_gr_multi(traj,n_segments=n_segments,**gofr_kw)
    class SteinAni:
        def __init__(self,tstein,plt_gofr_kw,callax):
            self.callax=callax
            self.tstein=tstein
            self.plt_gofr_kw=plt_gofr_kw
            self.fig,self.axs,_,_,_,_,_,_,_,_,_,_=do_plots_gr_sh(*self.tstein[0],**self.plt_gofr_kw)
        def __call__(self,i):
            self.axs.clear()
            self.fig,self.axs,_,_,_,_,_,_,_,_,_,_=do_plots_gr_sh(*self.tstein[i],**self.plt_gofr_kw,fig_ax=(self.fig,self.axs))
            self.callax(self.axs)
            return [self.axs]
        
    stani=SteinAni(gofr,plt_gofr_kw,callax)
    ani = matplotlib.animation.FuncAnimation(
        stani.fig, stani, interval=200, blit=True, save_count=n_segments)
    return HTML(ani.to_jshtml())

def elastic_c(kbT, box_matrix):
    "returns elastic constants in GPa"
    is_diagonal=False
    box_ave=box_matrix.mean(axis=0)
    if np.count_nonzero(box_ave-np.diag(np.diagonal(box_ave))) == 0:
        is_diagonal=True
    vol_mean=np.linalg.det(box_ave)
    G=np.einsum('tji,tjl->til',box_matrix,box_matrix)
    box_inv=np.linalg.inv(box_matrix)
    box_ave_inv=box_inv.mean(axis=0)
    strain=0.5*(np.einsum('ji,tjl,lm->tim',box_ave_inv,G,box_ave_inv)-np.eye(3))
    strain_0 = strain.mean(axis=0)
    nt=strain.shape[0]
    #pick=[0,1,2,4,5,8] # xx, xy, xz, yy, yz, zz
    pick=[0,4,8,1,2,5] # xx, yy, zz, xy, xz, yz
    if is_diagonal:
        pick=pick[:3]
    CS_minus1=np.einsum('ti,tk->ik',(strain-strain_0).reshape((nt,9))[:,pick],(strain-strain_0).reshape((nt,9))[:,pick])*vol_mean/kbT/nt*1e9
    CS=np.linalg.inv(CS_minus1)
    with np.printoptions(precision=3, suppress=False,linewidth=150):
        print('xx, yy, zz, xy, xz, yz : GPa')
        print (CS)
        print('============')
        print('1/GPa')
        print(np.linalg.det(CS_minus1))
        print(CS_minus1)
    if is_diagonal: #fill with zeros the 6x6 matrix not involved in the calculation
        CS2=np.zeros((6,6))
        CS_minus12=np.zeros((6,6))
        CS2[:3,:3]=CS
        CS_minus12[:3,:3]=CS_minus1
        CS=CS2
        CS_minus1=CS_minus12
    return CS, CS_minus1

def hist2gofr(gr_N,gr_dr,gr_0,gofr):
    rs_m=np.arange(gr_N)*gr_dr+gr_0
    rs_p=(np.arange(gr_N)+1)*gr_dr+gr_0
    vols=4*np.pi/3*(rs_p**3-rs_m**3)
    return gofr/vols

def peak_width(gr,NH,NH_thre,gr_r2i,min_spread=0.1):
    '''
    get the width of the first peak in the gr at a given height
    '''
    idx_23=0
    while gr[0,NH,idx_23] < NH_thre:
        idx_23 += 1
    idx_spread_low=idx_23
    while gr[0,NH,idx_23] > NH_thre*0.9:
        if idx_23 == gr.shape[2]-1 or (idx_23 -idx_spread_low > gr_r2i(min_spread) and gr[0,NH,idx_23] < NH_thre):
            break
        idx_23 += 1
    idx_spread_hi=idx_23
    return idx_spread_low,idx_spread_hi

def analyze_peak(gr,r0,peak_fractional_height,gr_r2i,gr_i2r):
    '''
    find the peak around r0 in the g(r) array and get some info
    '''
    maxs=np.argmax(gr,axis=2)
    peaks=gr_i2r(maxs)
    
    #find nearest peak to r0
    NH=np.argmin((peaks-r0)**2)
    #get spread at 'peak_fractional_height' of peak height
    NH_peak_idx=maxs[0,NH]
    NH_peak_val=gr[0,NH,NH_peak_idx]
    NH_thre=NH_peak_val*peak_fractional_height
    #get width of peak
    idx_spread_low,idx_spread_hi=peak_width(gr,NH,NH_thre,gr_r2i)
    spread_low=gr_i2r(idx_spread_low)
    spread_hi=gr_i2r(idx_spread_hi)
    w23h=spread_hi-spread_low
    return idx_spread_low,idx_spread_hi,spread_low,spread_hi,w23h,NH_thre,NH_peak_val,NH_peak_idx,peaks,NH

def get_conv_functs(gr_0,gr_dr):
    def gr_r2i(r):
        return int((r-gr_0)/gr_dr)
    def gr_i2r(i):
        return gr_0+i*gr_dr
    return np.vectorize(gr_r2i), np.vectorize(gr_i2r)
    

def do_compute_gr_multi(t,gr_kw={'tskip':10},n_segments=1,gr_0=0.5,gr_end=3.8,gr_N=150):
    '''
    computes the g(r) pair correlation function and find the peaks of it
    '''
    nts=t.getNtimesteps()
    param_gr=(nts,gr_0,gr_end,gr_N)
    #first calculate g(r), get the N-H peak, estabilish the range of sh correlation f
    gr_dr=(gr_end-gr_0)/gr_N

    #get utility functions for converting from index to r value and back
    gr_r2i,gr_i2r=get_conv_functs(gr_0,gr_dr)
    
    gofr=analyze_gofr(t,0,nts,gr_0,gr_end,gr_N,
                      tmax=1,**gr_kw,n_segments=n_segments) #only g(r)
    #from the histogram generate the g(r) -- divide by the volumes of the spherical shells
    res=[]
    for i in range(n_segments):
        gr=hist2gofr(gr_N,gr_dr,gr_0,gofr if n_segments==1 else gofr[i])
        
            
        #find nearest peak to 1.0 (N-H bond)
        idx_spread_low,idx_spread_hi,spread_low,spread_hi,w23h,NH_thre,NH_peak_val,NH_peak_idx,peaks,NH = analyze_peak(gr,1.0,0.5,gr_r2i,gr_i2r)
        
        #print(f'width at 1/2 of height: {w23h}')
        #calculate sh correlations of hydrogen peak
        sh_low=spread_low-w23h*0.1
        sh_hi=spread_hi+w23h*0.1
        NH_peak=(NH_thre,spread_low,spread_hi,peaks,NH,w23h,NH_peak_val)
        res.append( (None,param_gr,gofr if n_segments==1 else gofr[i],(None,None),None,NH_peak))
    return res


def do_compute_gr_sh(t,times,do_sh=True,neigh=[],analyze_sh_kw={'tskip':50},gr_kw={'tskip':10}):
    '''
    computes the g(r) pair correlation function and the spherical harmonics correlation function computed around the N-H peak of the g(r), at around 1.0 Angstrom
    '''
    gr_0=0.5
    gr_end=3.8
    gr_N=150
    nts=t.getNtimesteps()
    param_gr=(nts,gr_0,gr_end,gr_N)
    DT_PS=times[1]-times[0]
    #first calculate g(r), get the N-H peak, estabilish the range of sh correlation f
    gr_dr=(gr_end-gr_0)/gr_N

    #get utility functions for converting from index to r value and back
    gr_r2i,gr_i2r=get_conv_functs(gr_0,gr_dr)
    
    gofr=analyze_gofr(t,0,nts,gr_0,gr_end,gr_N,
                      tmax=1,**gr_kw) #only g(r)
    #from the histogram generate the g(r) -- divide by the volumes of the spherical shells
    gr=hist2gofr(gr_N,gr_dr,gr_0,gofr)
    
        
    #find nearest peak to 1.0 (N-H bond)
    idx_spread_low,idx_spread_hi,spread_low,spread_hi,w23h,NH_thre,NH_peak_val,NH_peak_idx,peaks,NH = analyze_peak(gr,1.0,0.5,gr_r2i,gr_i2r)
    
    #print(f'width at 1/2 of height: {w23h}')
    #calculate sh correlations of hydrogen peak
    sh_low=spread_low-w23h*0.1
    sh_hi=spread_hi+w23h*0.1
    if do_sh:
        sh=analyze_sh(t,0,nts,sh_low,sh_hi,1,tmax=int(0.5/DT_PS),sann=neigh,**analyze_sh_kw)
    else:
        sh=None
    NH_peak=(NH_thre,spread_low,spread_hi,peaks,NH,w23h,NH_peak_val)
    param_sh = (sh_low,sh_hi)
    return times,param_gr,gofr,param_sh,sh,NH_peak

def do_plots_gr_sh(times,param_gr,gofr,param_sh,sh,NH_peak,fig_ax=None,plot_gr_kw={}):

    nts,gr_0,gr_end,gr_N=param_gr
    sh_low,sh_hi=param_sh
    NH_thre,spread_low,spread_hi,peaks,NH,w23h,NH_peak_val=NH_peak

    gr_dr=(gr_end-gr_0)/gr_N
    #get utility functions for converting from index to r value and back
    gr_r2i,gr_i2r=get_conv_functs(gr_0,gr_dr)
    
    fig_gr,ax_gr=plot_gofr(gr_0,gr_end,gofr,fig_ax=fig_ax,**plot_gr_kw)
    ax_gr.grid()
    if sh_low is not None: 
        ax_gr.axvline(sh_low,color='r')
    if sh_hi is not None:
        ax_gr.axvline(sh_hi,color='r')
    
    #annotate N-N peak
    
    #from the histogram generate the g(r) -- divide by the volumes of the spherical shells
    gr=hist2gofr(gr_N,gr_dr,gr_0,gofr)
    def annotate_peak(r0,h,ax_gr):
        _,_,NN_spread_low,NN_spread_hi,w2NN,NN_thre,NN_peak_val,NN_peak_idx,peaks,NN = analyze_peak(gr,r0,h,gr_r2i,gr_i2r)
        ax_gr.hlines(NN_thre,NN_spread_low,NN_spread_hi)
        ax_gr.annotate(f'{w2NN:.2f}',(peaks[0,NN]*0.95,NN_thre*1.05))
        ax_gr.annotate(f'{peaks[0,NN]:.2f}',(peaks[0,NN]*0.95,NN_peak_val*1.02))
        return w2NN,peaks[0,NN]
    #find nearest peak to 2.6 (N-N bond)
    wNN,rNN=annotate_peak(2.6,0.5,ax_gr)
    #find nearest peak to 1.7 (H-H bond)
    wHH,rHH=annotate_peak(1.7,3.0/4,ax_gr)
    #find nearest peak to 1.0 (N-H bond)
    wNH,rNH=annotate_peak(1.0,0.5,ax_gr)

    if sh is not None:
        try:
            fig_sh,ax_sh,axins_sh,fit_sh=plot_sh(sh_low,sh_hi,times,sh,0,1,0,log=False,pre_fit=0.4)
        except:
            try:
                fig_sh,ax_sh,axins_sh,fit_sh=plot_sh(sh_low,sh_hi,times,sh,0,1,0,log=False,pre_fit=0.1)
            except:
                fig_sh,ax_sh,axins_sh,fit_sh=plot_sh(sh_low,sh_hi,times,sh,0,1,0,log=False)
        ax_sh.grid()
    else:
        fig_sh,ax_sh,axins_sh,fit_sh=None,None,None,None
    return fig_gr,ax_gr,fig_sh,ax_sh,axins_sh,fit_sh,wNN,rNN,wHH,rHH,wNH,rNH

def do_compute_msd(t_unw,times,msd_kw={}):
    nts=t_unw.getNtimesteps()
    msd=analyze_msd(t_unw,0,nts,**msd_kw)
    #compute slope
    DT_PS=times[1]-times[0]
    msd_start=max(min(int(0.3/DT_PS),msd.shape[0]-100),0)
    coeffs=np.polyfit(times[msd_start:msd.shape[0]]-times[0],msd[msd_start:,0,:],1)
    return msd,times,msd_start,coeffs

def do_plots_msd(msd,times,msd_start,coeffs):
    DT_PS=times[1]-times[0]
    fig,ax=plot_msd(times,msd,0)
    t0=times[msd_start]-times[0]
    tf=times[msd.shape[0]-1]-times[0]
    lines=[[(t0,coeffs[0,0]*t0+coeffs[1,0]),(tf,coeffs[0,0]*tf+coeffs[1,0])],
          [(t0,coeffs[0,1]*t0+coeffs[1,1]),(tf,coeffs[0,1]*tf+coeffs[1,1])]]
    lc=mc.LineCollection(lines,zorder=-1)
    ax.add_collection(lc)

    ax.annotate(f'D={coeffs[0,0]/6:.2f}',lines[0][1])
    ax.annotate(f'D={coeffs[0,1]/6:.2f}',lines[1][1])
    return fig,ax

def inspect(traj, only_cell=False,plot_traj=True,plot=True,
            do_sh=True,do_density=False, plot_sh_h=True,
            show_traj_dyn=False,dt=None,
            neigh=DEFAULT_NEIGH,
            plot_st_kw={'transpose':True,'xmax':.20,'ymax':.60},
            analyze_sh_kw={'tskip':50},
            compute_steinhardt_kw={'skip':10},
            msd_kw={},
            gr_kw={'tskip':10},
            nthreads=4):
    results={}
    analyze_sh_kw['nthreads']=nthreads
    compute_steinhardt_kw['nthreads']=nthreads
    msd_kw['nthreads']=nthreads
    gr_kw['nthreads']=nthreads


    natoms=traj.numsites
    #conversion factor from eV to K
    k_b=8.617333262145e-5 #eV/K
    eV_to_K=2/(3*k_b*natoms)

    #print('traj arraynames = {}'.format(traj.get_arraynames()))
    if 'ionic_temperature' in traj.get_arraynames():
        temp=traj.get_array('ionic_temperature')
        T_mean=temp.mean()
        results['T']=T_mean
    else:
        T_mean=float('nan')
    if 'pressure' in traj.get_arraynames():
        press=traj.get_array('pressure')
        press_mean=press.mean()
        results['P']=press_mean
    else:
        press_mean = float('nan')
    
    if 'times' in traj.get_arraynames():
        t=traj.get_array('times')
    else:
        t=np.arange(temp.shape[0])*dt
    def t_to_timestep(x):
        return np.interp(x,t,np.arange(t.shape[0]))
    def timestep_to_t(x):
        return np.interp(x,np.arange(t.shape[0]),t)
        
    
    DT_FS=(t[1]-t[0])*1e3
    cell_0 =traj.get_array('cells')
    volume=np.linalg.det(cell_0)
    #print('cell(t=0) = {}'.format(cell_0))
    #print('cell_volume(t=0) = {} A^3'.format(volume))
    cell_0=traj.get_array('cells').mean(axis=0)
    volume=np.linalg.det(cell_0)
    results['V']=volume
    if 'electronic_kinetic_energy' in traj.get_arraynames():
        results['ekinc']=traj.get_array('electronic_kinetic_energy').mean()
    if 'energy_constant_motion' in traj.get_arraynames(): 
        results['e_constant']=traj.get_array('energy_constant_motion').mean()
    results['DT_FS']=DT_FS

    plt_fname_pre=FIGURE_PATH+'/inspect_result'
    if 'ionic_temperature' in traj.get_arraynames() and 'pressure' in traj.get_arraynames():
        display(HTML(f'<h1>{T_mean:.0f}K {press_mean:.0f}GPa</h1>'))
        print('DT_FS = {:.3f}fs, T = {:.2f}K, P = {:.2f}GPa'.format(DT_FS,T_mean,press_mean))
        plt_fname_pre=FIGURE_PATH+f'/{T_mean:.0f}K_{press_mean:.0f}GPa_{DT_FS:.3f}_'
    print('traj pk = {}'.format(traj.pk))
    plt_fname_suff='.pdf'
    if traj.pk != None:
        pickle_dump_name=f'traj_{traj.pk}'
    else:
        pickle_dump_name=f'traj_{T_mean:.2f}K{press_mean:.2f}GPa'


    if plot_traj and plot:
        print('cell_average = {}'.format(cell_0))
        print('cell_average_volume = {} A^3'.format(volume))
        plt_key(traj,'electronic_kinetic_energy',eV_to_K,ylabel='K')
        plt_key(traj,'energy_constant_motion')
        plt_key(traj,'ionic_temperature',ylabel='K')
        plt_key(traj,'pressure',ylabel='GPa')
        plt.show()
    atraj,atraj_unw=get_analisi_traj_from_aiida(traj)


    if show_traj_dyn:
        traj_dyn=show_atraj(traj,atraj,wrap=False)
        
    
    if not only_cell:

        
        if plot and do_density:
            res=atomic_density(atraj)
            plot_=density_field(*res)
        
        #histogram of steinhardt parameters
        if plot_sh_h: 
            sh_h_pickle=pickle_dump_name+'_sh_h.pickle'
            sh_h = pickle_or_unpickle(sh_h_pickle)
            if sh_h is None:
                sh_h=compute_steinhardt(atraj,neigh=neigh,averaged=True,**compute_steinhardt_kw)
                fig_shh,axs_shh=plt_steinhardt(sh_h,**plot_st_kw)
                fig_shh.savefig(plt_fname_pre+'steinhardt'+plt_fname_suff)
            

        #g of r / spherical correlations
        grsh_pickle=pickle_dump_name+ ('_gofr_sh.pickle' if do_sh else '_gofr.pickle')
        gofrsh = pickle_or_unpickle(grsh_pickle)
        if gofrsh is None:
            gofrsh=do_compute_gr_sh(atraj,t,do_sh=do_sh,neigh=neigh,analyze_sh_kw=analyze_sh_kw,gr_kw=gr_kw)
            pickle_or_unpickle(grsh_pickle,analisi = gofrsh)
                
        #msd
        msd_pickle=pickle_dump_name+'_msd.pickle'
        msd = pickle_or_unpickle(msd_pickle)
        if msd is None:
            msd=do_compute_msd(atraj_unw,t,msd_kw=msd_kw)
            pickle_or_unpickle(msd_pickle,analisi = msd)
        
        wNN,rNN,wHH,rHH,wNH,rNH=None,None,None,None,None,None
        if plot:
            fig_gr,ax_gr,fig_sh,ax_sh,axins_sh,fit_sh,wNN,rNN,wHH,rHH,wNH,rNH=do_plots_gr_sh(*gofrsh)
            fig_gr.show()
            fig_gr.savefig(plt_fname_pre+'gr'+plt_fname_suff)
            if fig_sh is not None:
                fig_sh.show()
                fig_sh.savefig(plt_fname_pre+'sh'+plt_fname_suff)
            fig_msd,ax_msd=do_plots_msd(*msd)
            fig_msd.show()
            fig_msd.savefig(plt_fname_pre+'msd'+plt_fname_suff)
        
        results['msd']=msd[3][0,:]
        
        NH_thre,spread_low,spread_hi,peaks,NH,w23h,NH_peak_val=gofrsh[5]
        results['NH_peak_width']=w23h
        results['NH_peak_pos']=peaks[0,NH]
        results['NN_peak_width']=wNN
        results['NN_peak_pos']=rNN
        results['HH_peak_width']=wHH
        results['HH_peak_pos']=rHH

    cells_transition=atraj.get_box_copy()
    
    if (cells_transition!=cells_transition[0]).any():
        if plot:
            fig, ax = plt.subplots(nrows=2,figsize=(10,8),dpi=300)
            lines=ax[0].plot(t,cells_transition[:,3:6]*2)
            for i,c in enumerate(['x','y','z']):
                lines[i].set_label(c)
            ax[0].plot(t,(2*cells_transition[:,3:6]).prod(axis=1)**.33333,label=r'$volume^{\frac{1}{3}}$')
            ax[0].grid()
            ax[0].set_title('cell parameters: cell size (up) and cell tilt (down)')
            ax[1].set_xlabel('t (ps)')
            ax[0].set_ylabel('(A)')
            ax[0].secondary_xaxis('top',functions=(t_to_timestep,timestep_to_t))
            ax[0].legend()
            if cells_transition.shape[1]>6:
                lines=ax[1].plot(t,cells_transition[:,6:])
                ax[1].grid()
                ax[1].set_ylabel('(A)')
                ax[1].secondary_xaxis('top',functions=(t_to_timestep,timestep_to_t))
                for i,c in enumerate(['xy','xz','yz']):
                    lines[i].set_label(c)
                ax[1].legend()
            fig.show()
            fig.savefig(plt_fname_pre+'cell'+plt_fname_suff)
        try:
            uma=1.66e-27 #kg
            mcell=(sum( [ 1 for i in traj.get_attribute('symbols') if i=='H']) + sum( [ 14 for i in traj.get_attribute('symbols') if i=='N']))*uma
            def vs_C(traj,nsteps):
                density=(mcell/np.linalg.det(traj.get_array('cells')[:nsteps]*1e-10)).mean()
                CS,CSm=elastic_c(traj.get_array('ionic_temperature')[:nsteps].mean()*1.38064852e-23,
                                 traj.get_array('cells')[:nsteps]*1e-10)
                vs=((CS[0,0]+4.0/3*CS[3,3])*1e9/density)**.5
                print ('sqrt((C_{xx,xx} + C_{xy,xy})/density) [m/s] ' ,vs)
                return CS,CSm,vs,density
            CSs,CSms,vss,ts,densities=([],[],[],[],[])
            for i in range(10):
                stop_step=traj.numsteps*(i+1)//10
                ts.append(t[stop_step-1])
                CSi,CSmi,vsi,density=vs_C(traj,stop_step)
                CSs.append(CSi)
                CSms.append(CSmi)
                vss.append(vsi)
                densities.append(density)
            CSs=np.array(CSs)
            CSms=np.array(CSms)
            vss=np.array(vss)
            ts=np.array(ts)
            if plot:
                fig, ax = plt.subplots(figsize=(10,8),dpi=300)
                ax.plot(ts,vss)
                ax.set_xlabel('t (ps)')
                ax.set_ylabel('vs (m/s)')
                ax.grid()
                fig.show()
                fig.savefig(plt_fname_pre+'vs'+plt_fname_suff)
            results['Vs']=vss
            results['CS_GPa']=CSs
            results['density']=densities[-1]
        except:
            print('error in estimating elastic constant')
    return atraj,atraj_unw, results

def multiinspect(nodes,plot=False,prefix='',inspect_kw={}):
    all_res=[]
    for node in nodes:
        _,_,res = inspect(node,plot=plot,**inspect_kw)
        plt.show()
        all_res.append(res)
    return all_res, print_all(all_res,prefix=prefix)
def print_all(all_res,prefix=''):
    Ts=[]
    Ps=[]
    MSDs=[]
    Vss=[]
    DTs=[]
    CS_GPas=[]
    rhos=[]
    grp_pos=[]
    grp_width=[]
    for res in all_res:
        if 'T' in res:
            Ts.append(res['T'])
        if 'P' in res:
            Ps.append(res['P'])
        MSDs.append(res['msd'])
        DTs.append(res['DT_FS'])
        if 'Vs' in res:
            Vss.append(res['Vs'][-1])
            CS_GPas.append(res['CS_GPa'][-1])
        if 'density' in res:
            rhos.append(res['density'])
        grp_pos.append([res['NH_peak_pos'],res['NN_peak_pos'],res['HH_peak_pos']])
        grp_width.append([res['NH_peak_width'],res['NN_peak_width'],res['HH_peak_width']])
    MSDs=np.array(MSDs)
    CS_GPas=np.array(CS_GPas)
    rhos=np.array(rhos)
    grp_pos=np.array(grp_pos)
    grp_width=np.array(grp_width)
    alist=[('NN peak r',grp_pos[:,1]),('NN peak width',grp_width[:,1]),('msd0',MSDs[:,0]),('msd1',MSDs[:,1])] + ([('Vs',Vss)] if len(Vss) > 0 else [])
    try:
        for k,arr in alist:
            fig, ax = plt.subplots(nrows=1,figsize=(10,8),dpi=300)
            for i in range(len(Ts)):
                 ax.annotate(f'dt={DTs[i]:.3f}', xy=(Ts[i], arr[i]),
                     xytext=(-1, 1), textcoords="offset points",
                     horizontalalignment="right")
            ax.scatter(Ts,arr,label=k)
            ax.set_xlabel('T (K)')
            ax.set_ylabel(k)
            ax.grid()
            fname=FIGURE_PATH+'/summary_'+prefix+k.replace(' ', '_')+'.pdf'
            fig.show()
            fig.savefig(fname)
        #some sqrt(elastic constants/density)
        def plt_el(Ts,CS,density):
            fig, ax = plt.subplots(nrows=1,figsize=(10,8),dpi=300)
#            for arr,label in [(((CS[:,0,0]+CS[:,1,1]+CS[:,2,2])/(3.0*density))**.5,'longitudinal'),(((CS[:,3,3]+CS[:,4,4]+CS[:,5,5])/(3.0*density))**.5,'shear')]:
            for arr,label in [((CS[:,0,0]/density)**.5,'longitudinal'),((CS[:,3,3]/density)**.5,'shear')]:
                for i in range(len(Ts)):
                     ax.annotate(f'dt={DTs[i]:.3f}', xy=(Ts[i], arr[i]),
                         xytext=(-1, 1), textcoords="offset points",
                         horizontalalignment="right")
                ax.scatter(Ts,arr,label=label)
            ax.set_xlabel('T (K)')
            ax.set_ylabel('v (m/s)')
            ax.legend()
            ax.grid()
        if len(CS_GPas) > 0:
            plt_el(Ts,CS_GPas*1e9,rhos)
    except:
        pass
    
    return Ts,Ps,MSDs,Vss


