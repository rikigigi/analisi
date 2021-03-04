from aiida import load_profile
from aiida.orm import Code, load_node
from aiida.orm import *
from aiida.engine import submit
import re
load_profile()


from aiida_QECpWorkChain.workflow import *
from aiida_QECpWorkChain import write_cp_traj

import numpy as np
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000000
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pickle

try:
    import thermocepstrum as tc
except:
    print('WARINING: cannot import thermocepstrum')


import pyanalisi as pa

print(pa.info())

def plt_key(traj,key,conv=1.0,title='',ylabel=''):
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

def get_analisi_traj_from_aiida(traj):
    pos=traj.get_array('positions')
    vel=traj.get_array('velocities')
    cel=traj.get_array('cells')
    types=get_types_id_array(traj.get_attribute('symbols'))
    params= [pos, vel, types,  cel]
    atraj=pa.Trajectory(*params,True, True)
    atraj_unw=pa.Trajectory(*params,True, False)
    return atraj, atraj_unw

import k3d
from ipywidgets import interact, interactive, fixed
import ipywidgets as widgets

import scipy.ndimage

def analyze_msd(traj_unw,start,stop,nthreads=4,tskip_msd=10):
    l=stop-start
    msd_l=l//2
    
    #mean square displacement
    msd=pa.MeanSquareDisplacement(traj_unw,tskip_msd,msd_l,nthreads,True,False,False)
    msd.reset(l-msd_l)
    print('calculating msd...',flush=True)
    msd.calculate(start)
    return msd

def analyze(traj,traj_unw,start,stop,nthreads=4,tskip_msd=10,tskip_sh=200,t_sh=None):
    l=stop-start
    msd_l=l//2
    
    msd=analyze_msd(traj_unw,start,stop,nthreads,tskip_msd)
    
    #gofr (with time)
    startr=0.5
    endr=3.8
    nbin=100
    gofr=pa.Gofrt(traj,startr,endr,nbin,1,nthreads,tskip_msd,False)
    gofr.reset(l-1)
    print('calculating g(r)...',flush=True)
    gofr.calculate(start)
    
    #spherical harmonics correlations
    startr=0.7
    endr=3.5
    nbin=4
    dr=(endr-startr)/nbin
    sh=pa.SphericalCorrelations(traj
                                     ,startr #minimum of radial distance
                                     ,endr #maximum of radial distance
                                     ,nbin # number of distance bins
                                     ,msd_l if t_sh is None else t_sh # maximum length of correlation functions in timesteps
                                     ,nthreads # number of threads
                                     ,tskip_sh # time skip to average over the trajectory
                                     ,20 # buffer for sh
                                     ,False)
    sh.reset(l-msd_l) # number of timesteps to analyze
    print('calculating spherical harmonics correlation functions... (go away and have a coffee)',flush=True)
    sh.calculate(start) #calculate starting at this timestep
    print('all done',flush=True)
    ress= np.array(msd,copy=True),np.array(gofr,copy=True),np.array(sh,copy=True)
    return ress

import scipy as sp
import scipy.optimize 
import math

np.ALLOW_THREADS=4

sp.ALLOW_THREADS=4

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
    print(s)
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
    print(s)
    print(s_aux)
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
    print(s)
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

def plot_msd(times,res,cm,title='',res_var=None):
    fig,ax =plt.subplots(figsize=(10,8),dpi=300)
    ax=fig.add_axes([0,0,1,1])
    for i in range(res.shape[-1]):
        plt_err(ax,times[:res.shape[0]]-times[0],res[:,cm,i],res_var[:,cm,i] if res_var is not None else res_var,label='type='+str(i))
    ax.legend()
    ax.set_title('{}MSD'.format(title))
    ax.set_xlabel('time (ps)')
    ax.set_ylabel('$\AA^2/ps$')
    return ax
    
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
            maxt=-1.0,fitmax=-1.0,fitmin=-1.0,pre_fit=-1.0,log=True,rescale_array=True,scale_exp_coeff=0.1):
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
    return ax,axins, fits


def plot_gofr(startr,endr,res,title='',res_var=None):
    fig,ax =plt.subplots(figsize=(10,8),dpi=300)
    r=np.arange(res.shape[-1])*(endr-startr)/res.shape[-1]+startr
    for i in range(0,res.shape[1]//2):
        plt_err(ax,r,res[0,i,:]/(r**2),(res_var[0,i,:]/(r**4)) if res_var is not None else res_var)
    ax.set_xlabel('$\AA$')
    ax.set_ylabel('number of atoms in the shell / $r^2$')
    ax.set_title('{}$g(r)$'.format(title))
    return ax

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

def density_field(res_g, res_g1):
    plot = k3d.plot(grid=(-0.5,-0.5,-0.5,0.5,0.5,0.5))
    obj = k3d.volume(res_g)
    obj2 = k3d.volume(res_g1)
    #points = k3d.points(g0_sites,point_size=0.02)
    #points1 = k3d.points(g1_sites,point_size=0.005)
    plot += obj
    plot += obj2
    #plot += points
    #plot += points1
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
    def global_clipping_planes():
        if _sphere:
            return sph_tangents_cover(_w_x,7,_x,_y,_z)
        else:
            return [[1,0,0,_w_x-_x],[-1,0,0,_w_x+_x],
                [0,1,0,_w_y-_y],[0,-1,0,_w_y+_y],
                [0,0,1,_w_z-_z],[0,0,-1,_w_z+_z]
               ]
    button=widgets.Button(description='Center on 1')
    @button.on_click
    def button1_on_click(b):
        global _x
        global _y
        global _z
        print('finding maximum in spherical domain with r={}, center=({},{},{})'.format(_w_x,_x,_y,_z))
        idx, (_z,_y,_x)=get_max_in_spherical_domain(res_g,_w_x,_z,_y,_x)
        #_x,_y,_z=-_x,-_y,-_z
        print('center=({},{},{}), idx={} shape={}'.format(_x,_y,_z,idx,res_g.shape))
        print('value={}'.format(res_g[idx]))
        plot.clipping_planes=global_clipping_planes()
    display(button)#,button2,button3)
    @interact(spherical=widgets.Checkbox(value=False))
    def g(spherical):
        global _sphere
        _sphere=spherical
        plot.clipping_planes=global_clipping_planes()
    @interact(x=widgets.FloatSlider(value=0.14,min=-.5,max=.5,step=0.01))
    def g(x):
        global _x
        _x=x
        plot.clipping_planes=global_clipping_planes()
    @interact(w_x=widgets.FloatSlider(value=0.05,min=0.0,max=.5,step=0.01))
    def g(w_x):
        global _w_x
        _w_x=w_x
        plot.clipping_planes=global_clipping_planes()
    @interact(y=widgets.FloatSlider(value=0.0,min=-.5,max=.5,step=0.01))
    def g(y):
        global _y
        _y=y
        plot.clipping_planes=global_clipping_planes()
    @interact(w_y=widgets.FloatSlider(value=0.5,min=0.0,max=.5,step=0.01))
    def g(w_y):
        global _w_y
        _w_y=w_y
        plot.clipping_planes=global_clipping_planes()
    @interact(z=widgets.FloatSlider(value=0.0,min=-.5,max=.5,step=0.01))
    def g(z):
        global _z
        _z=z
        plot.clipping_planes=global_clipping_planes()
    @interact(w_z=widgets.FloatSlider(value=0.5,min=0.0,max=.5,step=0.01))
    def g(w_z):
        global _w_z
        _w_z=w_z
        plot.clipping_planes=global_clipping_planes()
        
def force_ratio_histogram(wf,print=print,ax=[]):
    pwcalcjobs=[]
    for c in [x for x in wf.called if str(x.process_class) == "<class 'aiida_quantumespresso.calculations.pw.PwCalculation'>"]:
        pwcalcjobs.append(c)
        print(c.pk)
    res,axs,figs,n_ax=analyze_forces_ratio(pwcalcjobs,minpk=wf.pk,ax_=ax,create_fig = lambda : plt.subplots(dpi=300))
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

def print_cp_with_traj(wf,print=print):
    def length_traj(traj):
        t=traj.get_array('times')
        return t[-1]-t[0]
    for c in [x for x in wf.called if str(x.process_class) == "<class 'aiida_quantumespresso.calculations.cp.CpCalculation'>"]:
        print ("===========")
        print(c.pk,'traj pk={}: {:.3}ps'.format(c.outputs.output_trajectory.pk,length_traj(c.outputs.output_trajectory)) if 'output_trajectory' in c.outputs else 'No output_trajectory')
        print(c.inputs.parameters.get_dict()['IONS'])
        print(c.inputs.parameters.get_dict()['CELL'] if 'CELL' in c.inputs.parameters.get_dict() else '')

def atomic_density(atraj,cell,dr=0.1):
    hist=pa.PositionHistogram(atraj,np.array(np.rint(cell/dr),dtype=int).diagonal().tolist(),1,1)
    hist.reset(atraj.getNtimesteps())
    hist.calculate(0)
    res=np.array(hist)
    return res

def fft_density(res,coeff=np.array([1.0,1.0])):
    resg=np.einsum('ijkl,i',np.fft.fftn(res,axes=(1,2,3)),coeff)
    resg[0,0,0]=0
    resg=np.fft.fftshift(resg)
    resgn=np.real(resg*np.conj(resg))
    resgn/=resgn.max()
    return resgn