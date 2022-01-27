import numpy as np

def curved(nsteps,r,a1,a2,tmax):
    t = np.linspace(0,a2*tmax,nsteps,endpoint=False)
    a = np.linspace(0,a1*tmax,nsteps,endpoint=False)
    x=r*np.cos(t)*np.sin(a)
    y=r*np.sin(t)*np.sin(a)
    z=r*np.cos(a)
    return np.c_[x,y,z]

def cell_traj(a,b,c,n):
    l=[]
    for i in range(n):
        l.append([a,b,c])
    return np.copy(np.array(l,dtype=np.double).transpose((0,2,1)),order='C' )

    
