import numpy as np
import pyanalisi
filepath_tests = '/home/bertossa/analisi/tests'
pos = np.load(filepath_tests + '/data/positions.npy')  
vel = np.load(filepath_tests + '/data/velocities.npy') 
box = np.load(filepath_tests + '/data/cells.npy')      
types = np.zeros(pos.shape[1], dtype = np.int32)            
types[-8:]=1                                           
params= [pos, vel, types, box]
traj=pyanalisi.Trajectory(*params,True, True)
test=pyanalisi.ShpericalCorrelations(traj,0.5,5.0,2,10,4,1,False)
test.reset(100)
test.calculate(0)
print(np.array(test,copy=False))

