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
#traj.write_lammps_binary("lammps.bin")
test=pyanalisi.ShpericalCorrelations(traj,0.5,2.5,2,100,4,3,False)
test.reset(400)
test.calculate(3000)
print(np.array(test,copy=False))

