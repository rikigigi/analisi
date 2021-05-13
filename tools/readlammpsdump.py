# coding=utf-8

import numpy as np
import sys

def read_natoms_dump(infile):
	with open(infile,'r') as inf :
		l=inf.readline()
		l=inf.readline()
		l=inf.readline()
		natoms=int(inf.readline())
	return natoms
def read_cols_to_read(infile,pos=True,vel=True,force=False):
	result =[]
	pos_col=[None,None,None]
	vel_col=[None,None,None]
	for_col=[None,None,None]
	with open(infile,'r') as inf :
		for i in range(9):
			l=inf.readline()
		l.split()
		l1 = l[2:]
		if(pos):
			for i,item in enumerate(l1):
				if(item[0]=='x'):
					pos_col[0]=i
				elif(item[0]=='y'):
					pos_col[1]=i
				elif(item[0]=='z'):
					pos_col[2]=i
		if(vel):
			for i,item in enumerate(l1):
				if(item=='vx'):
					vel_col[0]=i
				elif(item=='vy'):
					vel_col[1]=i
				elif(item=='vz'):
					vel_col[2]=i
		if(force):
			for i,item in enumerate(l1):
				if(item=='fx'):
					force_col[0]=i
				elif(item=='fy'):
					force_col[1]=i
				elif(item=='fz'):
					force_col[2]=i
	
	return pos_col,vel_col,force_col
def get_types(ifile,natoms=None):
	types = []
	types_col = None
	if (natoms is None) : natoms = read_natoms_dump(ifile)
	with open(infile,'r') as inf :
		for i in range(9):
			l=inf.readline()
		l.split()
		for icol,item in l:
			if item == 'type':
				types_col = icol
		for iatoms in range(natoms):
			l=inf.readline().split()
			types.append(int(l[types_col]))	
	return np.array(types,dtype=np.int64)

def readlammpsdump(infile,nstep,iskip=0,every=1,natoms=None,pos=True,vel=True, box=True , force = False):
	if(natoms is None):
		natoms = read_natoms_dump(infile)
	pos_col,vel_col,force_col = read_cols_to_read(infile,pos=True,vel=True,force=False)
	traj ={}
	nstep_tosave = nstep//every
	if (pos):
		traj['p']=np.zeros(nstep_tosave,natoms,3)
	if (vel):
		traj['v']=np.zeros(nstep_tosave,natoms,3)
	if (force) :
		traj['f']=np.zeros(nstep_tosave,natoms,3)
	#traj['box']=np.zeros(nstep_tosave,3,3)
	traj['box']=np.zeros(nstep_tosave,6)
	traj['step']=np.zeros(nstep_tosave,dtype=np.int64)
	traj['types']= get_types(infile,natoms)
	with open(infile,'r') as inf :
		for i in np.arange(iskip*(natoms+9)):
			l=inf.readline()
		itimestep = -1
		for i in np.arange(1,nstep):
			iread=i%every	
			if(iread==0):
				itimestep +=1
				l=inf.readline() #ITEM: TIMESTEP
				l=inf.readline() #timestep
				traj['step'][itimestep]=int(l)
				l=inf.readline() #ITEM: NUMBER OF ATOMS
				l=inf.readline() #n atoms
				l=inf.readline() # ITEM: BOX

				l=inf.readline().split()
				box[itimestep,0] = float(l[0])
				box[itimestep,1] = float(l[1])
				l=inf.readline().split()
				box[itimestep,2] = float(l[0])
				box[itimestep,3] = float(l[1])
				l=inf.readline().split()
				box[itimestep,1,4] = float(l[1])-float(l[0])
				box[itimestep,2,5] = float(l[1])-float(l[0])
				l=inf.readline() # ITEM: ATOMS
				for iatom in range(natoms):
					l=[float(j) for j in inf.readline().split()]
					if(pos):
						traj['p'][itimestep,iatom,:] = l[pos_col]
					if(vel):
						traj['v'][itimestep,iatom,:] = l[vel_col]
					if (force) :
						traj['f'][itimestep,iatom,:] = l[force_col]
			else:
				for ii in np.arange(natoms+9):
					l=inf.readline()
				

	return traj
