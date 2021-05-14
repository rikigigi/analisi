#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import sys

def read_natoms_dump(infile):
	"""
	get the number of atoms from a lammps dump file
	
	Parameters
	---------
	infile : str
		input file
	
	Returns
	-------
	natoms : int
		number of atoms
	"""
	with open(infile,'r') as inf :
		l=inf.readline()
		l=inf.readline()
		l=inf.readline()
		natoms=int(inf.readline())
	return natoms
def read_cols_to_read(infile,pos=True,vel=True,force=False):
	"""
	read the cols of positions, velocities and forces
	Parameters
	---------
	infile : str
		input file
	pos : bool
		read the positions
	vel : bool 
		read the velocities
	force :
		read the forces
	
	Returns
	------
	pos_col,vel_col,force_col : list(3)
		list of the indeces of the columns in the lammps file 
		if the specific cols do not exist then returns None
	
	"""
	result =[]
	pos_col=[None,None,None]
	vel_col=[None,None,None]
	force_col=[None,None,None]
	with open(infile,'r') as inf :
		for i in range(9):
			l=inf.readline()
		l1 = l.split()[2:]
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
	"""
	get the types array
	
	Parameters
	---------
	ifile : str
		input file

	natoms : int
		number of atoms. if None it reads the number from the traj file
	
	Returns
	-------
	types : np.array(dtype=np.int32)
	"""
	types = []
	types_col = None
	if (natoms is None) : natoms = read_natoms_dump(ifile)
	with open(ifile,'r') as inf :
		for i in range(9):
			l=inf.readline()
		l=l.split()[2:]
		#print(l)
		for icol,item in enumerate(l):
			if item == 'type':
				types_col = icol
		#print(icol)
		for iatoms in range(natoms):
			l=inf.readline().split()
			#print(l)
			types.append(int(l[types_col]))	
	return np.array(types,dtype=np.int32)

def readlammpsdump(infile,nstep,iskip=0,every=1,natoms=None,pos=True,vel=True, box=True , force = False):
	"""
	read lammps dump text file
	Parameters
	---------
	infile : str 
		input file
	nstep : int
		total number of step to read
	iskip : int
		skip first iskip steps
	every : int
		save "every" read steps. The output will have nstep//every step lenght
	natoms : int
		number of atoms, if None it read from the file
	pos : bool
		read the positions
	vel : bool 
		read the velocities
	force :
		read the forces
	Returns
	-------
	traj : dict
		dictonary with the trajectory 
	"""
	if(natoms is None):
		natoms = read_natoms_dump(infile)
	pos_col,vel_col,force_col = read_cols_to_read(infile,pos=pos,vel=pos,force=force)
	traj ={}
	nstep_tosave = nstep//every
	if (pos):
		traj['p']=np.zeros((nstep_tosave,natoms,3))
	if (vel):
		traj['v']=np.zeros((nstep_tosave,natoms,3))
	if (force) :
		traj['f']=np.zeros((nstep_tosave,natoms,3))
	#traj['box']=np.zeros(nstep_tosave,3,3)
	traj['box']=np.zeros((nstep_tosave,6))
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
				traj['box'][itimestep,0] = float(l[0])
				traj['box'][itimestep,1] = float(l[1])
				l=inf.readline().split()
				traj['box'][itimestep,2] = float(l[0])
				traj['box'][itimestep,3] = float(l[1])
				l=inf.readline().split()
				traj['box'][itimestep,4] = float(l[0])
				traj['box'][itimestep,5] = float(l[1])
				l=inf.readline() # ITEM: ATOMS
				for iatom in range(natoms):
					l=[float(j) for j in inf.readline().split()]
					if(pos):
						traj['p'][itimestep,iatom,:] = [l[jj] for jj in pos_col]
					if(vel):
						traj['v'][itimestep,iatom,:] = [l[jj] for jj in vel_col]
					if (force) :
						traj['f'][itimestep,iatom,:] = [l[jj] for jj in force_col]
			else:
				for ii in np.arange(natoms+9):
					l=inf.readline()
				

	return traj
if __name__ == "__main__" : 
	import argparse
	parser = argparse.ArgumentParser(description = 'Convert a lammps dump trajectory file to a binary.',formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-p', '--prefix',
			type = str,
			required = True,
			help = 'filename.')
	parser.add_argument('--nstep',
			type = int, 
			default = None,
			help = 'Number of steps to convert.')
	parser.add_argument('--natoms',
			type = int, 
			default = None,
			help = 'Number of atoms.')
	parser.add_argument('--every',
			type = int, 
			default = 1,
			help = 'Save 1 step every --every steps.')
	parser.add_argument('--sskip',
			type = int, 
			default = 0,
			help = 'skip first sskip steps.')
	parser.add_argument('--wrap',
			action='store_true',
			default = False,
			help = 'wrap coordinates.')
	
	args = parser.parse_args()
	
	directory = args.directory
	prefix = args.prefix
	natoms = args.natoms
	nstep = args.nstep
	every = args.every
	skip = args.sskip
	
	import pyanalisi
	#print(types)
	print('Reading {}...'.format( prefix))
	data = readlammpsdump(infile='{}'.format( prefix),nstep=nstep,iskip=skip,every=every,natoms=natoms,pos=True,vel=True, box=True , force = False)
	print('Done.')
	analisi_traj = pyanalisi.Trajectory(data['p'], data['v'], data['types'], data['box'],False, args.wrap)
	print('Writing output files...')
	analisi_traj.write_lammps_binary('{}.bin'.format(  prefix), 0, -1)
	print('Done.')
