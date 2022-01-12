#!/usr/bin/env python3
##### coding=utf-8
import numpy as np
import argparse
import warnings
from subprocess import check_output
import re

BOHR = 0.529177249    # Bohr constant in Angstrom
#TAU  = 4.8378e-5  # tau_PW constant in ps
TAU  = 0.5*4.8378e-5  # tau_CP constant in ps
HARTREE = 27.211386245988 #eV
eV=1.60217662 #10^-19 J
bar=1.0e-6/eV #eV/Ang^3

def read_file_cp_pos_vel(prefix, natoms, nstep=None,skip=0,every=1,tometal=True):
	"""
	read QE cp.x file 
	Parameters 
	----------
	prefix : str
		prefix of the files. ex prefix.pos prefix.vel prefix.cel
	natoms : int
		number of atoms
	nstep : int
		total number of step to read after initial skip
	skip : int
		skip initial 'skip' steps
	every : int
		save "every" read steps. The output will have nstep//every step length
	tometal : bool
		convert to metal
	
	Returns
	--------
	data : dict
		dictonary with the trajectory
	"""	
	data = {}
	nstep_total = nstep
	nstep = nstep//every
	print('nstep_total = ',nstep_total)
	print('nstep_tosave = ',nstep)
	#filethe.readline()  # skip first line
	
	filepos = open(prefix + '.pos', 'r')
	data['pos']  = np.zeros((nstep,natoms,3), dtype = np.float64)
	
	filevel = open(prefix + '.vel', 'r')
	data['vel']  = np.zeros((nstep,natoms,3), dtype = np.float64)
	
	filecel = open(prefix + '.cel', 'r')
	data['cell'] = np.zeros((nstep, 3,3), dtype = np.float64)
	
	istep = 0
	while(istep < skip):
		istep+=1
		linecel = filecel.readline()
		linecel = filecel.readline()
		linecel = filecel.readline()
		linecel = filecel.readline()
		#print(linethe)
		linepos = filepos.readline()
		linevel = filevel.readline()
		for iatom in range(natoms):
			linepos = filepos.readline()
			linevel = filevel.readline()

	istep=0
	#while (istep < nstep_total):
	for istep_total in range(nstep_total):
		iread = istep_total%every
		if(iread==0):
			#print(linethe)
			linepos = filepos.readline()
			linevel = filevel.readline()
			linecel = filecel.readline()
			#if (len(linethe)==0) or (len(linepos)==0) or (len(linevel)==0) or (len(linecel)==0):  # EOF
			if  (len(linepos)==0) or (len(linevel)==0) or (len(linecel)==0):  # EOF
						raise RuntimeError("End Of file")
	
	
			# lettura posizioni
			#values = np.array(linepos.split(), dtype = np.float)
			values_pos = linepos.split()
			#lettura velocity
			values_vel = linevel.split()
			#lettura forza
			if len(values_pos) and len(values_vel):
				if (int(values_pos[0]) != int(values_vel[0]) ):
					#print(data['step'][istep], int(values_vel[0]))
					raise RuntimeError("Different timesteps between files of velocity and positions \n line={}  \n {} != {}".format(istep,int(values_pos[0]),int(values_vel[0])))
				for iatom in range(natoms):
					linepos = filepos.readline()
					values = np.array(linepos.split())
					data['pos'][istep,iatom,:] = values[:]
					linevel = filevel.readline()
					values = np.array(linevel.split())
					data['vel'][istep,iatom,:] = values[:]
	    
			#lettura cella
			#values = np.array(linecel.split(), dtype=np.float64)
			values = linecel.split()
			#print values, data['step'][istep]
			if len(values):
				if ( int(values_pos[0]) != int(values[0]) ):
					raise RuntimeError("Different timesteps between files of cell and positions \n line={}  \n {} != {}".format(istep,int(values_pos[0]),int(values[0])))
				for i in range(3):
					values = np.array(filecel.readline().split())
					data['cell'][istep, i,0] = values[0]
					data['cell'][istep, i,1] = values[1]
					data['cell'][istep, i,2] = values[2]
	
			istep += 1
		else:
			linecel = filecel.readline()
			linecel = filecel.readline()
			linecel = filecel.readline()
			linecel = filecel.readline()
			linepos = filepos.readline()
			linevel = filevel.readline()
			for iatom in range(natoms):
				linepos = filepos.readline()
				linevel = filevel.readline()
	if(tometal):
		conv_pos = BOHR
		conv_vel = BOHR/TAU
		conv_energy = HARTREE
		data['pos'] *= conv_pos
		data['vel'] *= conv_vel
		data['cell'] *=conv_pos
	
	return data



###########################################################################################################################################################################
### Parser
if __name__ == "__main__" : 
	parser = argparse.ArgumentParser(description = 'Convert a cp.x trajectory file to a binary.',formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-p', '--prefix',
			type = str,
			required = True,
			help = 'Prefix of the filename.')
	parser.add_argument('-s', '--species',
			type = str,
			required = True,
			help = 'A single line file with the succession of the atom indexes. \n \
	for example a simulation of a water molecule with the atoms in the order HHO will have a file as: 0 0 1 ')
	parser.add_argument('--nstep',
			type = int, 
			default = None,
                        required = True,
			help = 'Number of steps to convert.')
	parser.add_argument('--natoms',
			type = int, 
			required = True,
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
	
	prefix = args.prefix
	natoms = args.natoms
	species = args.species
	nstep = args.nstep
	every = args.every
	skip = args.sskip
	
	import pyanalisi
	types = np.loadtxt(species,dtype=np.int32).reshape(-1)
	#print(types)
	print('Reading {}...'.format( prefix))
	data = read_file_cp_pos_vel('{}'.format( prefix),natoms, nstep = nstep, skip=skip , every=every)
	print('Done.')
	analisi_traj = pyanalisi.Trajectory(data['pos'], data['vel'], types, data['cell'],pyanalisi.BoxFormat.CellVectors, args.wrap)
	print('Writing output files...')
	analisi_traj.write_lammps_binary('{}.bin'.format( prefix), 0, -1)
	print('Done.')
