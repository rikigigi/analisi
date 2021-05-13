#!/bin/env python
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

def read_file_pos_vel(prefix, natoms, nstep=None,skip=0,every=1,tometal=True):
	"""
	Legge i file di output di quantum espresso (cartella dove sono posizionati i vari restart)
	Per esempio, se il prefisso is KCl-512:
	namepos is KCl-512.pos
	namevel is KCl-512.vel
	nstep is il numero di timestep da leggere (NON is determinato automaticamente!)
	natoms is il numero di atomi nella simulazione (NON is determinato automaticamente!)
	
	ritorna una lista, una per ogni timestep, con il contenuto:
	    [timestep,tempo]       [ posizioni ]               [velocita']
	dove [timestep,tempo] is una coppia di numeri, posizioni is un array numpy con le posizioni,
	e velocita' is un array numpy
	"""	
	def file_length( filename ):
	  i = -1
	  blank = 0
	  with open(filename) as f:
	    for i, l in enumerate(f,1):
	      if len(l) == 1:
	        blank += 1
	      pass
	  return i - blank
	
	if nstep is None:
		nlines = int(check_output(["wc", "-l", prefix + '.evp']).decode("utf8").split()[0])
		nstep = nlines -1
		#nstep = file_length(prefix + '.evp') - 1
		print("nstep not set: using all the steps in .evp file: nstep = {}".format(nstep))

# ToDo: possibilita' di leggere i file dal fondo
#	if nstep < 0:
#		reverse = True
#		nstep = -nstep
#	else:
#		reverse = False
#
	data = {}
	nstep_total = nstep
	nstep = nstep//every
	print('nstep_total = ',nstep_total)
	print('nstep_tosave = ',nstep)
	filethe = open(prefix + '.evp', 'r')
	data['step']  = np.zeros(nstep, dtype = np.int64)
	data['time']  = np.zeros(nstep, dtype = np.float64)
	data['ekinc'] = np.zeros(nstep, dtype = np.float64)
	data['Tcell'] = np.zeros(nstep, dtype = np.float64)
	data['Tion']  = np.zeros(nstep, dtype = np.float64)
	data['econt'] = np.zeros(nstep, dtype = np.float64)
	data['epot'] = np.zeros(nstep, dtype = np.float64)
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
		linethe = filethe.readline()
		if (linethe.split()[0] == '#') :
			print("Comment found in {}.evp at line {}. Please check that this is correct.".format(prefix, istep+1))
			linethe = filethe.readline()
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
	for istep_total in range(nstep_total)
		iread = istep_total%every
		if(iread==0):
			linethe = filethe.readline()
			#print(linethe)
			linepos = filepos.readline()
			linevel = filevel.readline()
			if read_force: linefor = filefor.readline()
			linecel = filecel.readline()
			if (len(linethe)==0) or (len(linepos)==0) or (len(linevel)==0) or (len(linecel)==0):  # EOF
				if read_force:
					if (len(linefor)==0):
						raise RuntimeError("End Of file")
	
			# controllo per commenti 
			if (linethe.split()[0] == '#') :
				print("Comment found in {}.evp at line {}. Please check that this is correct.".format(prefix, istep+1))
				linethe = filethe.readline()
			# lettura thermo
			values = np.array(linethe.split(), dtype = np.float)
			if len(values):
				#print istep, values[0], len(data['step'])
				data['step'][istep]  = int(values[0])
				data['time'][istep]  = values[1]
				if istep == 1:
					deltat = data['time'][1]-data['time'][0] 
				data['ekinc'][istep] = values[2]
				data['Tcell'][istep] = values[3]
				data['Tion'][istep]  = values[4]
				data['epot'][istep]  = values[5]
				data['econt'][istep] = values[8]
			else:
				istep -= 1
	
			# lettura posizioni
			#values = np.array(linepos.split(), dtype = np.float)
			values = linepos.split()
			#lettura velocity
			values_vel = linevel.split()
			#lettura forza
			if read_force:
				values_for = linefor.split()
				#values = np.array(linefor.split(), dtype=np.float)
				#print values,data[0][istep]
				if len(values_for):
					if (data['step'][istep] != int(values[0]) ):
						print(data['step'][istep], int(values[0]))
						raise RuntimeError("Different timesteps between files of forces and thermo")
					for iatom in range(natoms):
						linefor = filefor.readline()
						values = np.array(linefor.split())
						data['for'][istep,iatom,:] = values[:]
			if len(values) and len(values_vel):
				if (data['step'][istep] != int(values[0]) ):
					print(data['step'][istep], int(values[0]))
					raise RuntimeError("Different timesteps between files of positions and thermo")
				if (data['step'][istep] != int(values_vel[0]) ):
					print(data['step'][istep], int(values_vel[0]))
					raise RuntimeError("Different timesteps between files of velocity and thermo")
				if read_force:
					values_for = linefor.split()
					if len(values_for):
						if (data['step'][istep] != int(values_for[0]) ):
							print(data['step'][istep], int(values_for[0]))
							raise RuntimeError("Different timesteps between files of forces and thermo")
				for iatom in range(natoms):
					linepos = filepos.readline()
					values = np.array(linepos.split())
					data['pos'][istep,iatom,:] = values[:]
					linevel = filevel.readline()
					values = np.array(linevel.split())
					data['vel'][istep,iatom,:] = values[:]
					if read_force:
						linefor = filefor.readline()
						values = np.array(linefor.split())
						data['for'][istep,iatom,:] = values[:]
	    
	    
			#lettura cella
			#values = np.array(linecel.split(), dtype=np.float64)
			values = linecel.split()
			#print values, data['step'][istep]
			if len(values):
				if (data['step'][istep] != int(values[0]) ):
					print(data['step'][istep], int(values[0]))
					raise RuntimeError("Different timesteps between files of cell and thermo")
				for i in range(3):
					values = np.array(filecel.readline().split())
					data['cell'][istep, 3*i,0] = values[0]
					data['cell'][istep, 3*i,1] = values[1]
					data['cell'][istep, 3*i,2] = values[2]
	
			istep += 1
		else:
			linethe = filethe.readline()
			if (linethe.split()[0] == '#') :
				print("Comment found in {}.evp at line {}. Please check that this is correct.".format(prefix, istep+1))
				linethe = filethe.readline()
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
	if(tometal):
		conv_pos = BOHR
		conv_vel = BOHR/TAU
		conv_energy = HARTREE
		data['pos'] *= conv_pos
		data['vel'] *= conv_vel
		data['epot'] *=conv_energy
		data['cell'] *=conv_pos
	
	return data



###########################################################################################################################################################################
### Parser

parser = argparse.ArgumentParser(description = 'Convert a cp.x trajectory file to the LAMMPS trajectory format. The units are Angstrom and Angstrom/picosecond.')
parser.add_argument('-d', '--directory',
		type = str,
		required = False,
		help = 'Directory with the .pos and .vel files.',
		default = './tmp')
parser.add_argument('-p', '--prefix',
		type = str,
		required = True,
		help = 'Prefix of the filename.')
parser.add_argument('-s', '--species',
		nargs = '*',
		type = str,
		required = True,
		help = 'Sequence of atomic species in the simulation (in the same order as in the ATOMIC_SPECIES card in the cp.x input).')
parser.add_argument('-n', '--natm',
		nargs = '*',
		type = int,
		required = True,
		help = 'Number of atoms per species (in the same order as in the ATOMIC_SPECIES card in the cp.x input).')
parser.add_argument('-c', '--charge',
		nargs = '*',
		type = float,
		required = False,
		help = 'Oxidation number per species (in the same order as in the ATOMIC_SPECIES card in the cp.x input).')
parser.add_argument('--nstep',
		type = int, 
		default = None,
		help = 'Number of steps to convert.')
parser.add_argument('--tskip',
		type = int, 
		default = 1,
		help = 'Write 1 every tskip steps.')
parser.add_argument('--sskip',
		type = int, 
		default = 0,
		help = 'skip first sskip steps.')
parser.add_argument('--xyz',
		help = 'Write the coordinates in a .xyz file.',
		action = 'store_true',
		required = False,
		default = False)
parser.add_argument('--vel',
		help = 'Write also the velocities in a .xyz file.',
		action = 'store_true',
		required = False,
		default = False)
parser.add_argument('--analisi',
		help = 'Output data in analisi format.',
		action = 'store_true',
		default = False)
parser.add_argument('--vcm',
		help = 'Write the per-species velocity of the centre of mass.',
		action = 'store_true',
		required = False,
		default = False)
parser.add_argument('--raw',
		help = 'Write the .raw files needed for deepMD.',
		action = 'store_true',
		required = False,
		default = False)
parser.add_argument('--shuffle',
		help = 'In writing the .raw files needed for deepMD, shuffle the steps.',
		action = 'store_true',
		required = False,
		default = False)

args = parser.parse_args()

directory = args.directory
prefix = args.prefix
species = args.species
natm = args.natm
nstep = args.nstep
xyz = args.xyz
vel = args.vel
analisi = args.analisi
charge = args.charge
tskip = args.tskip
vcm = args.vcm
raw = args.raw
shuffle = args.shuffle
skip = args.sskip

if shuffle:
	if not raw:
		raise ValueError('--shuffle is only for --raw output.')

if isinstance(species, list):
	if not isinstance(natm, list):
		raise ValueError('--natm should have the same dimension of --species!')
	else:
		if len(species) != len(natm):
			raise ValueError('--natm should have the same dimension of --species!')
if isinstance(natm, list):
	if not isinstance(species, list):
		raise ValueError('--natm should have the same dimension of --species!')
	else:
		if len(species) != len(natm):
			raise ValueError('--natm should have the same dimension of --species!')

if isinstance(natm, list):
	natm_tot = np.sum(natm)
else:
	natm_tot = natm

print('Reading {}/{}...'.format(directory, prefix))
leggi = read_file_pos_vel('{}/{}'.format(directory, prefix), natm_tot, nstep = nstep, skip=skip)
print('Done.')
print('Writing output files...')
scrivi = write_xyz('{}.lammpstrj'.format(prefix), leggi, natm, species, xyz = xyz, vel = vel, charge=charge, tskip=tskip, vcm=vcm, raw=raw, shuffle=shuffle)
print('Done.')
