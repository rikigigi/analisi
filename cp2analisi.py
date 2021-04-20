# coding=utf-8

import numpy as np
import sys

BOHR = 0.529177249    # Bohr constant in Angstrom
TAU  = 0.5*4.8378e-5  # tau_CP constant in ps

def read_file_pos_vel(prefix, natoms, nstep=None):
    """
    Legge i file di output di quantum espresso (cartella dove sono posizionati i vari restart)
    Per esempio, se il prefisso è KCl-512:
    namepos è KCl-512.pos
    namevel è KCl-512.vel
    nstep è il numero di timestep da leggere (NON è determinato automaticamente!)
    natoms è il numero di atomi nella simulazione (NON è determinato automaticamente!)
    
    ritorna una lista, una per ogni timestep, con il contenuto:
        [timestep,tempo]       [ posizioni ]               [velocità]
    dove [timestep,tempo] è una coppia di numeri, posizioni è un array numpy con le posizioni,
    e velocità è un array numpy
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
       nstep = file_length(prefix + '.evp') - 1
       print("nstep = ", nstep)

    data = {}
    data['step']  = np.zeros(nstep, dtype=np.int64)
    data['time']  = np.zeros(nstep, dtype=np.float64)
    data['ekinc'] = np.zeros(nstep, dtype=np.float64)
    data['Tcell'] = np.zeros(nstep, dtype=np.float64)
    data['Tion']  = np.zeros(nstep, dtype=np.float64)
    data['econt'] = np.zeros(nstep, dtype=np.float64)
    data['pos']  = np.zeros((nstep,natoms,3), dtype=np.float64)
    data['vel']  = np.zeros((nstep,natoms,3), dtype=np.float64)
    data['cell'] = np.zeros((nstep,6), dtype=np.float64)

    filethe = open(prefix + '.evp')
    filethe.readline()  # skip first line
    filepos = open(prefix + '.pos')
    filevel = open(prefix + '.vel')
    filecel = open(prefix + '.cel')
    istep = 0
    while (istep < nstep):
        linethe = filethe.readline()
        linepos = filepos.readline()
        linevel = filevel.readline()
        linecel = filecel.readline()
        if (len(linethe)==0) or (len(linepos)==0) or (len(linevel)==0) or (len(linecel)==0):  # EOF
            raise RuntimeError("End Of file")

        # lettura thermo
        values = np.array(linethe.split())
        if len(values):
          #print istep, values[0], len(data['step'])
          data['step'][istep]  = values[0]
          data['time'][istep]  = values[1]
          data['ekinc'][istep] = values[2]
          data['Tcell'][istep] = values[3]
          data['Tion'][istep]  = values[4]
          data['econt'][istep] = values[8]
        else:
            istep -= 1

        # lettura posizioni
        values = linepos.split()
        #print linepos
        #print values, data['step'][istep]
        if len(values):
            if (data['step'][istep] != int(values[0]) ):
                raise RuntimeError("Different timesteps between files of positions and thermo")
            for iatom in range(natoms):
                linepos = filepos.readline()
                values = np.array(linepos.split())
                data['pos'][istep,iatom,:] = values[:]
            
        #lettura velocità
        values = linevel.split()
        #print values,data[0][istep]
        if len(values):
            if (data['step'][istep] != int(values[0]) ):
                raise RuntimeError("Different timesteps between files of velocity and thermo")
            for iatom in range(natoms):
                linevel = filevel.readline()
                values = np.array(linevel.split())
                data['vel'][istep,iatom,:] = values[:]
        
        #lettura cella
        values = linecel.split()
        #print values, data['step'][istep]
        if len(values):
            if (data['step'][istep] != int(values[0]) ):
                raise RuntimeError("Different timesteps between files of cell and thermo")
            for i in range(3):
                values = np.array(filecel.readline().split())
                data['cell'][istep,2*i] = values[i]

        istep += 1
    return data

#def wrap_cubic_timestep(positions,l_cube,cellAxis,natoms):
#    position_extended=np.array
#    for iatom in range(natoms):
#        for icoord in range(3)
#        if positions[iatom,icoord]<0 or positions[iatom,icoord]>


if len(sys.argv)<4:
    print("Uso: {} [prefix_file_in] [file_out] [type atom 1] [type atom 2] ...".format(sys.argv[0]))
    exit(-1)

tatoms=[]
for n in range(3,len(sys.argv)):
    tatoms.append(int(sys.argv[n]))
tatoms=np.array(tatoms,dtype='int32')
tot_natoms=tatoms.shape[0]

import pyanalisi

data=read_file_pos_vel(sys.argv[1],tot_natoms)
analisi_traj = pyanalisi.Trajectory(data['pos'],
                                    data['vel'],
                                    tatoms,
                                    data['cell'],
                                    False,
                                    False)
analisi_traj.write_lammps_binary(sys.argv[2]+'.bin', 0,-1)

