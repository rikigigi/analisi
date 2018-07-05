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
       print "nstep = ", nstep

    data = {}
    data['step']  = np.zeros(nstep, dtype=np.int64)
    data['time']  = np.zeros(nstep, dtype=np.float64)
    data['ekinc'] = np.zeros(nstep, dtype=np.float64)
    data['Tcell'] = np.zeros(nstep, dtype=np.float64)
    data['Tion']  = np.zeros(nstep, dtype=np.float64)
    data['econt'] = np.zeros(nstep, dtype=np.float64)
    data['pos']  = np.zeros((nstep,natoms,3), dtype=np.float64)
    data['vel']  = np.zeros((nstep,natoms,3), dtype=np.float64)
    data['cell'] = np.zeros((nstep,3), dtype=np.float64)

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
            for iatom in xrange(natoms):
                linepos = filepos.readline()
                values = np.array(linepos.split())
                data['pos'][istep,iatom,:] = values[:]
            
        #lettura velocità
        values = linevel.split()
        #print values,data[0][istep]
        if len(values):
            if (data['step'][istep] != int(values[0]) ):
                raise RuntimeError("Different timesteps between files of velocity and thermo")
            for iatom in xrange(natoms):
                linevel = filevel.readline()
                values = np.array(linevel.split())
                data['vel'][istep,iatom,:] = values[:]
        
        #lettura cella
        values = linecel.split()
        #print values, data['step'][istep]
        if len(values):
            if (data['step'][istep] != int(values[0]) ):
                raise RuntimeError("Different timesteps between files of cell and thermo")
            for i in xrange(3):
                values = np.array(filecel.readline().split())
                data['cell'][istep,i] = values[i]

        istep += 1
    return data

#def wrap_cubic_timestep(positions,l_cube,cellAxis,natoms):
#    position_extended=np.array
#    for iatom in range(natoms):
#        for icoord in range(3)
#        if positions[iatom,icoord]<0 or positions[iatom,icoord]>

def write_xyz(outfile, data, natoms_per_type, type_names=None):
    """
    Scrive un file nel formato leggibile dal programma per la successiva conversione nel formato binario.
    cp.x nell'output separa gli atomi per tipi. Questa funzione assume che la prima metà sono di tipo "1"
    e la seconda metà di tipo "2".
    outfile è il nome del file da scrivere.
    data è il risultato della chiamata a read_file_pos_vel
    l è la dimensione della cella cubica scritta nell'output. """
    
    out_file = open(outfile, "w")
    #out_file.write("This Text is going to out file\nLook at it and see\n")
    nsteps = data['pos'].shape[0]
    natoms = data['pos'].shape[1]
    if (natoms != sum(natoms_per_type)):
        raise ValueError('Sum of number of atoms per type does not match the total number of atoms.')
    if type_names is None:
        type_names = map(str, np.arange(1, len(natoms_per_type)+1))
    else:
        if (len(natoms_per_type) != len(type_names)):
            raise ValueError('Number of type_names not compatible with natoms_per_type.')
    for itimestep in xrange(nsteps):
#        out_file.write("ITEM: TIMESTEP\n")
#        out_file.write("{}\n".format(int(round(data['step'][itimestep]))))
#        out_file.write("ITEM: NUMBER OF ATOMS\n")
        out_file.write("{}\n".format(natoms))
#        out_file.write('ITEM: BOX BOUNDS pp pp pp\n')
        out_file.write('{} {}\n'.format(0, data['cell'][itimestep,0] * BOHR))
        out_file.write('{} {}\n'.format(0, data['cell'][itimestep,1] * BOHR))
        out_file.write('{} {}\n'.format(0, data['cell'][itimestep,2] * BOHR))
#        out_file.write('ITEM: ATOMS id type x y z vx vy vz\n')
        cumnattype = np.cumsum(np.append(0,natoms_per_type))
        for attype, nattype in enumerate(natoms_per_type):
            firstat = cumnattype[attype]
            lastat  = cumnattype[attype+1]
            for i, idat in enumerate(xrange(firstat, lastat)):
               out_file.write('{} {} {} {} {} {} {} {}\n'.format(idat+1, type_names[attype], \
                                data['pos'][itimestep,idat,0]*BOHR,     data['pos'][itimestep,idat,1]*BOHR,     data['pos'][itimestep,idat,2]*BOHR, \
                                data['vel'][itimestep,idat,0]*BOHR/TAU, data['vel'][itimestep,idat,1]*BOHR/TAU, data['vel'][itimestep,idat,2]*BOHR/TAU))
#            np.savetxt(out_file, np.vstack((np.arange(firstat+1,lastat+1), [type_names[attype]]*nattype, \
#                                           data['pos'][itimestep,firstat:lastat,:].T * BOHR, \
#                                           data['vel'][itimestep,firstat:lastat,:].T * BOHR / TAU)).T, \
#                       fmt='%d %s %f %f %f %f %f %f')
    out_file.close()
    return

if len(sys.argv)<4:
    print "Uso: {} [prefix_file_in] [file_out] [natoms1] [natoms2] ...".format(sys.argv[0])
    exit(-1)

natoms=[]
for n in range(3,len(sys.argv)):
    natoms.append(int(sys.argv[n]))
natoms=np.array(natoms)
tot_natoms=np.sum(natoms)

data=read_file_pos_vel(sys.argv[1],tot_natoms)
write_xyz(sys.argv[2],data,natoms)


