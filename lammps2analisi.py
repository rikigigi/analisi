# coding=utf-8

import numpy as np
import sys


def read_write_xyz(infile,outfile,cols_pos=[3,4,5],cols_vel=None):
    """
    """
   
    in_file = open(infile, "r") 
    out_file = open(outfile, "w")
    #out_file.write("This Text is going to out file\nLook at it and see\n")
    while True:
        line=in_file.readline()
        if not line: break
        line=in_file.readline()
#        out_file.write("ITEM: TIMESTEP\n")
#        out_file.write("{}\n".format(int(round(data['step'][itimestep]))))
        line=in_file.readline()
        natoms=in_file.readline()
#        out_file.write("ITEM: NUMBER OF ATOMS\n")
        out_file.write(natoms)
        line=in_file.readline()
#        out_file.write('ITEM: BOX BOUNDS pp pp pp\n')
        line=in_file.readline()
        out_file.write(line)
        line=in_file.readline()
        out_file.write(line)
        line=in_file.readline()
        out_file.write(line)
        line=in_file.readline()
#        out_file.write('ITEM: ATOMS id type x y z vx vy vz\n')
        for iatom in range(int(natoms)):
             line = in_file.readline()
             line_s=line.split()
             print line
             line_o=line_s[0]
             for col in cols_pos:
                 line_o += " {}".format(line_s[col])
             if cols_vel != None:
                 for col in cols_vel:
                     line_o+=" {}".format(line_s[col])
             else:
                 line_o+=" 0 0 0"
             out_file.write(line_o+'\n')
    out_file.close()
    return

if len(sys.argv)<3:
    print "Uso: {} [prefix_file_in] [file_out] [col_pos_x col_pos_y col_pos_z] [col_vel_x col_vel_y col_vel_z] ...".format(sys.argv[0])
    exit(-1)

col_pos=[3,4,5]
col_vel=None
if len(sys.argv) == 6 or len(sys.argv) == 9:
    col_pos=[int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])]
if len(sys.argv) == 9 :
    col_vel=[int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8])]

read_write_xyz(sys.argv[1],sys.argv[2],col_pos,col_vel)


