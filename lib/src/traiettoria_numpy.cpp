#include "traiettoria_numpy.h"
#include <fstream>
#include "lammps_struct.h"
#include "buffer_utils.h"

Traiettoria_numpy::Traiettoria_numpy(pybind11::buffer &&buffer_pos_, pybind11::buffer &&buffer_vel_, pybind11::buffer &&buffer_types_, pybind11::buffer &&buffer_box_, bool matrix_box, bool pbc_wrap) :
    buffer_pos{buffer_pos_},buffer_vel{buffer_vel_},buffer_types{buffer_types_},buffer_box{buffer_box_}, matrix_box{matrix_box}
{
    loaded_timesteps=0;
    wrap_pbc=pbc_wrap;
    pybind11::buffer_info info_pos{buffer_pos.request()};
    pybind11::buffer_info info_vel{buffer_vel.request()};

    if (info_pos.ndim != 3) {
        throw std::runtime_error("Wrong number of dimension of position array (must be 3)");
    }
    if (info_vel.ndim != 3) {
        throw std::runtime_error("Wrong number of dimension of velocities array (must be 3)");
    }
    if (info_pos.shape[2]!=3)
        throw std::runtime_error("Wrong number of cartesian components in the third dimension of positions array");
    if (info_vel.shape[2]!=3)
        throw std::runtime_error("Wrong number of cartesian components in the third dimension of velocities array");

    for (unsigned long i=0;i<info_pos.shape.size();++i) {
        if (info_pos.shape[i] != info_vel.shape[i])
            throw std::runtime_error("Shape of positions and velocities array is different");
    }
    //check types
    pybind11::buffer_info info_types{buffer_types.request()};
    if (info_types.ndim != 1 )
        throw std::runtime_error("Wrong number of dimension of types array (must be 1)");
    if (info_types.shape[0]!=info_pos.shape[1])
        throw std::runtime_error("Wrong size of the type array");

    //check boxes
    pybind11::buffer_info info_box{buffer_box.request()};
    if (matrix_box){
        if (info_box.ndim != 3)
            throw std::runtime_error("Wrong number of dimensions of box array (must be 3)");
        if (info_box.shape[0]!=info_pos.shape[0] || info_box.shape[1]!=3 || info_box.shape[2] != 3)
            throw std::runtime_error("Wrong shape of box array");
    } else {
        if (info_box.ndim != 2)
            throw std::runtime_error("Wrong number of dimensions of box array (must be 2)");
        if (info_box.shape[0]!=info_pos.shape[0] || info_box.shape[1]!=6 )
            throw std::runtime_error("Wrong shape of box array: must be (nsteps, 6)");

    }

    //check formats of stuff
    if (info_box.format != pybind11::format_descriptor<double>::format())
        throw std::runtime_error("Format of box array should be double");
    if (info_types.format != pybind11::format_descriptor<int>::format())
        throw std::runtime_error("Format of types array should be int");
    if (info_vel.format != pybind11::format_descriptor<double>::format())
        throw std::runtime_error("Format of velocities array should be double");
    if (info_pos.format != pybind11::format_descriptor<double>::format())
        throw std::runtime_error("Format of positions array should be double");

    //check unsupported stride
    if (!has_no_stride<double>(info_box))
        throw std::runtime_error("Unsupported stride in box array");
    if (!has_no_stride<int>(info_types))
        throw std::runtime_error("Unsupported stride in types array");
    if (!has_no_stride<double>(info_vel))
        throw std::runtime_error("Unsupported stride in vel array");
    if (!has_no_stride<double>(info_pos))
        throw std::runtime_error("Unsupported stride in pos array");

    // set positions and velocities and boxes array
    natoms=info_pos.shape[1];
    n_timesteps=info_pos.shape[0];
    if (wrap_pbc){
        buffer_posizioni=new double[info_pos.shape[0]*info_pos.shape[1]*info_pos.shape[2]];
        posizioni_allocated=true;
        std::cerr << "WARNING: pbc lazily implemented only for a non rotated orthogonal cell"<<std::endl;
        for (int i=0;i<n_timesteps;++i) {
            for (int iatom=0;iatom<natoms;++iatom) {
                for (int idim=0;idim<3;++idim) {
                    double l=static_cast<double*>(info_box.ptr)[i*9+idim*3+idim];
                    double x=static_cast<double*>(info_pos.ptr)[i*3*natoms+iatom*3+idim];
                    buffer_posizioni[i*3*natoms+iatom*3+idim]=x-round(x/l)*l;
                }
            }
        }
    }else {
        posizioni_allocated=false;
        buffer_posizioni=static_cast<double*>(info_pos.ptr);
    }
    buffer_velocita=static_cast<double*>(info_vel.ptr);

    if (!matrix_box){
        buffer_scatola=static_cast<double *>(info_box.ptr);
    } else {
        buffer_scatola = new double[info_box.shape[0]*6];
        for (int i=0;i<info_box.shape[0];++i){
            //check it is an orthogonal (and non rotated) cell
            double * cur_box=static_cast<double*>(info_box.ptr)+9*i;
            for (int l=0;l<3;++l){
                for (int m=l+1;m<3;++m) {
                    if (cur_box[l*3+m] != 0.0 || cur_box[m*3+l] != 0.0){
                        throw std::runtime_error("Non orthogonal cell (or rotated orthogonal one) not implemented");
                    }
                }
                buffer_scatola[i*6+2*l]=0.0;
                buffer_scatola[i*6+2*l+1]=cur_box[3*l+l];
            }
        }
    }
    //TODO: convert everything in the program to general lattice matrix format

    buffer_tipi=static_cast<int*>(info_types.ptr);
    buffer_tipi_id=new int[natoms];
    get_ntypes(); //init type ids

    //
    //calculate center of mass velocity and position (without pbc)
    calc_cm_pos_vel(static_cast<double*>(info_pos.ptr),buffer_posizioni_cm);
    calc_cm_pos_vel(static_cast<double*>(info_vel.ptr),buffer_velocita_cm);

    loaded_timesteps=n_timesteps;
}

void
Traiettoria_numpy::calc_cm_pos_vel(double * a, double * & cm){
    if (cm==nullptr){
        cm=new double[3*n_timesteps*ntypes];
    } else {
        return;
    }
    int *cont=new int[ntypes];
    for (int i=0;i<n_timesteps;++i) {
        for (int itype=0;itype<ntypes;++itype) {
            cont[itype]=0;
            for (int icoord=0;icoord<3;icoord++){
                cm[(i*ntypes+itype)*3+icoord]=0.0;
            }
        }
        for (int iatom=0;iatom<natoms;++iatom){
            int itype=get_type(iatom);
            cont[itype]++;
            for (int icoord=0;icoord<3;icoord++){
                cm[(i*ntypes+itype)*3+icoord]+=(a[(i*natoms+iatom)*3+icoord]-cm[(i*ntypes+itype)*3+icoord])/double(cont[itype]);
            }
        }
    }
    delete [] cont;
}
void
Traiettoria_numpy::dump_lammps_bin_traj(const std::string &fname, int start_ts, int stop_ts){
    if (start_ts<0 || start_ts>=n_timesteps){
        throw std::runtime_error("You must provide a starting timestep between 0 and the number of timesteps!");
    }
    if (stop_ts<=0)
        stop_ts=n_timesteps;
    std::ofstream out(fname,std::ofstream::binary);
    for (int t=start_ts;t<stop_ts;++t){
        Intestazione_timestep head;
        head.natoms=natoms;
        for (unsigned int i=0;i<6;++i)
            head.scatola[i]=scatola(t)[i];
        head.timestep=t;
        head.triclinic=false;
        head.condizioni_al_contorno[0]=0;
        head.condizioni_al_contorno[1]=0;
        head.condizioni_al_contorno[2]=0;
        head.condizioni_al_contorno[3]=0;
        head.condizioni_al_contorno[4]=0;
        head.condizioni_al_contorno[5]=0;
        head.dimensioni_riga_output=NDOUBLE_ATOMO;
        head.nchunk=1;
        int n_data=head.natoms*head.dimensioni_riga_output;
        //write timestep header
        out.write((char*) &head,sizeof(Intestazione_timestep));
        out.write((char*) &n_data,sizeof(int));

        for (unsigned int iatom=0;iatom<head.natoms;++iatom) {
            double data[NDOUBLE_ATOMO];
            data[0]=iatom;
            data[1]=get_type(iatom);
            for (int i=0;i<3;++i){
                data[2+i]=posizioni(t,iatom)[i];
                data[5+i]=velocita(t,iatom)[i];
            }
            static_assert (NDOUBLE_ATOMO==8, "You have to change the file writing (what do I have to write?) this if you change NDOUBLE_ATOMO" );
            out.write((char*) data,sizeof(double)*NDOUBLE_ATOMO);
        }
    }
}

Traiettoria_numpy::~Traiettoria_numpy() {
    if (matrix_box) {
        delete [] buffer_scatola;
    }
    if (posizioni_allocated)
        delete [] buffer_posizioni;
    delete [] buffer_tipi_id;
    delete [] masse;
    delete [] cariche;
    delete [] buffer_posizioni_cm;
    delete [] buffer_velocita_cm;

}
