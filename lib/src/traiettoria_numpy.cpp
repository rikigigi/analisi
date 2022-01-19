#include "traiettoria_numpy.h"
#include <fstream>
#include "lammps_struct.h"
#include "buffer_utils.h"
#include "triclinic.h"

Traiettoria_numpy::Traiettoria_numpy(pybind11::buffer &&buffer_pos_,
                                     pybind11::buffer &&buffer_vel_,
                                     pybind11::buffer &&buffer_types_,
                                     pybind11::buffer &&buffer_box_,
                                     TraiettoriaBase::BoxFormat matrix_box,
                                     bool pbc_wrap) :
    buffer_pos{buffer_pos_},buffer_vel{buffer_vel_},buffer_types{buffer_types_},buffer_box{buffer_box_}
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
    if (matrix_box == BoxFormat::Cell_vectors){
        if (info_box.ndim != 3)
            throw std::runtime_error("Wrong number of dimensions of box array (must be 3) for cell matrix format");
        if (info_box.shape[0]!=info_pos.shape[0] || info_box.shape[1]!=3 || info_box.shape[2] != 3)
            throw std::runtime_error("Wrong shape of box array: must be (nsteps, 3, 3) for matrix format");
    } else {
        if (info_box.ndim != 2)
            throw std::runtime_error("Wrong number of dimensions of box array (must be 2 for lammps ortho/triclinic cell format)");
        if (info_box.shape[0]!=info_pos.shape[0] ){
            throw std::runtime_error("Wrong shape of box array: first dimension must be nsteps");
        }
        if (matrix_box == BoxFormat::Lammps_ortho && info_box.shape[1]!=6 ) {
            throw std::runtime_error("Wrong shape of box array: must be (:, 6) for orthogonal lammps cell format");
        }
        if (matrix_box == BoxFormat::Lammps_triclinic && info_box.shape[1]!=9 ) {
            throw std::runtime_error("Wrong shape of box array: must be (:, 9) for triclinic lammps cell format");
        }
    }
    //check formats of stuff
    if (info_box.format != pybind11::format_descriptor<double>::format())
        throw std::runtime_error("Format of box array should be double");
    if (info_types.format != pybind11::format_descriptor<int>::format())
        throw std::runtime_error("Format of types array should be int ("+ pybind11::format_descriptor<int>::format()+") while it was "+info_types.format);
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
    buffer_scatola_stride = 6;
    triclinic=false;
    std::vector<std::pair<ssize_t,TriclinicLammpsCell<double> >> cells_qr; // first step, rotation matrix
    if (matrix_box== BoxFormat::Lammps_ortho || matrix_box==BoxFormat::Lammps_triclinic){ //nothing to do; maybe here we could do a copy, but maybe not
        if (matrix_box==BoxFormat::Lammps_triclinic) {
            buffer_scatola_stride = 9;
            triclinic=true;
        }
        buffer_scatola = new double[info_box.shape[0]*buffer_scatola_stride];
        box_allocated=true;
        //copy and do a permutation of everything
        for (ssize_t i=0;i<n_timesteps;++i) {
            std::memcpy(buffer_scatola+i*buffer_scatola_stride,&static_cast<double*>(info_box.ptr)[i*buffer_scatola_stride],buffer_scatola_stride*sizeof (double));
            lammps_to_internal(buffer_scatola+i*buffer_scatola_stride);
        }
        std::cerr << "Input format is lammps"<<std::endl;
    } else if (matrix_box==BoxFormat::Cell_vectors){
        //decide if we have to use triclinic or orthogonal cell; goes through all cells and convert/test to lammps format
        for (ssize_t i=0;i<info_box.shape[0];++i){
            double * cur_box=static_cast<double*>(info_box.ptr)+9*i;
            //check for triclinic
            auto lammps_cell=TriclinicLammpsCell(cur_box);
            if (!lammps_cell.isDiagonal()) {
                triclinic=true;
            }
            if (cells_qr.size()==0 || !lammps_cell.is_same_cell(cells_qr.back().second.getCell())) {
                cells_qr.push_back({i,std::move(lammps_cell)});
            }
        }
        ssize_t stride;
        if (triclinic){
            std::cerr << "Detected non orthorombic simulation cell. Using triclinic format"<<std::endl;
            velocities_allocated=true;
            buffer_velocita=new double[info_pos.shape[0]*info_pos.shape[1]*info_pos.shape[2]];
            stride=9;
        } else {
            stride=6;
        }
        buffer_scatola = new double[info_box.shape[0]*stride];
        ssize_t cur_idx=0;
        for (ssize_t i=0;i<n_timesteps;++i){
            //copy everything in the box array and rotate velocities and positions
            cells_qr[cur_idx].second.set_lammps_cell(buffer_scatola+i*stride,triclinic);
            if (cells_qr[cur_idx].first-1 == i) cur_idx++;
        }
    } else {
        throw std::runtime_error("Invalid input cell format");
    }

    if (!triclinic) { //no need for a copy of the velocities
        velocities_allocated=false;
        buffer_velocita=static_cast<double*>(info_vel.ptr);
        box_format=BoxFormat::Lammps_ortho;
    } else {
        buffer_scatola_stride = 9;
        box_format=BoxFormat::Lammps_triclinic;
    }
    if (wrap_pbc || triclinic){ // when I need also a copy of the positions
        buffer_posizioni=new double[info_pos.shape[0]*info_pos.shape[1]*info_pos.shape[2]];
        posizioni_allocated=true;
    }else {
        posizioni_allocated=false;
        buffer_posizioni=static_cast<double*>(info_pos.ptr);
        velocities_allocated=false;
        buffer_velocita=static_cast<double*>(info_vel.ptr);
    }
    //now everything is allocated/moved. Do the work of translation to the lammps (wapped) format

    if (triclinic &&matrix_box==BoxFormat::Cell_vectors ) { // in this case we have to rotate velocities and positions to be consistent with the new cell
        //do a rotation of vel and pos
        ssize_t cur_idx=0;
        for (ssize_t i=0;i<n_timesteps;++i){
            // rotate velocities and positions
            for (ssize_t n=0;n<natoms;++n){
                double *pos = buffer_posizioni+i*3*natoms+n*3;
                double *vel = buffer_velocita+i*3*natoms+n*3;
                std::memcpy(pos,&static_cast<double*>(info_pos.ptr)[i*3*natoms+n*3],3*sizeof (double));
                std::memcpy(vel,&static_cast<double*>(info_vel.ptr)[i*3*natoms+n*3],3*sizeof (double));
                cells_qr[cur_idx].second.rotate_vec(pos);
                cells_qr[cur_idx].second.rotate_vec(vel);
            }
            if (cells_qr[cur_idx].first-1 == i) cur_idx++;
        }
    }
    buffer_tipi=static_cast<int*>(info_types.ptr);
    buffer_tipi_id=new int[natoms];
    get_ntypes(); //init type ids

    //
    //calculate center of mass velocity and position (without pbc)
    calc_cm_pos_vel(static_cast<double*>(info_pos.ptr),buffer_posizioni_cm);
    calc_cm_pos_vel(static_cast<double*>(info_vel.ptr),buffer_velocita_cm);


    if (wrap_pbc) {
        for (ssize_t i=0;i<n_timesteps;++i){
            if (triclinic) {
                pbc_wrap<true>(i);
            } else {
                pbc_wrap<false>(i);
            }
        }
    }



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
        Intestazione_timestep_triclinic head;
        head.natoms=natoms;
        for (unsigned int i=0;i<6;++i)
            head.scatola[i]=scatola(t)[i];
        internal_to_lammps(head.scatola);
        head.timestep=t;
        head.triclinic=triclinic;
        if (triclinic) {
            for (unsigned int i=0;i<3;++i)
                head.xy_xz_yz[i]=scatola(t)[6+i];
        }
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
        head.write(out);
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
    if (box_allocated)
        delete [] buffer_scatola;
    if (velocities_allocated)
        delete [] buffer_velocita;
    if (posizioni_allocated)
        delete [] buffer_posizioni;
    delete [] buffer_tipi_id;
    delete [] masse;
    delete [] cariche;
    delete [] buffer_posizioni_cm;
    delete [] buffer_velocita_cm;

}
