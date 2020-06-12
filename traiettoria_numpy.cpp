#include "traiettoria_numpy.h"

template <typename T>
bool has_no_stride(const pybind11::buffer_info & buf) {
    ssize_t stride=sizeof (T);
    for (int i=buf.ndim-1; i>=0; --i) {
        if (buf.strides[i]!=stride)
            return false;
        stride*=buf.shape[i];
    }
    return true;
}

Traiettoria_numpy::Traiettoria_numpy(pybind11::buffer buffer_pos, pybind11::buffer buffer_vel, pybind11::buffer buffer_types, pybind11::buffer buffer_box, bool lammps_box, bool pbc_wrap) : lammps_box{lammps_box}
{
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

    for (int i=0;i<info_pos.shape.size();++i) {
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
    if (info_box.ndim != 3)
        throw std::runtime_error("Wrong number of dimensions of box array (must be 3)");
    if (info_box.shape[0]!=info_pos.shape[0] || info_box.shape[1]!=3 || info_box.shape[2] != 3)
        throw std::runtime_error("Wrong shape of box array");

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
    if (!has_no_stride<double>(info_types))
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
        buffer_posizioni=static_cast<double*>(info_pos.ptr);
    }
    buffer_velocita=static_cast<double*>(info_vel.ptr);

    if (!lammps_box){
        buffer_scatola=static_cast<double *>(info_box.ptr);
    } else {
        buffer_scatola = new double[info_box.shape[0]*6];
        for (int i=0;i<info_box.shape[0];++i){
            //check it is an orthogonal (and non rotated) cell
            double * cur_box=static_cast<double*>(info_box.ptr)+9*i;
            for (int l=0;l<3;++l){
                for (int m=l+1;m<3;++m) {
                    if (cur_box[l*3+m] != 0.0){
                        throw std::runtime_error("Non orthogonal cell (or rotated orthogonal one) not implemented");
                    }
                }
            }
        }
    }
    //TODO: convert everything in the program to general lattice matrix format

    buffer_tipi=static_cast<int*>(info_types.ptr);
    buffer_tipi_id=new int[natoms];
    get_ntypes(); //init type ids
}

Traiettoria_numpy::~Traiettoria_numpy() {
    if (lammps_box) {
        delete [] buffer_scatola;
    }
    if (posizioni_allocated)
        delete [] buffer_posizioni;
    delete [] buffer_tipi_id;
    delete [] masse;
    delete [] cariche;

}
