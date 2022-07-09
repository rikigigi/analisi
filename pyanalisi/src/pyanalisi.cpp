#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <pybind11/numpy.h>
#include "traiettoria.h"
#include "spettrovibrazionale.h"
#include "readlog.h"
#include "correlatorespaziale.h"
#include "gofrt.h"
#include "greenkuboNcomponentionicfluid.h"
#include "heatc.h"
#include "centerdiff.h"
#include "centerofmassdiff.h"
#include "msd.h"
#include "traiettoria_numpy.h"
#include "readlog_numpy.h"
#include "sphericalcorrelations.h"
#include "config.h"
#include "atomicdensity.h"
#include "steinhardt.h"
#include "neighbour.h"

namespace py = pybind11;

template <class SOP, class T>
void steinDef(py::module & m, std::string suffix){
    py::class_< SOP >(m,(std::string("SteinhardtOrderParameterHistogram")+suffix).c_str(),py::buffer_protocol())
            .def(py::init<T*,typename SOP::Rminmax_t,unsigned int, unsigned int, std::vector<unsigned int>,
                 unsigned int,unsigned int, bool,typename SOP::NeighListSpec,bool>(),
R"lol(
Parameters
----------
Trajectory instance
rminmax list
number of histograms (cut rmax - rmin in this number of part a fill a different istogram for each part)
size of each dimension of the histogram (you want a big number here)
list of l values to use for building the histogram
number of threads
time skip
if true do the histogram, otherwise calculate a trajectory-like object
provide a list of [(max_neighbours, rmax^2, rmax_2^2),...] if you want to use the SANN algorithm and use ibin=0 always. rminmax list is ignored
if true compute the averaged steinhardt order parameter (an additional loop over neighbours is performed). Implemented only with SANN
)lol")
            .def("reset",&SOP::reset)
            .def("calculate", &SOP::calcola)
            .def_buffer([](SOP & m) -> py::buffer_info {
        /*
        std::cerr <<"shape ("<< m.get_shape().size() << "): ";
        for (auto & n: m.get_shape()) std::cerr << n << " ";
        std::cerr <<std::endl<< m.lunghezza()<<std::endl;
        std::cerr << "allocated memory from"<<m.accesso_lista() << " to " <<m.accesso_lista() + m.lunghezza()<<std::endl;
        std::cerr <<"strides ("<<  m.get_stride().size() << "): ";
        for (auto & n: m.get_stride()) std::cerr << n << " ";
        std::cerr << std::endl;*/
        return py::buffer_info(
                    m.accesso_lista(),
                    sizeof(double),
                    py::format_descriptor<double>::format(),
                    m.get_shape().size(),
                    m.get_shape(),
                    m.get_stride()
                    );

    });
}

template <class T>
void define_atomic_traj(py::module & m, std::string typestr){
    py::class_<Gofrt<double,T> >(m,(std::string("Gofrt")+typestr).c_str(), py::buffer_protocol())
            .def(py::init<T*,double,double,unsigned int, unsigned int, unsigned int,unsigned int, unsigned int , bool>(),R"begend(
                                                                                                          calculates g(r) and, in general, g(r,t). Interesting dynamical property
                                                                                                          Parameters                                                                                                          ----------                                                                                                         Trajectory instance                                                                                                          rmin                                                                                                          rmax                                                                                                          nbin                                                                                                          maximum time lag                                                                                                          number of threads                                                                                                          time skip                                                    every time                                         debug flag                                                                                                          )begend")
            .def("reset",& Gofrt<double,T>::reset,R"begend()begend")
            .def("getNumberOfExtraTimestepsNeeded",&Gofrt<double,T>::numeroTimestepsOltreFineBlocco ,R"begend()begend")
            .def("calculate", & Gofrt<double,T>::calcola,R"begend()begend")
            .def_buffer([](Gofrt<double,T> & g) -> py::buffer_info {
        return py::buffer_info(
                    g.accesso_lista(),                               /* Pointer to buffer */
                    sizeof(double),                          /* Size of one scalar */
                    py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                    g.get_shape().size(),                                      /* Number of dimensions */
                    g.get_shape(),                 /* Buffer dimensions */
                    g.get_stride());
    });


    using SHC = SphericalCorrelations<10,double,T>;
    py::class_< SHC >(m,(std::string("SphericalCorrelations")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T*,typename SHC::rminmax_t,size_t,size_t,size_t,size_t,size_t, bool,typename SHC::NeighListSpec>(),R"lol(
                 Parameters
                 ----------
                 Trajectory instance
                 rminmax list
                 nbin
                 maximum time lag
                 number of threads
                 time skip
                 buffer size
                 debug flag
                 provide a list of [(max_neighbours, rmax^2, rmax_2^2),...] if you want to use the SANN algorithm and use ibin=0 always
)lol")
            .def("reset",&SHC::reset)
            .def("calculate", &SHC::calcola)
            .def_buffer([](SHC & m) -> py::buffer_info {
        return py::buffer_info(
                    m.accesso_lista(),
                    sizeof(double),
                    py::format_descriptor<double>::format(),
                    m.get_shape().size(),
                    m.get_shape(),
                    m.get_stride()
                    );

    });

    steinDef< Steinhardt<6,double,T>,T > (m,typestr);
    steinDef< Steinhardt<8,double,T>,T > (m,std::string("_8")+typestr);
    steinDef< Steinhardt<10,double,T>,T > (m,std::string("_10")+typestr);

    using MSD=MSD<T>;
    py::class_<MSD>(m,(std::string("MeanSquareDisplacement")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T *, unsigned int, unsigned int, unsigned int , bool, bool,bool>(),
                 R"begend(
                 Parameters
                 ----------
                 Trajectory instance
                 time skip for the average computation
                 max time
                 number of threads
                 calculate center of mass MSD
                 calculate for each atomic specie the atomic msd in the center of mass reference system of that specie
                 debug flag
                 )begend")
            .def("reset",&MSD::reset)
            .def("getNumberOfExtraTimestepsNeeded",&MSD::numeroTimestepsOltreFineBlocco)
            .def("calculate", &MSD::calcola)
            .def_buffer([](MSD & m) -> py::buffer_info {
                return py::buffer_info(
                            m.accesso_lista(),
                            sizeof(double),
                            py::format_descriptor<double>::format(),
                            m.get_shape().size(),
                            m.get_shape(),
                            m.get_stride()
                            );
            })
            .def("__enter__",[](MSD & m) -> MSD & {return m;})
            .def("__exit__",[](MSD & m, py::object * exc_type, py::object * exc_value, py::object * traceback){})
            ;
    using AD = AtomicDensity<T,long>;
    py::class_<AD>(m,(std::string("PositionHistogram")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T*,std::array<ssize_t, 3>,unsigned, unsigned>())
            .def("reset",&AD::reset)
            .def("calculate",&AD::calcola)
            .def_buffer([](AD & m)->py::buffer_info {
        return py::buffer_info(
                    m.accesso_lista(),
                    sizeof(long),
                    py::format_descriptor<long>::format(),
                    m.get_shape().size(),
                    m.get_shape(),
                    m.get_stride()
                    );
    })
            ;
    using VDOS=SpettroVibrazionale<T>;
    py::class_<VDOS>(m,(std::string("VibrationSpectrum")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T *, bool>(),
                 R"begend(
                 Parameters
                 ----------
                 Trajectory instance
                 debug flag
                 )begend")
            .def("reset",&VDOS::reset)
            .def("calculate", &VDOS::calcola)
            .def("spectrum", &VDOS::spettro)
            .def_buffer([](VDOS & m) -> py::buffer_info {
                return py::buffer_info(
                            m.accesso_lista(),
                            sizeof(double),
                            py::format_descriptor<double>::format(),
                            m.get_shape().size(),
                            m.get_shape(),
                            m.get_stride()
                            );
            })
            .def("__enter__",[](VDOS & m) -> VDOS & {return m;})
            .def("__exit__",[](VDOS & m, py::object * exc_type, py::object * exc_value, py::object * traceback){})
            ;

    using NEIGH = Neighbours<T,double>;
    py::class_<NEIGH>(m,(std::string("Neighbours")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T*,typename NEIGH::ListSpec>())
            .def("calculate_neigh", &NEIGH::update_neigh)
            .def("get_sann", [](NEIGH & n, size_t iatom, size_t jtype) {
        auto sannit = n.get_sann_r(iatom,jtype);
        const typename NEIGH::TType4 * foo=sannit.begin();

        return pybind11::array_t<double>(
        {{(long)sannit.size(),4}}, //shape
            {4*sizeof(double),sizeof(double)},
            (double*) foo
         );})
            .def("get_sann_idx", [](NEIGH & n, size_t iatom, size_t jtype) {
        auto sannit = n.get_sann(iatom,jtype);
        const size_t * foo=sannit.begin();

        return pybind11::array_t<size_t>(
        {{(long)sannit.size()}}, //shape
            {sizeof (size_t)},
            foo
         );})
            .def("get_neigh", [](NEIGH & n, size_t iatom, size_t jtype) {
        auto sannit = n.get_neigh_r(iatom,jtype);
        const typename NEIGH::TType4 * foo=sannit.begin();

        return pybind11::array_t<double>(
        {{(long)sannit.size(),4}}, //shape
            {4*sizeof(double),sizeof(double)},
            (double*) foo
         );})
            .def("get_neigh_idx", [](NEIGH & n, size_t iatom, size_t jtype) {
        auto it = n.get_neigh(iatom,jtype);
        const size_t * foo=it.begin();

        return pybind11::array_t<size_t>(
        {{(long)it.size()}}, //shape
            {sizeof(size_t)},
            foo
         );})
                ;

}


template <class T>
void gk(py::module & m, std::string typestr){

    using GK = GreenKuboNComponentIonicFluid<T,double,double>;
    py::class_<GK>(m,("GreenKubo"+typestr).c_str(), py::buffer_protocol())
            .def(py::init<T*, std::string, unsigned int, std::vector<std::string>, bool, unsigned int, unsigned int, bool, unsigned int, unsigned int, bool, unsigned int, unsigned int>(),
                 R"begend(
                 Parameters
                 ----------
                 ReadLog -> log object
                 string -> dump file output name
                 unsigned int -> skip in doing the trajectory average
                 string list -> header list (first one is the heat current header). You have to specify the header of the first column of each triplet
                 bool -> if true dump a file
                 unsigned int -> max lenght of correlation function
                 unsigned int -> number of threads
                 bool -> if true subtract the average from the correlation function before integrating
                 unsigned int -> calculate the average starting from this time
                 unsigned int -> number of segments (it does impact only cache performance)
                 bool -> do a benchmark
                 unsigned int -> benchmark parameter start
                 unsigned int -> benchmark parameter stop
                 )begend")
            .def("getNumberOfExtraTimestepsNeeded",&GK::numeroTimestepsOltreFineBlocco,R"begend(
                perform an efficient, multithreaded Green-Kubo calculations of the transport coefficient of a condensed matter system, given the currents.
)begend")
            .def("reset",&GK::reset,R"begend()begend")
            .def("calculate",&GK::calcola,R"begend()begend")
            .def_buffer([](GK & g) -> py::buffer_info {
                return py::buffer_info(
                g.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                g.get_shape().size(),                                      /* Number of dimensions */
                g.get_shape(),                 /* Buffer dimensions */
                g.get_stride());
            });

}

template <class Tk,class T>
T& trajectory_common_interfaces(T&&t) {
    return t.def("write_lammps_binary",&Tk::dump_lammps_bin_traj,
R"lol(
    Parameters
    ----------
    file name (string)
    starting timestep (int)
    end timestep (int)  -- if < 0 it will dump all the trajectory
)lol")
     .def("get_positions_copy", [](Tk & t) {
        double * foo=nullptr;
        if (t.posizioni(t.get_current_timestep(),0) == nullptr) {
            return pybind11::array_t<double>();
        }
        long nts=t.get_nloaded_timesteps();
        long nat=t.get_natoms();
        foo = new double[nts*nat*3];
        std::memcpy(foo,t.posizioni(t.get_current_timestep(),0),sizeof (double)*nts*nat*3);
        pybind11::capsule free_when_done(foo, [](void *f) {
         double *foo = reinterpret_cast<double *>(f);
         std::cerr << "freeing memory @ " << f << "\n";
         delete[] foo;
        });
        return pybind11::array_t<double>(
        {{nts,nat,3}}, //shape
            {3*nat*sizeof (double),3*sizeof(double),sizeof(double)},
            foo,
            free_when_done
         );})
    .def("get_velocities_copy", [](Tk & t) {
        double * foo=t.velocita(t.get_current_timestep(),0);
        if (foo == nullptr) {
            return pybind11::array_t<double>();
        }
        long nts=t.get_nloaded_timesteps();
        long nat=t.get_natoms();
        foo = new double[nts*nat*3];
        std::memcpy(foo,t.velocita(t.get_current_timestep(),0),sizeof (double)*nts*nat*3);
        pybind11::capsule free_when_done(foo, [](void *f) {
         double *foo = reinterpret_cast<double *>(f);
         std::cerr << "freeing memory @ " << f << "\n";
         delete[] foo;
        });
        return pybind11::array_t<double>(
        {{nts,nat,3}}, //shape
            {3*nat*sizeof (double),3*sizeof(double),sizeof(double)},
            foo,
            free_when_done
         );})
    .def("get_box_copy", [](Tk & t) {
        double * foo=t.scatola(t.get_current_timestep());
        if (foo == nullptr) {
            return pybind11::array_t<double>();
        }
        long nts=t.get_nloaded_timesteps();
        long nb=t.get_box_stride();
        foo = new double[nts*nb];
        std::memcpy(foo,t.scatola(t.get_current_timestep()),sizeof (double)*nts*nb);
        pybind11::capsule free_when_done(foo, [](void *f) {
         double *foo = reinterpret_cast<double *>(f);
         std::cerr << "freeing memory @ " << f << "\n";
         delete[] foo;
        });
        return pybind11::array_t<double>(
        {{nts,nb}}, //shape
            {nb*sizeof (double),sizeof(double)},
            foo,
            free_when_done
         );})
     .def("get_nloaded_timesteps",&Tk::get_nloaded_timesteps)
     .def("getNtimesteps",&Tk::get_ntimesteps,R"begend(
     returns estimated number of timesteps from the file size
)begend")
     .def("get_current_timestep",&Tk::get_current_timestep,"return the first timestep currently loaded in this object (meaningful for the lammps binary trajectory interface)")
     .def("getWrapPbc",&Tk::get_pbc_wrap,"return the pbc wrapping of the trajectory around the center of the cell flag")
     .def("minImage",[](Tk & t,size_t i, size_t j, size_t it, size_t ij){
        double *x = new double[4];
        x[3]=t.d2_minImage(i,j,it,ij,x);
        pybind11::capsule free_when_done(x, [](void *f) {
           double *x = reinterpret_cast<double *>(f);
           std::cerr << "freeing memory @ " << f << "\n"; 
           delete [] x;
        });
        return pybind11::array_t<double>(
         {{4}}, //shape
         {{sizeof(double)}},
         x,
         free_when_done
         );
     });
}

PYBIND11_MODULE(pyanalisi,m) {
#ifdef BUILD_MMAP
    trajectory_common_interfaces<Traiettoria,py::class_<Traiettoria>>(
    py::class_<Traiettoria>(m,"Traj", py::buffer_protocol()))
            .def(py::init<std::string>(),R"begend(
                 Parameters
                 ----------
                 name of the binary file to open
)begend")
            .def("setWrapPbc",&Traiettoria::set_pbc_wrap,R"begend(
                 wrap all the atomic coordinate inside the simulation box, red from the binary file

                 Parameters
                 ----------
                 True/False
)begend")
            .def("setAccessWindowSize",[](Traiettoria & t,int ts) {return (int) t.imposta_dimensione_finestra_accesso(ts);},R"begend(
                 sets the size of the read block. Must fit in memory
                 
                 Parameters
                 ----------
                 int -> size of the block
                
                 Returns
                 ----------
                 1 -> success
                 0/-1 -> failure
)begend")
            .def("setAccessStart",[](Traiettoria & t,int ts) {return (int) t.imposta_inizio_accesso(ts);},R"begend(
                 sets the first timestep to read, and read the full block

                 Parameters
                 ----------
                 int -> first timestep to read
)begend")
    //.def("getPosition",[](Traiettoria & t, int ts,int natoms){return std::array<>})
            .def("serving_pos", &Traiettoria::serving_pos)
            .def("toggle_pos_vel", &Traiettoria::toggle_pos_vel)
            .def_buffer([](Traiettoria & g) -> py::buffer_info {

                return py::buffer_info(
                g.serving_pos()? g.posizioni_inizio() : g.velocita_inizio(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                g.get_shape().size(),                                      /* Number of dimensions */
                g.get_shape(),                 /* Buffer dimensions */
                g.get_stride());
            })
            .def("get_lammps_id", [](Traiettoria & t) {
                int * foo=t.get_lammps_id();
                pybind11::capsule free_when_done(foo, [](void *f) {
                 int *foo = reinterpret_cast<int *>(f);
                 std::cerr << "freeing memory @ " << f << "\n";
                 delete[] foo;
                });
                size_t natoms=t.get_natoms();
                return pybind11::array_t<int>(
                    {natoms}, //shape
                    {sizeof(int)},
                    foo,
                    free_when_done
                 );})
            .def("get_lammps_type", [](Traiettoria & t) {
                int * foo=t.get_lammps_type();
                pybind11::capsule free_when_done(foo, [](void *f) {
                 int *foo = reinterpret_cast<int *>(f);
                 std::cerr << "freeing memory @ " << f << "\n";
                 delete[] foo;
                });
                size_t natoms=t.get_natoms();
                return pybind11::array_t<int>(
                    {natoms}, //shape
                    {sizeof(int)},
                    foo,
                    free_when_done
                 );})

                
    ;

    using RL = ReadLog<double>;
    py::class_<RL>(m,"ReadLog_file")
            .def(py::init<std::string,Traiettoria*,unsigned int, unsigned int, unsigned int, std::vector<std::string> >() ,R"begend(
                 perform an advanced read of a column formatted textfile with headers. Support multithreaded reading of chunks of lines.

                 Parameters
                 ----------
                 string -> filename
                 Traj -> Traj instance, to eventually calculate some currents
                 unsigned int -> skip
                 unsigned int -> number of threads
                 unsigned int -> batch size
                 string list -> eventual list of commands to calculate additional currents using the Traj instance
)begend")
            .def("getNtimesteps", &RL::n_timestep,R"begend(
                 return number of timesteps read
)begend");

#endif
    py::class_<ReadLog_numpy<double> >(m,"ReadLog")
            .def(py::init<py::buffer, std::vector<std::string> >(), py::keep_alive<1, 2>(),R"begend(
                 proxy class between python arrays and library

                 Parameters
                 ----------
                 string -> 2d array with data
                 string list -> header list
)begend")
            .def("getNtimesteps", &ReadLog_numpy<double>::n_timestep,R"begend(
                 return number of timesteps read
)begend");

    trajectory_common_interfaces<Traiettoria_numpy,py::class_<Traiettoria_numpy>>(
    py::class_<Traiettoria_numpy>(m,"Trajectory"))
            .def(py::init<py::buffer,py::buffer,py::buffer,py::buffer,Traiettoria_numpy::BoxFormat,bool,bool>(),
	    py::keep_alive<1, 2>(),
	    py::keep_alive<1, 3>(),
	    py::keep_alive<1, 4>(),
	    py::keep_alive<1, 5>(),
	    R"lol(
                Parameters
                ----------
                positions (double) python array (ntimesteps,natoms,3)
                velocities (double) python array (ntimesteps,natoms,3)
                types (int) python array (natoms)
                lattice vectors (double) (ntimestep,3,3)
                format of lattice vectors (BoxFormat)
                wrap atoms inside the cell using pbc
                save rotation matrix if triclinic format is used and a rotation is needed
)lol")



            .def("get_rotation_matrix", [](Traiettoria_numpy & t) {
        double * foo=t.get_rotation_matrix(0);
        if (foo == nullptr) {
            return pybind11::array_t<double>();
        }
        long nts=t.get_ntimesteps();
        foo = new double[nts*9];
        std::memcpy(foo,t.get_rotation_matrix(0),sizeof (double)*nts*9);
        pybind11::capsule free_when_done(foo, [](void *f) {
         double *foo = reinterpret_cast<double *>(f);
         std::cerr << "freeing memory @ " << f << "\n";
         delete[] foo;
        });
        return pybind11::array_t<double>(
        {{nts,3,3}}, //shape
            {3*3*sizeof (double),3*sizeof(double),sizeof(double)},
            foo,
            free_when_done
         );})
    ;
    py::enum_<Traiettoria_numpy::BoxFormat>(m, "BoxFormat", py::arithmetic())
            .value("Invalid", Traiettoria_numpy::BoxFormat::Invalid)
            .value("CellVectors", Traiettoria_numpy::BoxFormat::Cell_vectors)
            .value("LammpsOrtho", Traiettoria_numpy::BoxFormat::Lammps_ortho)
            .value("LammpsTriclinic", Traiettoria_numpy::BoxFormat::Lammps_triclinic);

#ifdef BUILD_MMAP
    define_atomic_traj<Traiettoria>(m,"_lammps");
    gk<ReadLog<> >(m,"_file");
#endif
    define_atomic_traj<Traiettoria_numpy>(m,"");
    gk<ReadLog_numpy<double> >(m,"");
    m.def("info",[]() -> std::string{
        return _info_msg ;
    });
    m.def("has_mmap",[]() -> bool {
#ifdef BUILD_MMAP
		    return true;
#else
		    return false;
#endif
    });
}
