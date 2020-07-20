#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "traiettoria.h"
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

namespace py = pybind11;

template <class T>
void define_atomic_traj(py::module & m, std::string typestr){
    py::class_<Gofrt<double,T> >(m,(std::string("Gofrt")+typestr).c_str(), py::buffer_protocol())
            .def(py::init<T*,double,double,unsigned int, unsigned int, unsigned int,unsigned int, bool>(),R"begend(
                                                                                                          calculates g(r) and, in general, g(r,t). Interesting dynamical property
                                                                                                          Parameters                                                                                                          ----------                                                                                                         Trajectory instance                                                                                                          rmin                                                                                                          rmax                                                                                                          nbin                                                                                                          maximum time lag                                                                                                          number of threads                                                                                                          time skip                                                                                            debug flag                                                                                                          )begend")
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
            .def(py::init<T*,double,double,unsigned int, unsigned int, unsigned int,unsigned int,unsigned int, bool>(),R"lol(
                 Parameters
                 ----------
                 Trajectory instance
                 rmin
                 rmax
                 nbin
                 maximum time lag
                 number of threads
                 time skip
                 buffer size
                 debug flag
)lol")
            .def("reset",&SHC::reset)
            .def("calculate", &SHC::calcola)
            .def_buffer([](SHC & m) -> py::buffer_info {
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


    using MSD=MSD<T>;
    py::class_<MSD>(m,(std::string("MeanSquareDisplacement")+typestr).c_str(),py::buffer_protocol())
            .def(py::init<T *, unsigned int, unsigned int, unsigned int , bool, bool,bool>())
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

PYBIND11_MODULE(pyanalisi,m) {
    py::class_<Traiettoria>(m,"Traj")
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
            .def("getNtimesteps",&Traiettoria::get_ntimesteps,R"begend(
                 returns estimated number of timesteps from the file size
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

    py::class_<ReadLog_numpy<double> >(m,"ReadLog")
            .def(py::init<py::buffer, std::vector<std::string> >() ,R"begend(
                 proxy class between python arrays and library

                 Parameters
                 ----------
                 string -> 2d array with data
                 string list -> header list
)begend")
            .def("getNtimesteps", &ReadLog_numpy<double>::n_timestep,R"begend(
                 return number of timesteps read
)begend");

    py::class_<CorrelatoreSpaziale>(m,"_CurrentCalc", py::buffer_protocol())
            .def(py::init<Traiettoria*,std::vector< std::array<double,3> >,double,unsigned int,unsigned int,bool>(),R"begend(
                 calculates  \dot \tilde e(k) / |k|  (can be useful to define some new current)
)begend")
            .def("reset",&CorrelatoreSpaziale::reset,R"begend()begend")
            .def("calculate",&CorrelatoreSpaziale::calcola,R"begend()begend")
            .def("printRes",[](CorrelatoreSpaziale & c){
                c.print(std::cout);
            },R"begend()begend")
            .def("__enter__",[](CorrelatoreSpaziale & c) -> CorrelatoreSpaziale & { return c;} )
            .def("__exit__",[](CorrelatoreSpaziale & c, py::object * exc_type, py::object * exc_value, py::object * traceback){})
            .def_buffer([](CorrelatoreSpaziale & c) -> py::buffer_info{
                return py::buffer_info(
                c.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                c.get_shape().size(),                                      /* Number of dimensions */
                c.get_shape(),                 /* Buffer dimensions */
                c.get_stride()
        );
            });


    py::class_<HeatC>(m,"_HeatCurrentCalc", py::buffer_protocol())
            .def(py::init<Traiettoria*,double,unsigned int,unsigned int>(),R"begend(
                 calculates  something ill defined
)begend")
            .def("reset",&HeatC::reset,R"begend()begend")
            .def("calculate",&HeatC::calcola,R"begend()begend")
            .def("__enter__",[](HeatC & c) -> HeatC & { return c;} )
            .def("__exit__",[](HeatC & c, py::object * exc_type, py::object * exc_value, py::object * traceback){})
            .def_buffer([](HeatC & c) -> py::buffer_info{
                return py::buffer_info(
                c.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                c.get_shape().size(),                                      /* Number of dimensions */
                c.get_shape(),                 /* Buffer dimensions */
                c.get_stride()
        );
            });


    py::class_<CenterDiff>(m,"_CenterDiff", py::buffer_protocol())
            .def(py::init<Traiettoria *,unsigned int,unsigned int, unsigned int,bool,bool>())
            .def("reset",&CenterDiff::reset)
            .def("getNumberOfExtraTimestepsNeeded", &CenterDiff::numeroTimestepsOltreFineBlocco)
            .def("calculate", & CenterDiff::calcola)
            .def("setStartingCenter",&CenterDiff::set_starting_center)
            .def_buffer([](CenterDiff & c) -> py::buffer_info {
                return py::buffer_info(
                            c.accesso_lista(),
                            sizeof(double),
                            py::format_descriptor<double>::format(),
                            c.get_shape().size(),
                            c.get_shape(),
                            c.get_stride()
                );
            })
            .def("__enter__",[](CenterDiff & c) -> CenterDiff & { return c;} )
            .def("__exit__",[](CenterDiff & c, py::object * exc_type, py::object * exc_value, py::object * traceback){});


    py::class_<CenterOfMassDiff>(m,"_CenterOfMassDiff", py::buffer_protocol())
            .def(py::init<Traiettoria *,unsigned int,unsigned int>())
            .def("reset",&CenterOfMassDiff::reset)
            .def("getNumberOfExtraTimestepsNeeded", &CenterOfMassDiff::numeroTimestepsOltreFineBlocco)
            .def("calculate", & CenterOfMassDiff::calcola)
            .def("setStartingCenter",&CenterOfMassDiff::set_starting_center)
            .def("setZeroThreshold",&CenterOfMassDiff::set_zero_threshold)
            .def("getResetError",&CenterOfMassDiff::get_reset_error)
            .def_buffer([](CenterOfMassDiff & c) -> py::buffer_info {
                return py::buffer_info(
                            c.accesso_lista(),
                            sizeof(double),
                            py::format_descriptor<double>::format(),
                            c.get_shape().size(),
                            c.get_shape(),
                            c.get_stride()
                );
            })
            .def("__enter__",[](CenterOfMassDiff & c) -> CenterOfMassDiff & { return c;} )
            .def("__exit__",[](CenterOfMassDiff & c, py::object * exc_type, py::object * exc_value, py::object * traceback){});


    py::class_<Traiettoria_numpy>(m,"Trajectory")
            .def(py::init<py::buffer,py::buffer,py::buffer,py::buffer,bool,bool>(),R"lol(
                Parameters
                ----------
                positions (double) python array (ntimesteps,natoms,3)
                velocities (double) python array (ntimesteps,natoms,3)
                types (int) python array (natoms)
                lattice vectors (double) (ntimestep,3,3) -- currently only diagonal matrices are supported
)lol")

            .def("getNtimesteps",&Traiettoria_numpy::get_ntimesteps,R"begend(
                 returns number of timesteps
)begend")
            .def("write_lammps_binary",&Traiettoria_numpy::dump_lammps_bin_traj,R"lol(
                 Parameters
                 ----------
                 file name (string)
                 starting timestep (int)
                 end timestep (int)  -- if < 0 it will dump all the trajectory
)lol");

    define_atomic_traj<Traiettoria>(m,"_lammps");
    define_atomic_traj<Traiettoria_numpy>(m,"");
    gk<ReadLog<> >(m,"_file");
    gk<ReadLog_numpy<double> >(m,"");
    m.def("info",[]() -> std::string{
        return _info_msg ;
    });
}
