#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "traiettoria.h"
#include "correlatorespaziale.h"

namespace py = pybind11;

PYBIND11_MODULE(pyanalisi,m) {
    py::class_<Traiettoria>(m,"Traj")
            .def(py::init<std::string>())
            .def("setWrapPbc",&Traiettoria::set_pbc_wrap)
            .def("getNtimesteps",&Traiettoria::get_ntimesteps)
            .def("setAccessWindowSize",[](Traiettoria & t,int ts) {return (int) t.imposta_dimensione_finestra_accesso(ts);})
            .def("setAccessStart",[](Traiettoria & t,int ts) {return (int) t.imposta_inizio_accesso(ts);});

    py::class_<CorrelatoreSpaziale>(m,"CurrentCalc", py::buffer_protocol())
            .def(py::init<Traiettoria*,std::vector< std::array<double,3> >,double,unsigned int,unsigned int,bool>())
            .def("reset",&CorrelatoreSpaziale::reset)
            .def("calculate",&CorrelatoreSpaziale::calcola)
            .def("printRes",[](CorrelatoreSpaziale & c){
                c.print(std::cout);
            })
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
}
