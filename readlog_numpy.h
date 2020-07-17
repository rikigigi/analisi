#ifndef READLOG_NUMPY_H
#define READLOG_NUMPY_H

#include <pybind11/pybind11.h>

template <class TFLOAT>
class ReadLog_numpy
{
public:
    ReadLog_numpy(pybind11::buffer data, std::vector<std::string> headers);
    TFLOAT * line(unsigned int index){return data_ + index*ncols;}
    unsigned int n_timestep(){return ntimestep;}
    std::pair<unsigned int,bool> get_index_of(std::string header) {
        for (unsigned int i=0;i<ncols;++i) {
            if (header==headers.at(i)) {
                return {i,true};
            }
        }
        return {0,false};
    }
private:
    unsigned int ntimestep, ncols;
    std::vector<std::string> headers;
    TFLOAT * data_;


};

#endif // READLOG_NUMPY_H
