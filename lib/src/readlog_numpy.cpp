#include "readlog_numpy.h"
#include "buffer_utils.h"
template <class TFLOAT>
ReadLog_numpy<TFLOAT>::ReadLog_numpy(pybind11::buffer data, std::vector<std::string> headers) : headers{headers}
{

    pybind11::buffer_info info_data{data.request()};
    if (info_data.ndim != 2) {
        throw std::runtime_error("Wrong number of dimensions of the log array (must be 2)");
    }

    if (info_data.format != pybind11::format_descriptor<TFLOAT>::format())
        throw std::runtime_error("Format of log array is wrong (should be some kind of float, try to change it if it does not work)");

    if (!has_no_stride<TFLOAT>(info_data))
        throw std::runtime_error("Unsupported stride in log array");

    if (info_data.shape[1] != headers.size())
        throw std::runtime_error("Number of headers is not equal to the number of columns");

    ntimestep=info_data.shape[0];
    ncols=info_data.shape[1];
    data_=static_cast<TFLOAT *> (info_data.ptr);


}

template class ReadLog_numpy<double>;
template class ReadLog_numpy<long double>;
