#ifndef BUFFER_UTILS_H
#define BUFFER_UTILS_H

#include "compiler_types.h"

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

#endif
