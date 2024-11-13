/**
 * function to write an array to a binary stream in the numpy array format
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <sstream>

namespace npy {
    void write_npy(std::ostream& out, const double * data, const size_t data_size, const std::vector<ssize_t>& shape, bool fortran_order = false) {
        // Step 1: Write the magic string
        const char magic_string[] = "\x93NUMPY";
        out.write(magic_string, 6);

        // Step 2: Write the major and minor version numbers
        uint8_t major_version = 1;
        uint8_t minor_version = 0;
        out.write(reinterpret_cast<const char*>(&major_version), 1);
        out.write(reinterpret_cast<const char*>(&minor_version), 1);

        // Step 3: Create the header dictionary
        std::ostringstream header_stream;
        header_stream << "{'descr': '<f8', 'fortran_order': " << (fortran_order ? "True" : "False") << ", 'shape': (";
        for (size_t i = 0; i < shape.size(); ++i) {
            header_stream << shape[i];
            if (i < shape.size() - 1) {
                header_stream << ", ";
            }
        }
        header_stream << "), }";

        std::string header = header_stream.str();
        size_t header_len = header.size();
        size_t padding_len = 64 - ((10 + header_len) % 64);
        header.append(padding_len, ' ');

        // Step 4: Write the header length
        uint16_t header_len_le = static_cast<uint16_t>(header_len + padding_len);
        out.write(reinterpret_cast<const char*>(&header_len_le), 2);

        // Step 5: Write the header data
        out.write(header.c_str(), header.size());

        // Step 6: Write the array data
        out.write(reinterpret_cast<const char*>(data), data_size * sizeof(double));
    }
}