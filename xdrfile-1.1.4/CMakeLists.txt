cmake_minimum_required(VERSION 3.9)
include(GNUInstallDirs)


set ( SRC src/trr2xtc.c  src/xdrfile.c  src/xdrfile_c_test.c  src/xdrfile_trr.c  src/xdrfile_xtc.c)

add_library(xdrfile STATIC ${SRC})
target_include_directories(xdrfile PUBLIC include/)
set_property(TARGET xdrfile PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(xdrfile PROPERTIES
           VERSION 1.1.4
)
