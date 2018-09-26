project(analisi)
include(CheckIncludeFileCXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
cmake_minimum_required(VERSION 2.8)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "No build type specified: assuming RELEASE build.")
endif(NOT CMAKE_BUILD_TYPE)

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
CHECK_CXX_COMPILER_FLAG("-flto" CXX_COMPILER_SUPPORTS_LTO)
CHECK_C_COMPILER_FLAG("-flto" C_COMPILER_SUPPORTS_LTO)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_library( BPO NAMES boost_program_options)

if (BPO)
    message (STATUS "boost_program_options library found.")
else()
    message (FATAL_ERROR "boost_program_options library not found!")
endif()

#fftw3
find_library (FFTW NAMES "fftw3")
find_library (FFTW_TH NAMES "fftw3_threads")
find_library (FFTW_OMP NAMES "fftw3_omp")

#xdrfile
find_library (XDRFILE NAMES "xdrfile" )

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
if (FFTW_OMP)
    set (FFTW_MT "${FFTW_OMP}" )
elseif (FFTW_TH)
    set (FFTW_MT "${FFTW_TH}" )
else()
    set (FFTW_MT "" )
endif()


if (XDRFILE)
    set (XDR_FILE "${XDRFILE}")
    message(STATUS "xdrfile library found!")
    set (HAVE_XDRFILE "")
else ()
    set (XDR_FILE "")
    message(STATUS "WARNING: xdrfile library not found (conversion from gromacs file disabled)")
endif()

target_link_libraries(${PROJECT_NAME} "${BPO}" "${FFTW}" "${FFTW_MT}" "${XDR_FILE}")
add_definitions(-DANALISI -DDEBUG)
#add_definitions(-DANALISI)

# look for environment variables like $FFTW_INCLUDE and $BOOST_INCLUDE (found in cineca's cluster nodes) and add to the include path

if(DEFINED ENV{FFTW_INCLUDE})
    message(STATUS "Found FFTW_INCLUDE=$ENV{FFTW_INCLUDE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I $ENV{FFTW_INCLUDE}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I $ENV{FFTW_INCLUDE}")
endif()

if(DEFINED ENV{BOOST_INCLUDE})
    message(STATUS "Found BOOST_INCLUDE=$ENV{BOOST_INCLUDE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I $ENV{BOOST_INCLUDE}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I $ENV{BOOST_INCLUDE}")
endif()

check_include_file_CXX ("eigen3/Eigen/Dense" HAVEeigen3EigenDense)
check_include_file_CXX ("Eigen/Dense" HAVEEigenDense)
check_include_file_CXX ("fftw3.h" HAVEfftw3)
check_include_file_CXX ("fftw.h" HAVEfftw)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
if (C_COMPILER_SUPPORTS_LTO AND CXX_COMPILER_SUPPORTS_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Enabled -flto link time optimization option in RELEASE code.")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
endif()



CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )
