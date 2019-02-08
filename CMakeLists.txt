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
    set (HAVE_XDRFILE 1)
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

set (EIGEN_DOWNLOADED_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-eigen-323c052e1731")

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/eigen-eigen-323c052e1731)
	include_directories(${EIGEN_DOWNLOADED_LIB_DIR})
    set(HAVEEigenDense 1)
    message(STATUS "Using downloaded eigen3 library (if you have problems delete directory ${EIGEN_DOWNLOADED_LIB_DIR})")
else()

check_include_file_CXX ("eigen3/Eigen/Dense" HAVEeigen3EigenDense)
check_include_file_CXX ("Eigen/Dense" HAVEEigenDense)


if (  HAVEeigen3EigenDense OR HAVEEigenDense )
message(STATUS "Using system eigen3 library")
else()
    message(STATUS "Downloading eigen3 library and unpacking in source directory")
    file (DOWNLOAD
        http://bitbucket.org/eigen/eigen/get/3.3.7.tar.bz2
        ${CMAKE_CURRENT_SOURCE_DIR}/eigen3.3.7.tar.bz2
        SHA512=c12bfd034d0a1112bf2df3f773ad98f36b2b53fbbfceac5f143ee1793322746fb8a6546d7db7da2a1000da6a535bd8cea0c4125b549afc90a6570743b02dcf3d
        SHOW_PROGRESS
       )
    execute_process(
	COMMAND ${CMAKE_COMMAND} -E tar xvf ${CMAKE_CURRENT_SOURCE_DIR}/eigen3.3.7.tar.bz2
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)

endif()
endif()

check_include_file_CXX ("fftw3.h" HAVEfftw3)
check_include_file_CXX ("fftw.h" HAVEfftw)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DEIGEN_INITIALIZE_MATRICES_BY_NAN -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
if (C_COMPILER_SUPPORTS_LTO AND CXX_COMPILER_SUPPORTS_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Enabled -flto link time optimization option in RELEASE code.")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
endif()



CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )

include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

add_subdirectory(libcg)
