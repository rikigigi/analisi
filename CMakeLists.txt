project(analisi)
include(CheckIncludeFileCXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
cmake_minimum_required(VERSION 2.8)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "No build type specified: assuming RELEASE build.")
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_TESTS "Build test programs" ON)
option(PYTHON_INTERFACE "Build python interface" ON)
option(SYSTEM_FFTW3 "Use system's fftw3 library" OFF)
option(SYSTEM_EIGEN3 "Use system's eigen3 library" OFF)
option(SYSTEM_BOOST "Use system's boost library" OFF)
#option(SYSTEM_PYBIND11 "Use system's pybind11 library" OFF)
option(SYSTEM_XDRFILE "Use system's xdrfile" OFF)

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++1z" COMPILER_SUPPORTS_CXX0X)

#set(SAFE_CMAKE_REQUIRED_LINK_OPTIONS "${CMAKE_REQUIRED_LINK_OPTIONS}")
#set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=address")
#CHECK_CXX_COMPILER_FLAG("-fsanitize=address" COMPILER_HAS_SANITIZER_ADDR)
#set(CMAKE_REQUIRED_LINK_OPTIONS "${SAFE_CMAKE_REQUIRED_LINK_OPTIONS}")
#unset(SAFE_CMAKE_REQUIRED_LINK_OPTIONS)
CHECK_CXX_COMPILER_FLAG(" -fsanitize=address " COMPILER_HAS_SANITIZER_ADDR)
CHECK_CXX_COMPILER_FLAG(" -fsanitize=leak " COMPILER_HAS_SANITIZER_LEAK)
CHECK_CXX_COMPILER_FLAG(" -fsanitize=undefined " COMPILER_HAS_SANITIZER_UNDEF)
#CHECK_CXX_COMPILER_FLAG(" -ffpe-trap=invalid,overflow " COMPILER_HAS_FTRAP)
if (COMPILER_HAS_SANITIZER_ADDR)
    set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -fsanitize=address")
else()
    message(WARNING "Seems like compiler does not support -fsanitize=address, but I'm setting it ")
    set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -fsanitize=address -fno-omit-frame-pointer ")
endif()
if (COMPILER_HAS_SANITIZER_LEAK)
    set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -fsanitize=leak")
else()
    message(WARNING "Compiler does not support -fsanitize=leak")
endif()
#if (COMPILER_HAS_FTRAP)
#    set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -ffpe-trap=invalid,overflow")
#else()
#    message(WARNING "Compiler does not support -ffpe-trap=invalid,overflow")
#endif()
if (COMPILER_HAS_SANITIZER_UNDEF)
    set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -fsanitize=undefined")
else()
    message(WARNING "Compiler does not support -fsanitize=undefined")
endif()
#CHECK_CXX_COMPILER_FLAG("-flto" CXX_COMPILER_SUPPORTS_LTO)
#CHECK_C_COMPILER_FLAG("-flto" C_COMPILER_SUPPORTS_LTO)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif()



if (SYSTEM_BOOST)
find_library( BPO NAMES boost_program_options)

if (BPO)
    message (STATUS "program_options library found: ${BPO}")
else()
    message (FATAL_ERROR "boost_program_options library not found!")
endif()
else()
add_subdirectory(boost)
endif()

#fftw3
if(SYSTEM_FFTW3)
find_library (FFTW NAMES "fftw3")
find_library (FFTW_TH NAMES "fftw3_threads")
find_library (FFTW_OMP NAMES "fftw3_omp")
check_include_file_CXX ("fftw3.h" HAVEfftw3)
check_include_file_CXX ("fftw.h" HAVEfftw)
   if (  HAVEfftw3 OR HAVEfftw )
      message(STATUS "Using system fftw3 library")
   else()
      message(FATAL_ERROR "cannot find fftw3 library's headers!")
   endif()
else()
#or use included one
set(ENABLE_OPENMP OFF CACHE BOOL "Enables multithread in FFTW3")
set(ENABLE_THREADS ON CACHE BOOL "Enables multithread in FFTW3")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")

add_subdirectory(fftw-3.3.8)
include_directories( fftw-3.3.8/api )
set(HAVEfftw3 1)
set(FFTW fftw3)
set_property(TARGET ${FFTW} PROPERTY POSITION_INDEPENDENT_CODE ON)
set(FFTW_TH ${FFTW}_threads)
endif()

#xdrfile
if(SYSTEM_XDRFILE)
find_library (xdrfile NAMES "xdrfile" )
if (xdrfile)
    set (XDR_FILE "${xdrfile}")
    message(STATUS "xdrfile library found!")
    set (HAVE_XDRFILE 1)
else ()
    set (XDR_FILE "")
    message(STATUS "WARNING: xdrfile library not found (conversion from gromacs file disabled)")
endif()
else()
add_subdirectory(xdrfile-1.1.4)
set (XDR_FILE xdrfile)
    message(STATUS "using packaged xdrfile library")
    set (HAVE_XDRFILE 1)
endif()
set (ANALISI_LIB
    lib/src/calcoliblocchi.cpp
    lib/src/chargefluxts.cpp
    lib/src/convertibinario.cpp
    lib/src/convolution.cpp
    lib/src/correlatorespaziale.cpp
    lib/src/cronometro.C
    lib/src/doubleround.cpp
    lib/src/gofrt.cpp
    lib/src/greenkubo2componentionicfluid.cpp
    lib/src/greenkuboNcomponentionicfluid.cpp
    lib/src/heatfluxts.cpp
    lib/src/interpolator.cpp
    lib/src/istogrammaatomiraggio.cpp
    lib/src/istogrammavelocita.cpp
    lib/src/modivibrazionali.cpp
    lib/src/mp.cpp
    lib/src/msd.cpp
    lib/src/posizioniequilibrio.cpp
    lib/src/readlog.cpp
    lib/src/rnd.c
    lib/src/spettrovibrazionale.cpp
    lib/src/testtraiettoria.cpp
    lib/src/traiettoria.cpp
    lib/src/heatc.cpp
    lib/src/calcolamultithread.cpp
    lib/src/centerdiff.cpp
    lib/src/centerofmassdiff.cpp
    lib/src/specialfunctions.cpp
    lib/src/sphericalcorrelations.cpp
    lib/src/atomicdensity.cpp
    )


if(PYTHON_INTERFACE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/pybind11/tools")
#check python
# Add a CMake parameter for choosing a desired Python version
if(NOT PYBIND11_PYTHON_VERSION)
  set(PYBIND11_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling modules")
endif()
set(Python_ADDITIONAL_VERSIONS 3.9 3.8 3.7 3.6 3.5 3.4)
find_package(PythonLibsNew ${PYBIND11_PYTHON_VERSION} REQUIRED)
if (${PYTHONLIBS_FOUND})
    set (ANALISI_LIB
	"${ANALISI_LIB}"
	lib/src/traiettoria_numpy.cpp
	lib/src/readlog_numpy.cpp )
endif()
endif()

add_library(${PROJECT_NAME}_lib STATIC ${ANALISI_LIB})
target_link_libraries(${PROJECT_NAME}_lib rt)
target_include_directories(${PROJECT_NAME}_lib PUBLIC lib/include )
set_property(TARGET ${PROJECT_NAME}_lib PROPERTY POSITION_INDEPENDENT_CODE ON)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC Threads::Threads)
add_executable(${PROJECT_NAME} 
	analisi/main.cpp
	)
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)


## test suite
if (BUILD_TESTS)
add_library(test_fixtures INTERFACE)
target_include_directories(test_fixtures INTERFACE tests/src/)
add_executable(test_sh tests/src/test_sh.cpp)
if (SYSTEM_BOOST)
set(BOOST_LIB "")
else()
set(BOOST_LIB boost)
endif()
target_link_libraries(test_sh ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)
add_executable(test_2loop tests/src/test_double_loop_splitter.cpp)
target_link_libraries(test_2loop ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)
add_executable(test_calcBuffer tests/src/test_calc_buffer.cpp)
target_link_libraries(test_calcBuffer ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)
add_executable(test_gk tests/src/test_gk.cpp)
target_link_libraries(test_gk ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)
add_executable(test_ph tests/src/test_position_histogram.cpp)
target_link_libraries(test_ph ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)
add_executable(test_lammps2020 tests/src/test_lammps2020.cpp)
target_link_libraries(test_lammps2020 ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures)

enable_testing()
add_test(NAME sh COMMAND test_sh)
add_test(NAME double_loop_splitter COMMAND test_2loop)
add_test(NAME calcBuffer COMMAND test_calcBuffer)
add_test(NAME gk COMMAND test_gk)
add_test(NAME atomic_density COMMAND test_ph)
add_test(NAME msd COMMAND test_lammps2020)
endif()

if(PYTHON_INTERFACE)
if (${PYTHONLIBS_FOUND})
#python interface
set (PYTHON_SUPPORT ${PYTHON_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
include_directories(pybind11/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory( pybind11 )
pybind11_add_module(pyanalisi 
	pyanalisi/src/pyanalisi.cpp
	)
target_link_libraries(pyanalisi PRIVATE ${PROJECT_NAME}_lib)
endif()
endif()


if (FFTW_TH)
    set (FFTW_MT "${FFTW_TH}" )
elseif (FFTW_OMP)
    set (FFTW_MT "${FFTW_OMP}" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
else()
    set (FFTW_MT "" )
endif()


target_link_libraries(${PROJECT_NAME} boost_program_options "${FFTW}" "${FFTW_MT}" "${XDR_FILE}")
target_link_libraries(${PROJECT_NAME}_lib "${FFTW}" "${FFTW_MT}" "${XDR_FILE}")
add_definitions(-DANALISI -DDEBUG)
#add_definitions(-DANALISI)

# look for environment variables like $FFTW_INCLUDE and $BOOST_INCLUDE (found in cineca's cluster nodes) and add to the include path

if(DEFINED ENV{FFTW_INC})
    message(STATUS "Found FFTW_INC=$ENV{FFTW_INC}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I $ENV{FFTW_INC}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I $ENV{FFTW_INC}")
endif()
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

if(SYSTEM_EIGEN3)
check_include_file_CXX ("eigen3/Eigen/Dense" HAVEeigen3EigenDense)
check_include_file_CXX ("Eigen/Dense" HAVEEigenDense)
   if (  HAVEeigen3EigenDense OR HAVEEigenDense )
      message(STATUS "Using system eigen3 library")
   else()
      message(FATAL_ERROR "cannot find eigen3 library!")
   endif()
else()
set (EIGEN_DOWNLOADED_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen3")
   include_directories(${EIGEN_DOWNLOADED_LIB_DIR})
   set(HAVEEigenDense 1)
   message(STATUS "Using packaged eigen3 library")
endif()




set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DEIGEN_INITIALIZE_MATRICES_BY_NAN ${COMPILER_SANITIZE_FLAGS}")
set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${COMPILER_SANITIZE_FLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
if (C_COMPILER_SUPPORTS_LTO AND CXX_COMPILER_SUPPORTS_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "Enabled -flto link time optimization option in RELEASE code.")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto")
endif()

set(PROJ_DIR "${CMAKE_CURRENT_LIST_DIR}")

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h )

include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

message(STATUS "python install library = ${PYTHON_SITE_PACKAGES}" )
install(TARGETS pyanalisi DESTINATION ${PYTHON_SITE_PACKAGES} COMPONENT python )

