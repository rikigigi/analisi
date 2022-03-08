project(analisi)
include(CheckIncludeFileCXX)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
cmake_minimum_required(VERSION 3.9)
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
option(BUILD_MMAP "Build command line interface and lammps binary interface. If OFF removes dependency over unix's mmap" ON)

if (NOT PYTHON_INTERFACE AND NOT BUILD_MMAP)
message(FATAL_ERROR "Cannot set both PYTHON_INTERFACE and BUILD_MMAP to OFF: in this way I'm building nothing but thin air")
endif()

if (NOT BUILD_MMAP AND BUILD_TESTS)
message(FATAL_ERROR "tests are implemented only if BUILD_MMAP is ON")
endif()

if ( USE_MPI ) 

find_package(MPI REQUIRED)
message(STATUS "MPI FOUND:")
message(STATUS "  Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

endif()

################################################################################
# Linux specific configuration
################################################################################
include(CheckCXXSourceCompiles)
set(LRT "")
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  # Try to detect if it is necessary to link against librt.
  # Note that glibc < 2.17 required librt to be linked to use clock_gettime()
  # and friends.
  set(CLOCK_GETTIME_REQUIRES_LIBRT_TEST_CODE
    "
    #include <time.h>
    int main() {
      timespec res;
      int result = clock_gettime(CLOCK_REALTIME, &res);
      return result == 0;
    }
    "
  )
  check_cxx_source_compiles(
    "${CLOCK_GETTIME_REQUIRES_LIBRT_TEST_CODE}"
    CLOCK_GETTIME_NO_REQUIRE_LIBRT
  )
  if (NOT CLOCK_GETTIME_NO_REQUIRE_LIBRT)
    # Try again with librt
    message(STATUS "Failed to link against clock_gettime(), trying with librt")
    set(CMAKE_REQUIRED_LIBRARIES_OLD "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} rt")
    check_cxx_source_compiles(
      "${CLOCK_GETTIME_REQUIRES_LIBRT_TEST_CODE}"
      CLOCK_GETTIME_REQUIRES_LIBRT
    )
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES_OLD}")
    if (CLOCK_GETTIME_REQUIRES_LIBRT)
      set(LRT "rt")
    else()
      message(FATAL_ERROR "Failed to link against clock_gettime()")
    endif()
  endif()
endif()


include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
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
    message(WARNING "Seems like compiler does not support -fsanitize=address")
    #set(COMPILER_SANITIZE_FLAGS " ${COMPILER_SANITIZE_FLAGS} -fsanitize=address -fno-omit-frame-pointer ")
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
set (MMAP_LIB
    lib/src/readlog.cpp
    lib/src/testtraiettoria.cpp
    lib/src/traiettoria.cpp
    lib/src/greenkubo2componentionicfluid.cpp
    lib/src/posizioniequilibrio.cpp
    lib/src/chargefluxts.cpp
    )
set (ANALISI_LIB
    lib/src/mp.cpp
    lib/src/calcoliblocchi.cpp
    lib/src/convertibinario.cpp
    lib/src/convolution.cpp
    lib/src/correlatorespaziale.cpp
    lib/src/cronometro.C
    lib/src/doubleround.cpp
    lib/src/gofrt.cpp
    lib/src/greenkuboNcomponentionicfluid.cpp
    lib/src/heatfluxts.cpp
    lib/src/interpolator.cpp
    lib/src/istogrammaatomiraggio.cpp
    lib/src/istogrammavelocita.cpp
    lib/src/modivibrazionali.cpp
    lib/src/msd.cpp
    lib/src/rnd.c
    lib/src/spettrovibrazionale.cpp
    lib/src/heatc.cpp
    lib/src/calcolamultithread.cpp
    lib/src/centerdiff.cpp
    lib/src/centerofmassdiff.cpp
    lib/src/specialfunctions.cpp
    lib/src/sphericalcorrelations.cpp
    lib/src/atomicdensity.cpp
    lib/src/traiettoriabase.cpp
    )

if (BUILD_MMAP)
set (ANALISI_LIB ${ANALISI_LIB} ${MMAP_LIB} )
endif()

if(PYTHON_INTERFACE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/pybind11/tools")
#check python
# Add a CMake parameter for choosing a desired Python version
if(NOT PYBIND11_PYTHON_VERSION)
  set(PYBIND11_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling modules")
endif()
#set(Python_ADDITIONAL_VERSIONS 3.9 3.8 3.7 3.6 3.5 3.4)
find_package(PythonLibsNew ${PYBIND11_PYTHON_VERSION} REQUIRED)
if (${PYTHONLIBS_FOUND})
    set (ANALISI_LIB
	"${ANALISI_LIB}"
	lib/src/traiettoria_numpy.cpp
	lib/src/readlog_numpy.cpp )
endif()
endif()

add_library(${PROJECT_NAME}_lib STATIC ${ANALISI_LIB})
target_include_directories(${PROJECT_NAME}_lib PUBLIC lib/include )
if ( USE_MPI ) 
   target_link_libraries(${PROJECT_NAME}_lib PUBLIC MPI::MPI_CXX)
endif()
set_property(TARGET ${PROJECT_NAME}_lib PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET ${PROJECT_NAME}_lib PROPERTY CXX_STANDARD 17)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC Threads::Threads ${LRT})
if (BUILD_MMAP)
   add_executable(${PROJECT_NAME} 
   	analisi/main.cpp
   	)
   target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib)
   set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
   if ( USE_MPI ) 
      target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
   endif()
endif()

## test suite
if (BUILD_TESTS)
add_library(test_fixtures INTERFACE)
target_include_directories(test_fixtures INTERFACE tests/src/)
if (SYSTEM_BOOST)
set(BOOST_LIB "")
else()
set(BOOST_LIB boost)
endif()

add_executable(all_tests tests/src/all_tests.cpp)
target_link_libraries(all_tests ${PROJECT_NAME}_lib ${BOOST_LIB} test_fixtures ${LRT} )

set_property(TARGET all_tests PROPERTY CXX_STANDARD 17)

enable_testing()
add_test(NAME all_tests COMMAND all_tests)
endif()

if(PYTHON_INTERFACE)
if (${PYTHONLIBS_FOUND})
#python interface
set (PYTHON_SUPPORT ${PYTHON_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} PRIVATE pybind11::module)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
include_directories(pybind11/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory( pybind11 )
pybind11_add_module(pyanalisi 
	pyanalisi/src/pyanalisi.cpp
	)
target_link_libraries(pyanalisi PRIVATE ${PROJECT_NAME}_lib ${LTR} )
set_property(TARGET pyanalisi PROPERTY CXX_STANDARD 17)
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

if (BUILD_MMAP)
target_link_libraries(${PROJECT_NAME} PUBLIC boost_program_options "${FFTW}" "${FFTW_MT}" "${XDR_FILE}" ${LTR})
endif()
target_link_libraries(${PROJECT_NAME}_lib PRIVATE "${FFTW}" "${FFTW_MT}" "${XDR_FILE}")
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

