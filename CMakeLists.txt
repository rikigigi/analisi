project(analisi)
include(CheckIncludeFileCXX)

cmake_minimum_required(VERSION 2.8)

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
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_library( BPO boost_program_options)
#fftw3
find_library (FFTW NAMES "fftw3")
find_library (FFTW_TH NAMES "fftw3_threads")
find_library (FFTW_OMP NAMES "fftw3_omp")

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
if (FFTW_OMP)
    set (FFTW_MT "${FFTW_OMP}" )
elseif (FFTW_TH)
    set (FFTW_MT "${FFTW_TH}" )
else()
    set (FFTW_MT "" )
endif()
target_link_libraries(${PROJECT_NAME} "${BPO}" "${FFTW}" "${FFTW_MT}")
add_definitions(-DANALISI -DDEBUG)
#add_definitions(-DANALISI)


check_include_file_CXX ("eigen3/Eigen/Dense" HAVEeigen3EigenDense)
check_include_file_CXX ("Eigen/Dense" HAVEEigenDense)



set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fopenmp -march=native  -fsanitize=address -fsanitize=leak -fsanitize=undefined")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fopenmp -march=native -fsanitize=address -fsanitize=leak -fsanitize=undefined")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fopenmp -march=native ")
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fopenmp -march=native")
#if (C_COMPILER_SUPPORTS_LTO AND CXX_COMPILER_SUPPORTS_LTO)
#    message(STATUS "Enabled -flto link time optimization option.")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
#endif()
set(CMAKE_CXX_FLAGS_DEBUG)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h )
