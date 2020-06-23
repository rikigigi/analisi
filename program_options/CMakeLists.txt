cmake_minimum_required(VERSION 3.9)

file (GLOB SRC src/*.cpp) 

include(GNUInstallDirs)
add_library(boost_program_options SHARED ${SRC})
set_target_properties(boost_program_options PROPERTIES
         VERSION 1.73.0
)
target_include_directories(boost_program_options PUBLIC headers/ )
