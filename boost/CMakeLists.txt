cmake_minimum_required(VERSION 3.9)
include(GNUInstallDirs)


add_library(boost INTERFACE)
target_include_directories(boost INTERFACE headers/)

file (GLOB SRC program_options/src/*.cpp) 

add_library(boost_program_options SHARED ${SRC})
#target_link_library(boost_program_options boost)
set_target_properties(boost_program_options PROPERTIES
         VERSION 1.73.0
)
target_include_directories(boost_program_options PUBLIC headers/ )
