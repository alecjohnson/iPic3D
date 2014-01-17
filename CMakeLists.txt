cmake_minimum_required(VERSION 2.8.8) 

#
# Make sure that toolchain file is provided
#
if (NOT DEFINED CMAKE_TOOLCHAIN_FILE) 
  message(FATAL_ERROR "Toolchain file required. Try: <build_dir>$ cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/<toolchain_file>")
endif()

#
# Project declaration
#
project(iPic3D)

#
# External Dependencies
#
find_package(HDF5 COMPONENTS HL C REQUIRED)
# We use the HDF5 header files only, so unset the libraries
if (__MIC__)
  set(HDF5_LIBRARY_DIRS)
  set(HDF5_LIBRARIES)
  set(HDF5_HL_LIBRARIES)
endif()

#
# Set exec path
#
set(EXECUTABLE_OUTPUT_PATH exec/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib)

#
# include and lib directories
#

# include_directories: files there are accessible to the project
include_directories(
        include
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${HDF5_LIBRARY_DIRS}
        ${MPI_LIBRARY_DIRS}
	${EXTRA_LIBS}
)

#
# Header file list
#

file(
        GLOB
        inc_files
        include/*.h
)

#
# Source file list
#

# As there's no HDF5 library for MIC, we only link 
# empty functions for the symbols needed 
if (__MIC__)
  set(IPIC3D_MIC_HDF5_SYMBOLS mic/hdf5_symbols.c)
endif()

file(
        GLOB
        src_files
        ${IPIC3D_MIC_HDF5_SYMBOLS}
        ConfigFile/src/*.cpp
        PSKOutput3D/*.cpp
        bc/*.cpp
        communication/*.cpp
        fields/*.cpp
        grids/*.cpp
        inputoutput/*.cpp
        mathlib/*.cpp
        mpidata/*.cpp
        particles/*.cpp
        performances/*.cpp
        processtopology/*.cpp
        solvers/*.cpp
        utility/*.cpp
        main/iPic3Dlib.cpp
)

#
# Macro definitions
#

set(TEST_B $ENV{BATSRUS})
if(DEFINED TEST_B)
  add_definitions( -DBATSRUS )
  message("WARNING: BATSRUS flag is active.")
else(DEFINED TEST_B)
  message("INFO: BATSRUS is not active.")
endif(DEFINED TEST_B)

#
# Executable declaration
#

# Particle solver
add_executable(
        iPic3D
        iPic3D.cpp
)

# Fields solver
add_executable(
        iPic3D_fields
        iPic3D_fields.cpp
)


#build iPic as a library also
#libiPic3Dlib.so in folder lib
add_library(
        iPic3Dlib  #name of the library
        SHARED	   #type of the library
        ${inc_files}    # stuff to build the library
        ${src_files}
)

#
# Link external libraries
#

target_link_libraries(
         iPic3Dlib
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
	 ${EXTRA_LIBS}
)

target_link_libraries(
         iPic3D
         iPic3Dlib
)

target_link_libraries(
         iPic3D_fields
         iPic3Dlib
)


## to save the executable in the folder where the CMakeLists.txt file is, i.e. CMAKE_CURRENT_SOURCE_DIR
#set_target_properties(iPic3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

## debug releases have a _d appended to the executable
set_target_properties(iPic3D PROPERTIES DEBUG_POSTFIX "_d")

#
# Some useful infos
#
message("INFO: CMAKE_SYSTEM_PROCESSOR : ${CMAKE_SYSTEM_PROCESSOR}")
message("INFO: CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
message("INFO: CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")
message("INFO: MPI_CXX_COMPILER : ${MPI_CXX_COMPILER}")
if (__MIC__)
  message("INFO: CMAKE_C_COMPILER : ${CMAKE_C_COMPILER}")
  message("INFO: CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
  message("INFO: MPI_C_COMPILER : ${MPI_C_COMPILER}")
else()
  message("INFO: HDF5_ROOT : $ENV{HDF5_ROOT}")
endif()
