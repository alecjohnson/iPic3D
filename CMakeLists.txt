cmake_minimum_required(VERSION 2.8.8)

#
# Set project name
#

project(iPic3D)

#
# Compilation options
#
option(IPIC_HDF5_OUTPUT "write output" ON)
if(NOT IPIC_HDF5_OUTPUT)
  add_definitions(-DNO_HDF5)
endif()

#
# Find packages (already managed by CMake)
#

if(IPIC_HDF5_OUTPUT)
find_package(HDF5 COMPONENTS HL C REQUIRED)
endif()
find_package(MPI REQUIRED)

#
# Options
#

option(IPIC_TESTS "Set up the code tests" OFF)

#
# Include sub-directories
#

add_subdirectory(main)
## if(IPIC_H5HUT_OUTPUT)
##   add_subdirectory(H5hut-io)
##   set(H5HUTIO_LIBRARIES H5hut-io)
## endif()

if(IPIC_HDF5_OUTPUT)
include_directories(
    include
    ${HDF5_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
)
else()
include_directories(
    include
    ${MPI_INCLUDE_PATH}
)
endif()

#
# Executable compilation
#

add_executable(
    iPic3D
    iPic3D.cpp
)

#
# Libraries linked
#

if(IPIC_HDF5_OUTPUT)
  target_link_libraries(
    iPic3D
    iPic3Dlib
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
  )
else()
target_link_libraries(
    iPic3D
    iPic3Dlib
    ${MPI_LIBRARIES}
)
endif()

#
# Code testing
#

if(IPIC_TESTS)
  enable_testing()
  set(IPIC_TESTS_DIR "${CMAKE_BINARY_DIR}/tests" CACHE STRING "Location of the source files for iPic3D")

  add_test(NAME GEM-test
           COMMAND ${CMAKE_COMMAND}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -DIPIC_SOURCE_DIR=${CMAKE_SOURCE_DIR}
           -DIPICEXEC=$<TARGET_FILE:iPic3D>
           -DMPIEXEC=${MPIEXEC}
           -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
           -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-GEM.txt)

  add_test(NAME uname-test
           COMMAND ${CMAKE_COMMAND}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-uname.txt)

endif(IPIC_TESTS)
