cmake_minimum_required(VERSION 3.16)
project(oscprob LANGUAGES CXX)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

include(GNUInstallDirs)

# Instructions to properly handle RPATH
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()

#
# Define our RPATH settings, for both build and install phases
#

# use, i.e. do not skip, the full RPATH in the _build_ tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, do not use the install RPATH already (will only be used when
# actually installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add to the install RPATH the (automatically determined) parts of the RPATH
# that point to directories outside the build tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# specify libraries directory relative to binaries one.
file(RELATIVE_PATH relDir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
     ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

if(APPLE)
  set(basePoint @loader_path)
else()
  set(basePoint $ORIGIN)
endif()

set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir}
                        ${basePoint}/../../lib)

include(dependencies.cmake)

# useful functions to ease declaration of libraries and executables for this
# project
include(cmake/functions/CMakeLists.txt)

# Add library
add_library(OscProb)
# This is needed if we want to compile this project from a dependency with
# fetchContent
add_library(OscProb::OscProb ALIAS OscProb)

target_sources(
  OscProb
  PRIVATE src/Absorption.cxx
          src/EarthModelBase.cxx
          src/EarthModelBinned.cxx
          src/EigenPoint.cxx
          src/NuPath.cxx
          src/PMNS_Base.cxx
          src/PMNS_Decay.cxx
          src/PMNS_Deco.cxx
          src/PMNS_DensityMatrix.cxx
          src/PMNS_Fast.cxx
          src/PMNS_Iter.cxx
          src/PMNS_LIV.cxx
          src/PMNS_NSI.cxx
          src/PMNS_NUNM.cxx
          src/PMNS_OQS.cxx
          src/PMNS_SNSI.cxx
          src/PMNS_Sterile.cxx
          src/PremModel.cxx
          src/complexsolver.cxx
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         inc
         ${CMAKE_BINARY_DIR}
         FILES
         inc/Absorption.h
         inc/EarthModelBase.h
         inc/EarthModelBinned.h
         inc/EigenPoint.h
         inc/NuPath.h
         inc/PMNS_Base.h
         inc/PMNS_Decay.h
         inc/PMNS_Deco.h
         inc/PMNS_DensityMatrix.h
         inc/PMNS_Fast.h
         inc/PMNS_Iter.h
         inc/PMNS_LIV.h
         inc/PMNS_NSI.h
         inc/PMNS_NUNM.h
         inc/PMNS_OQS.h
         inc/PMNS_SNSI.h
         inc/PMNS_Sterile.h
         inc/PremModel.h
         inc/complexsolver.h
         inc/Definitions.h
         ${CMAKE_BINARY_DIR}/prem_default.hpp)

# Link libraries
target_link_libraries(
  OscProb
  PRIVATE MatrixDecomp
  PUBLIC ROOT::RIO Eigen3::Eigen)

set_target_properties(OscProb PROPERTIES OUTPUT_NAME OscProb)
get_target_property(hdrs OscProb HEADER_SET)
add_root_dictionary(OscProb HEADERS ${hdrs} LINKDEF inc/LinkDef.h)

install(
  TARGETS OscProb
  EXPORT oscprob
  FILE_SET HEADERS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(MatrixDecomp)

include(cmake/packaging/CMakeLists.txt)
