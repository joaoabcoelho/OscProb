# This version is relatively old, but it is the one shipped in the image used to
# compile KOFI. As soon as this requirement will not be needed anymore, we can
# increase the version so that this simplifies the code in dependencies.cmake
cmake_minimum_required(VERSION 3.22.1)
project(oscprob LANGUAGES CXX)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(OSCPROB_ENABLE_TESTING "Enable OscProb tests" OFF)

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

if(OSCPROB_ENABLE_TESTING)
  include(CTest)
endif()

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
          src/PMNS_Maltoni.cxx
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
          src/PMNS_Avg.cxx
          src/PremModel.cxx
          src/complexsolver.cxx)
target_include_directories(
  OscProb
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/inc>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:include>)

set(hdrs
    Absorption.h
    EarthModelBase.h
    EarthModelBinned.h
    EigenPoint.h
    NuPath.h
    PMNS_Base.h
    PMNS_Maltoni.h
    PMNS_Decay.h
    PMNS_Deco.h
    PMNS_DensityMatrix.h
    PMNS_Fast.h
    PMNS_Iter.h
    PMNS_LIV.h
    PMNS_NSI.h
    PMNS_NUNM.h
    PMNS_OQS.h
    PMNS_SNSI.h
    PMNS_Sterile.h
    PMNS_Avg.h
    PremModel.h
    Definitions.h)

set(hdrs_full_path ${hdrs})
list(TRANSFORM hdrs_full_path PREPEND "${PROJECT_SOURCE_DIR}/inc/")

list(APPEND hrds prem_default.hpp)
list(APPEND hdrs_full_path ${PROJECT_BINARY_DIR}/prem_default.hpp)

set_target_properties(OscProb PROPERTIES PUBLIC_HEADER "${hdrs_full_path}"
                                         OUTPUT_NAME OscProb)

# Link libraries
target_link_libraries(
  OscProb
  PRIVATE MatrixDecomp
  PUBLIC ROOT::Core Eigen3::Eigen)

add_root_dictionary(OscProb HEADERS ${hdrs} LINKDEF inc/LinkDef.h)

add_subdirectory(MatrixDecomp)

if(OSCPROB_ENABLE_TESTING)
  add_subdirectory(test)
  # file(CREATE_LINK ${PROJECT_SOURCE_DIR}/inc ${PROJECT_BINARY_DIR}/inc
  # SYMBOLIC) file(CREATE_LINK ${PROJECT_SOURCE_DIR}/PremTables
  # ${PROJECT_BINARY_DIR}/stage/PremTables SYMBOLIC) file(CREATE_LINK
  # ${PROJECT_SOURCE_DIR}/EarthTables ${PROJECT_BINARY_DIR}/stage/EarthTables
  # SYMBOLIC) file(COPY ${PROJECT_SOURCE_DIR}/tutorial DESTINATION
  # ${PROJECT_BINARY_DIR}/stage/)
  install(DIRECTORY tutorial DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()

install(
  TARGETS OscProb
  EXPORT oscprob
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(NOT EXISTS "${CMAKE_INSTALL_INCLUDEDIR}/Eigen/Core")
  # The Eigen version 3.4.0 has some issues with CMake. This ensures that the
  # headers are installed where expected
  if(eigen3_POPULATED)
    set(eigen_inc ${eigen3_SOURCE_DIR})
  else()
    get_target_property(eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
  endif()
  message(STATUS "Explicitly installing Eigen/Core")
  install(DIRECTORY "${eigen_inc}/Eigen"
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()

include(cmake/packaging/CMakeLists.txt)

if(NOT PROJECT_IS_TOP_LEVEL)
  # do not show those options in the basic options as they are probably of
  # little interest for project including us
  mark_as_advanced(OSCPROB_ENABLE_TESTING)
endif()
