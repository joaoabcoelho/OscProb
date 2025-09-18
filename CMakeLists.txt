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

# include(CTest)

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
          src/complexsolver.cxx)
target_include_directories(
  OscProb
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/inc>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:include>)
list(
  APPEND
  PUB_HDRS
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
  inc/Definitions.h
  ${PROJECT_BINARY_DIR}/prem_default.hpp)
set_target_properties(OscProb PROPERTIES PUBLIC_HEADER "${PUB_HDRS}" OUTPUT_NAME
                                                                     OscProb)

# Link libraries
target_link_libraries(
  OscProb
  PRIVATE MatrixDecomp
  PUBLIC ROOT::RIO Eigen3::Eigen)

add_root_dictionary(OscProb HEADERS ${PUB_HDRS} LINKDEF inc/LinkDef.h)

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

include(cmake/packaging/CMakeLists.txt)

if(NOT PROJECT_IS_TOP_LEVEL)
  # do not show those options in the basic options as they are probably of
  # little interest for project including us
  mark_as_advanced(OSCPROB_ENABLE_TESTING)
endif()
