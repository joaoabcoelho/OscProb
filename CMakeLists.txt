# This version is relatively old, but it is the one shipped in the image used to
# compile KOFI. As soon as this requirement will not be needed anymore, we can
# increase the version so that this simplifies the code in dependencies.cmake
cmake_minimum_required(VERSION 3.22.1)
project(oscprob LANGUAGES CXX)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(OSCPROB_ENABLE_TESTING "Enable OscProb tests" OFF)

# Default configurations for RPATH etc
include(cmake/defaults.cmake)

include(dependencies.cmake)

# useful functions to ease declaration of libraries and executables for this
# project
include(cmake/functions.cmake)

add_subdirectory(src)
add_subdirectory(MatrixDecomp)

install(DIRECTORY tutorial DESTINATION ${CMAKE_INSTALL_PREFIX})

if(OSCPROB_ENABLE_TESTING)
  include(CTest)
  add_subdirectory(test)
endif()

# Add target to simplify test of installed software We add it here since it is
# not really a unit test. We simply want to mimic the behaviour of make test
add_custom_target(
  test-installed
  DEPENDS OscProb
  COMMAND ${CMAKE_COMMAND} -E env ${ROOT_root_CMD} -b -q -l
          ${CMAKE_INSTALL_PREFIX}/tutorial/LoadOscProb.C TestMethods.C
  COMMAND ${CMAKE_COMMAND} -E env ${ROOT_root_CMD} -b -q -l
          ${CMAKE_INSTALL_PREFIX}/tutorial/LoadOscProb.C StressTest.C
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)

# Install the software
include(cmake/install.cmake)

# Properly package the software
include(cmake/package.cmake)

if(NOT PROJECT_IS_TOP_LEVEL)
  # do not show those options in the basic options as they are probably of
  # little interest for project including us
  mark_as_advanced(OSCPROB_ENABLE_TESTING)
endif()
