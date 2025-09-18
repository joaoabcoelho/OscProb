include_guard()

include(FetchContent)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO)

# This code could be replaced with:
# FetchContent_Declare(Eigen3...FIND_PACKAGE_ARGS), but this feature was added
# only on version 3.24. Since the KOFI project is currently compiled with an old
# image shipping cmake 3.22.1 we use this old way to decide if the fetch is
# needed or not. This hack can be improved as soon as we can upgrade the minimum
# cmake version.
find_package(Eigen3 QUIET)
if(NOT TARGET Eigen3::Eigen)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  # The CMakeLists of eigen3 version 3.4.0 is not well optimized and attempts in
  # compiling lots of targets that we do not use. This code allows at least to
  # disable the testing, while still being able to compile our tests if
  # required. Notice that these lines can be removed when a new eigen version
  # will be used, since the CMakeLists are better optimized in newer versions.
  if(BUILD_TESTING)
    set(OLD_BUILD_TESTING ${BUILD_TESTING})
  endif()
  set(BUILD_TESTING OFF)
  set(EIGEN_BUILD_DOC OFF)
  FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        EXCLUDE_FROM_ALL)
  FetchContent_MakeAvailable(Eigen3)
  if(OLD_BUILD_TESTING)
    set(BUILD_TESTING ${OLD_BUILD_TESTING})
  else()
    unset(BUILD_TESTING)
  endif()
endif()
