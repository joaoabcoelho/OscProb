include_guard()

include(FetchContent)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Tree)

set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3147391d946bb4b6c68edd901f2add6ac1f31f8c
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  EXCLUDE_FROM_ALL FIND_PACKAGE_ARGS)
FetchContent_MakeAvailable(Eigen3)
