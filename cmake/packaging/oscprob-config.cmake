include(CMakeFindDependencyMacro)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
find_dependency(Eigen3)
include(${CMAKE_CURRENT_LIST_DIR}/oscprob.cmake)
