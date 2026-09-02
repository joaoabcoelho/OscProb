include_guard()

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

install(
  EXPORT oscprob
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oscprob
  NAMESPACE OscProb::)

# Add a symbolic link to the standard install include directory. This is mostly
# needed for dependencies that expects the include dir to be called inc as for
# the source.
install(
  CODE "
    message(STATUS \"Creating symlink inc -> include\")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink include inc
        WORKING_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}\"
    )
")

install(FILES ${CMAKE_CURRENT_LIST_DIR}/oscprob-config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oscprob)
