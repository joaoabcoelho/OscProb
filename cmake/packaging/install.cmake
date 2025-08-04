include_guard()

set(PREMDIR ${CMAKE_INSTALL_PREFIX}/PremTables)
set(MODEL3DDIR ${CMAKE_INSTALL_PREFIX}/EarthTables)
set(PREMFILE ${PREMDIR}/prem_default.txt)
set(PREM3DFILE ${MODEL3DDIR}/earth_binned_default.txt)

configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/prem_default.hpp.in
  ${PROJECT_BINARY_DIR}/prem_default.hpp @ONLY USE_SOURCE_PERMISSIONS)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/PremTables
        DESTINATION ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/EarthTables
        DESTINATION ${CMAKE_INSTALL_PREFIX})

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
