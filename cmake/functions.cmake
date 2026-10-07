include_guard()

function(configure_tables TABLES_DIR HEADER_DIR)
  set(PREMDIR ${TABLES_DIR}/PremTables)
  set(MODEL3DDIR ${TABLES_DIR}/EarthTables)
  set(PREMFILE ${PREMDIR}/prem_default.txt)
  set(PREM3DFILE ${MODEL3DDIR}/earth_binned_default.txt)

  configure_file(${PROJECT_SOURCE_DIR}/cmake/prem_default.hpp.in
                 ${HEADER_DIR}/prem_default.hpp @ONLY USE_SOURCE_PERMISSIONS)

  if(NOT "${TABLES_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}")
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/PremTables
            DESTINATION ${TABLES_DIR})

    install(DIRECTORY ${PROJECT_SOURCE_DIR}/EarthTables
            DESTINATION ${TABLES_DIR})
  endif()

endfunction()
