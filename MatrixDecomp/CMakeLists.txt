add_library(MatrixDecomp)

target_sources(MatrixDecomp PRIVATE zheevc3.cxx zheevh3.cxx zheevq3.cxx
                                    zhetrd3.cxx)

target_include_directories(
  MatrixDecomp PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                      $<INSTALL_INTERFACE:include/MatrixDecomp>)

list(APPEND PUB_HDRS_MD zheevc3.h zheevh3.h zheevq3.h zhetrd3.h)
set_target_properties(MatrixDecomp PROPERTIES PUBLIC_HEADER "${PUB_HDRS_MD}"
                                              OUTPUT_NAME MatrixDecomp)
target_link_libraries(MatrixDecomp PUBLIC ROOT::RIO)

add_root_dictionary(MatrixDecomp HEADERS ${PUB_HDRS_MD} LINKDEF LinkDef.h)

install(
  TARGETS MatrixDecomp
  EXPORT oscprob
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MatrixDecomp)
