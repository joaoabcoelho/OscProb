add_library(MatrixDecomp)
target_sources(
  MatrixDecomp
  PRIVATE zheevc3.cxx zheevh3.cxx zheevq3.cxx zhetrd3.cxx
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${PROJECT_SOURCE_DIR}
         FILES
         ${CMAKE_CURRENT_SOURCE_DIR}/zheevc3.h
         ${CMAKE_CURRENT_SOURCE_DIR}/zheevh3.h
         ${CMAKE_CURRENT_SOURCE_DIR}/zheevq3.h
         ${CMAKE_CURRENT_SOURCE_DIR}/zhetrd3.h)

target_link_libraries(MatrixDecomp PUBLIC ROOT::RIO)

set_target_properties(MatrixDecomp PROPERTIES OUTPUT_NAME MatrixDecomp)
get_target_property(hdrs MatrixDecomp HEADER_SET)
add_root_dictionary(MatrixDecomp HEADERS ${hdrs} LINKDEF LinkDef.h)

install(
  TARGETS MatrixDecomp
  EXPORT oscprob
  FILE_SET HEADERS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
