set(files colormod.h Utils.h MakeTestSamples.C StressTest.C TestMethods.C)

# add_library(Test) target_sources(Test PRIVATE Utils.h)
# set_target_properties(Test PROPERTIES LINKER_LANGUAGE CXX)
# target_link_libraries(Test PRIVATE OscProb)

# file(COPY ${files} DESTINATION ${PROJECT_BINARY_DIR}/stage/test)
# file(CREATE_LINK ${CMAKE_CURRENT_LIST_DIR}/data
# ${PROJECT_BINARY_DIR}/stage/test/data SYMBOLIC)

# add_test( NAME TestMethods COMMAND ${ROOT_root_CMD} -b -q
# ${PROJECT_BINARY_DIR}/stage/tutorial/LoadOscProb.C
# ${PROJECT_BINARY_DIR}/stage/test/TestMethods.C WORKING_DIRECTORY
# ${PROJECT_SOURCE_DIR}/test) add_test( NAME StressTest COMMAND ${ROOT_root_CMD}
# -b -q ${PROJECT_BINARY_DIR}/stage/tutorial/LoadOscProb.C
# ${PROJECT_BINARY_DIR}/stage/test/StressTest.C WORKING_DIRECTORY
# ${PROJECT_SOURCE_DIR}/test)

install(FILES ${files} DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_PREFIX}/test/)
