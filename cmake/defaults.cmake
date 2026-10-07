include(GNUInstallDirs)

# Instructions to properly handle RPATH
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
endif()

#
# Define our RPATH settings, for both build and install phases
#

# use, i.e. do not skip, the full RPATH in the _build_ tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, do not use the install RPATH already (will only be used when
# actually installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add to the install RPATH the (automatically determined) parts of the RPATH
# that point to directories outside the build tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# specify libraries directory relative to binaries one.
file(RELATIVE_PATH relDir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
     ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

if(APPLE)
  set(basePoint @loader_path)
else()
  set(basePoint $ORIGIN)
endif()

set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir}
                        ${basePoint}/../../lib)
