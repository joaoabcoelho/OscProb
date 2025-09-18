include_guard()

#
# add_root_dictionary generates one dictionary to be added to a target.
#
# this function is not meant to be called by the user, but by other CMake
# function, e.g. oscProbAddRootDictionary, as the target name is not the user
# friendly one, but the internal namespaced-and-all one.
#
# Besides the dictionary source itself a PCM file is also generated. It will be
# installed alongside the target's library file
#
# arguments :
#
# * 1st parameter (required) is the target the dictionary should be added to
#
# * HEADERS (required) is a list filepaths needed for the dictionary definition.
#
# * LINKDEF (required) the LINKDEF file needed by rootcling.
#
# HEADERS and LINKDEF paths can be absolute or relative. If they are relative
# they are assumed to be relative to the current CMakeLists.txt (aka
# CMAKE_CURRENT_LIST_DIR).
#
# The target must be of course defined _before_ calling this function (i.e.
# add_library(target ...) has been called), and target_include_directories
# _must_ have be called as well, in order to be able to compute the list of
# include directories needed to _compile_ the dictionary
#
# Note also that the generated dictionary is added to PRIVATE SOURCES list of
# the target.
#

function(add_root_dictionary target)
  cmake_parse_arguments(PARSE_ARGV 1 A "" "LINKDEF" "HEADERS")
  if(A_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR "Unexpected unparsed arguments: ${A_UNPARSED_ARGUMENTS}")
  endif()

  set(required_args "LINKDEF;HEADERS")
  foreach(required_arg IN LISTS required_args)
    if(NOT A_${required_arg})
      message(FATAL_ERROR "Missing required argument: ${required_arg}")
    endif()
  endforeach()

  # # convert all relative paths to absolute ones. LINKDEF must be the last one.
  # foreach(h IN LISTS A_HEADERS A_LINKDEF) cmake_path(IS_RELATIVE h
  # is_relative) if(${is_relative}) cmake_path( ABSOLUTE_PATH h BASE_DIRECTORY
  # ${CMAKE_CURRENT_LIST_DIR} NORMALIZE OUTPUT_VARIABLE habs) list(APPEND
  # headers ${habs}) else() list(APPEND headers ${h}) endif() endforeach()

  get_property(
    outputName
    TARGET ${target}
    PROPERTY OUTPUT_NAME)
  # set(dictionaryFile G__${outputName}.cc)

  root_generate_dictionary(G__${outputName}_dict ${A_HEADERS} MODULE ${target}
                           LINKDEF ${A_LINKDEF})

  # add_custom_command( OUTPUT ${dictionaryFile} G__${outputName}_rdict.pcm
  # COMMAND ${CMAKE_BINARY_DIR}/rootcling_wrapper.sh --rootcling_cmd
  # $<TARGET_FILE:ROOT::rootcling> --dictionary_file ${dictionaryFile}
  # --include_dirs
  # $<LIST:REMOVE_DUPLICATES,$<TARGET_PROPERTY:${target},INCLUDE_DIRECTORIES>>
  # --headers "${headers}" --pcm_destination $<TARGET_FILE_DIR:${target}>
  # DEPENDS "${headers}" VERBATIM)

  # target_sources(${target} PRIVATE ${dictionaryFile})

  install(FILES $<TARGET_FILE_DIR:${target}>/lib${outputName}_rdict.pcm
          TYPE LIB)

endfunction()
