# Add library
add_library(OscProb)
# This is needed if we want to compile this project from a dependency with
# fetchContent
add_library(OscProb::OscProb ALIAS OscProb)

target_sources(
  OscProb
  PRIVATE Absorption.cxx
          EarthModelBase.cxx
          EarthModelBinned.cxx
          EigenPoint.cxx
          NuPath.cxx
          PMNS_Base.cxx
          PMNS_Maltoni.cxx
          PMNS_Decay.cxx
          PMNS_Deco.cxx
          PMNS_DensityMatrix.cxx
          PMNS_Fast.cxx
          PMNS_Iter.cxx
          PMNS_LIV.cxx
          PMNS_NSI.cxx
          PMNS_NUNM.cxx
          PMNS_OQS.cxx
          PMNS_SNSI.cxx
          PMNS_Sterile.cxx
          PMNS_Avg.cxx
          PMNS_SiderealLIV.cxx
          PremModel.cxx
          complexsolver.cxx)
target_include_directories(
  OscProb
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/inc>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:include>)

set(hdrs
    Absorption.h
    EarthModelBase.h
    EarthModelBinned.h
    EigenPoint.h
    NuPath.h
    PMNS_Base.h
    PMNS_Maltoni.h
    PMNS_Decay.h
    PMNS_Deco.h
    PMNS_DensityMatrix.h
    PMNS_Fast.h
    PMNS_Iter.h
    PMNS_LIV.h
    PMNS_NSI.h
    PMNS_NUNM.h
    PMNS_OQS.h
    PMNS_SNSI.h
    PMNS_Sterile.h
    PMNS_Avg.h
    PMNS_SiderealLIV.h
    PremModel.h
    Definitions.h)

set(hdrs_full_path ${hdrs})
list(TRANSFORM hdrs_full_path PREPEND "${PROJECT_SOURCE_DIR}/inc/")

# This configures prem_default.hpp
configure_tables(${CMAKE_INSTALL_PREFIX} ${PROJECT_BINARY_DIR})

list(APPEND hrds prem_default.hpp)
list(APPEND hdrs_full_path ${PROJECT_BINARY_DIR}/prem_default.hpp)

set_target_properties(OscProb PROPERTIES PUBLIC_HEADER "${hdrs_full_path}")

# Link libraries
target_link_libraries(
  OscProb
  PRIVATE MatrixDecomp
  PUBLIC ROOT::Core Eigen3::Eigen)

root_generate_dictionary(G__OscProb_dict ${hdrs} MODULE OscProb LINKDEF
                         "${PROJECT_SOURCE_DIR}/inc/LinkDef.h")
