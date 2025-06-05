set(SOURCE_FILES
        dataloader.cpp
        hypergraph.cpp
        preprocessing.cpp
        td.cpp
        twbound.cpp
        SCsolver.cpp
        td_utils.hpp
)

set(GLPK_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glpk_lib)
set(GLPK_INCLUDE_DIR ${GLPK_ROOT_DIR}/include)
set(GLPK_LIB_DIR ${GLPK_ROOT_DIR}/lib)

add_library(td SHARED ${SOURCE_FILES})
#include_directories(${GLPK_INCLUDE_DIR})
target_include_directories(td PUBLIC
        ${GLPK_INCLUDE_DIR}
)

target_link_directories(td PUBLIC
        ${GLPK_LIB_DIR}
)

find_library(GLPK_LIB
        NAMES glpk glpk.dylib
        PATHS ${GLPK_LIB_DIR}
        NO_DEFAULT_PATH
)
if(NOT GLPK_LIB)
    message(FATAL_ERROR "GLPK library not found in: ${GLPK_LIB_DIR}")
endif()
target_link_libraries(td PUBLIC ${GLPK_LIB})

set_target_properties(td PROPERTIES BUILD_RPATH "${GLPK_LIB_DIR}" INSTALL_RPATH "${GLPK_LIB_DIR}")

message(STATUS "GLPK include dir: ${GLPK_INCLUDE_DIR}")
message(STATUS "GLPK library dir: ${GLPK_LIB_DIR}")