cmake_minimum_required(VERSION 3.12)
project(in_memory_join)

set(CMAKE_CXX_STANDARD 17)
IF(NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
ENDIF()

include_directories(
        ${PROJECT_SOURCE_DIR}/config
        ${PROJECT_SOURCE_DIR}/utility
        ${PROJECT_SOURCE_DIR}/relation
        ${PROJECT_SOURCE_DIR}/plan
        ${PROJECT_SOURCE_DIR}/matching
)

add_subdirectory(config)
add_subdirectory(utility)
add_subdirectory(relation)
add_subdirectory(plan)
add_subdirectory(matching)
add_subdirectory(executable)
add_subdirectory(utility/td)

# link to GLPK
include_directories(${PROJECT_SOURCE_DIR}/utility/td/glpk_lib/include)
link_directories(${PROJECT_SOURCE_DIR}/utility/td/glpk_lib/lib)

add_executable(ASDMatch.out main.cpp)
target_link_libraries(ASDMatch.out matching)
