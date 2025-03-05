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
IF(NOT APPLE)
    include_directories($ENV{HOME}/.local/include)
    link_directories($ENV{HOME}/.local/lib)
ELSE()
    include_directories(/opt/homebrew/Cellar/glpk/5.0/include)
    link_directories(/opt/homebrew/Cellar/glpk/5.0/lib)
ENDIF()

add_executable(ASDMatch.out main.cpp)
target_link_libraries(ASDMatch.out matching)
