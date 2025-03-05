set(MATCH_SRC
        join.cpp adaptive.cpp)

add_library(matching SHARED ${MATCH_SRC})
target_link_libraries(matching plan)