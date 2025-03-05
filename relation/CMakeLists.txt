set(RELATION_SRC
        candidate_space.cpp
        graph.cpp
        relation.cpp
        dynamic_array.cpp dynamic_array.h)

add_library(relation SHARED ${RELATION_SRC})
target_link_libraries(relation utility)