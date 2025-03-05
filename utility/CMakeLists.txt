set(UTILITY_SRC
    command_parser.cpp utils.cpp command.cpp compute_set_intersection.cpp td/preprocessing.h clique.cpp)


add_library(utility SHARED
${UTILITY_SRC})

target_link_libraries(utility td)