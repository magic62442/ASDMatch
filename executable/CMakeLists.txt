set (PREPARE_SRC
        prepare.cpp)
add_executable(prepare.out ${PREPARE_SRC})
target_link_libraries(prepare.out relation plan)