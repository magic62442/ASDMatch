set(PLAN_SRC
       decomposition.cpp estimator.cpp optimizer.cpp subset_structure.cpp subset_structure.h)

add_library(plan SHARED ${PLAN_SRC})
target_link_libraries(plan relation)