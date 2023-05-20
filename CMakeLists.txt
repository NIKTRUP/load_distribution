cmake_minimum_required(VERSION 3.24)
project(load_distribution)


set(CMAKE_CXX_STANDARD 17)

add_executable(load_distribution main.cpp
        include/membrane_load_distribution.h
        include/write.h
        include/functions.h
        include/log_duration.h
        include/string_load_distribution.h
        include/test_framework.h
        include/tests.h
        src/functions.cpp
        src/membrane_load_distribution.cpp
        src/string_load_distribution.cpp
        src/tests.cpp
        )

#find_package(CUDA)
#if(CUDA_CXX_FOUND)
#    target_link_libraries(load_distribution PUBLIC CUDA::CUDA_CXX)
#endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(load_distribution PUBLIC OpenMP::OpenMP_CXX)
endif()
