cmake_minimum_required(VERSION 3.2)

project(GMMRegistration)

set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

# Find ITK
find_package(ITK 4.9 REQUIRED)
include(${ITK_USE_FILE})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(gmmPointSetRegistration gmmPointSetRegistration.cxx)
target_link_libraries(gmmPointSetRegistration ${ITK_LIBRARIES})

add_executable(gmmTransformMesh gmmTransformMesh.cxx)
target_link_libraries(gmmTransformMesh ${ITK_LIBRARIES})
