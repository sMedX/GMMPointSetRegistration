cmake_minimum_required(VERSION 3.3)

project(GMMRegistration)

set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

# Find ITK
find_package(ITK 4.9 REQUIRED)
include(${ITK_USE_FILE})

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(GMM_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/thirdparty
    CACHE INTERNAL "" FORCE
)

set(GMM_LIBRARIES "" CACHE INTERNAL "" FORCE)

add_subdirectory(${CMAKE_SOURCE_DIR}/gmm)
add_subdirectory(${CMAKE_SOURCE_DIR}/transform)
add_subdirectory(${CMAKE_SOURCE_DIR}/mesh)
add_subdirectory(${CMAKE_SOURCE_DIR}/initializer)
add_subdirectory(${CMAKE_SOURCE_DIR}/utils)
add_subdirectory(${CMAKE_SOURCE_DIR}/application)
