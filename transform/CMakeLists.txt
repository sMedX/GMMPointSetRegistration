set(_name Transform)
project(${_name})

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGenerateRandomTransform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkTransformAndAddNoiseMeshFilter.h
)

add_library(${_name} INTERFACE)
target_sources(${_name} INTERFACE ${HEADERS})

set(GMM_INCLUDE_DIRS ${GMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL "" FORCE)

set(GMM_LIBRARIES ${GMM_LIBRARIES} ${_name} CACHE INTERNAL "" FORCE)
