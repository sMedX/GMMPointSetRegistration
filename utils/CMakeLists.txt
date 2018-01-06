set(_name GMMUtils)
project(${_name})

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/itkIOutils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/argsCustomParsers.h
)

add_library(${_name} INTERFACE)
target_sources(${_name} INTERFACE ${HEADERS})

set(GMM_INCLUDE_DIRS ${GMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL "" FORCE)

set(GMM_LIBRARIES ${GMM_LIBRARIES} ${_name} CACHE INTERNAL "" FORCE)
