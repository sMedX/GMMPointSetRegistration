set(_name GMM)
project(${_name})

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMICPPointSetToPointSetMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMKCPointSetToPointSetMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMMLEPointSetToPointSetMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetMetricBase.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetRegistrationMethod.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMRegstrationUtils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMML2RigidPointSetToPointSetMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkInitializeMetric.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkInitializeTransform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkNormalizePointSet.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkPointSetPropertiesCalculator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkPointSetToPointSetMetrics.h
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMICPPointSetToPointSetMetric.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMKCPointSetToPointSetMetric.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMMLEPointSetToPointSetMetric.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetMetric.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetMetricBase.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMMPointSetToPointSetRegistrationMethod.hxx
    ${CMAKE_CURRENT_SOURCE_DIR}/itkGMML2RigidPointSetToPointSetMetric.hxx
)

add_library(${_name} INTERFACE)
target_sources(${_name} INTERFACE ${HEADERS})

set(GMM_INCLUDE_DIRS ${GMM_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE INTERNAL "" FORCE)

set(GMM_LIBRARIES ${GMM_LIBRARIES} ${_name} CACHE INTERNAL "" FORCE)
