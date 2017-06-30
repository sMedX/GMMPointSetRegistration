project(apps)

add_executable(gmmPointSetRegistration gmmPointSetRegistration.cxx)
target_link_libraries(gmmPointSetRegistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmPointSetRegistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(itkICPPointSetRegistration itkICPPointSetRegistration.cxx)
target_link_libraries(itkICPPointSetRegistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(itkICPPointSetRegistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(gmmTransformMesh gmmTransformMesh.cxx)
target_link_libraries(gmmTransformMesh ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmTransformMesh PUBLIC ${GMM_INCLUDE_DIRS})
