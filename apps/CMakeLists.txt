project(apps)

add_executable(gmmPointSetRegistration gmmPointSetRegistration.cxx)
target_link_libraries(gmmPointSetRegistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmPointSetRegistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(icpPointSetRegistration icpPointSetRegistration.cxx)
target_link_libraries(icpPointSetRegistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(icpPointSetRegistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(gmmTransformMesh gmmTransformMesh.cxx)
target_link_libraries(gmmTransformMesh ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmTransformMesh PUBLIC ${GMM_INCLUDE_DIRS})
