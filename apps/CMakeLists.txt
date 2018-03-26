project(apps)

add_executable(gmmregistration gmmPointSetRegistration.cxx)
target_link_libraries(gmmregistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmregistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(icpregistration icpPointSetRegistration.cxx)
target_link_libraries(icpregistration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(icpregistration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(gmmtransform gmm-transform.cxx)
target_link_libraries(gmmtransform ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(gmmtransform PUBLIC ${GMM_INCLUDE_DIRS})
