project(apps)

add_executable(registration registration.cxx)
target_link_libraries(registration ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(registration PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(transform transform.cxx)
target_link_libraries(transform ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(transform PUBLIC ${GMM_INCLUDE_DIRS})

add_executable(statistics statistics.cxx)
target_link_libraries(statistics ${ITK_LIBRARIES} ${GMM_LIBRARIES})
target_include_directories(statistics PUBLIC ${GMM_INCLUDE_DIRS})
