set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )

if (CMAKE_BUILD_TYPE MATCHES DEBUG)
  message ("Debug mode")
  add_executable(main3d-d main3d.cpp)
  add_executable(evaluate3d-d evaluate3d.cpp)
  add_executable(generate3DSorter-d generate3DSorter.cpp)

  target_link_libraries(main3d-d mopsa)
  target_link_libraries(evaluate3d-d mopsa)
  target_link_libraries(generate3DSorter-d mopsa)
else()
  message ("Release mode")
  add_executable(main3d main3d.cpp)
  add_executable(evaluate3d evaluate3d.cpp)
  add_executable(generate3DSorter generate3DSorter.cpp)

  target_link_libraries(main3d mopsa)
  target_link_libraries(evaluate3d mopsa)
  target_link_libraries(generate3DSorter mopsa)
endif()


