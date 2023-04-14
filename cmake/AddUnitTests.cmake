
function(mopsa_add_unittest test_source)

  get_filename_component(binname ${test_source} NAME_WE)

  message("[+] AddUnitTests ${test_source} ${binname}")

  add_executable(${binname} ${test_source})
  target_include_directories(${binname} PUBLIC ${PROJECT_SOURCE_DIR} "${gtest_SOURCE_DIR}/include")
  target_link_libraries(${binname} mopsa gtest_main gtest)

  #set_target_properties(readverilogtest PROPERTIES
    ##RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/unittests_bin"
  #)

  gtest_discover_tests(
    ${binname}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )

endfunction()
