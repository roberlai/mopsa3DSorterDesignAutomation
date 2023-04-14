
add_library(mopsa)

## BOOT
find_package(Boost)
target_include_directories(mopsa PUBLIC ${Boost_INCLUDE_DIR} )
target_link_libraries(mopsa ${Boost_LIBRARIES})

## OMP
find_package(OpenMP)
target_link_libraries(mopsa PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(mopsa PUBLIC ${PROJECT_SOURCE_DIR})

target_sources(mopsa PRIVATE headerdef.hpp)
target_include_directories(mopsa PUBLIC ${PROJECT_SOURCE_DIR}/mopsa)

add_subdirectory(logger)
add_subdirectory(util)
add_subdirectory(chip)
add_subdirectory(geometry)
add_subdirectory(sim)
add_subdirectory(gen)
