
# GTest dependency
find_package(Threads REQUIRED)

# add pbcopper, pbcopper test, & GoogleTest include dirs
include_directories(
    ${pbcopper_INCLUDE_DIRS}
    ${pbcopper_TestsDir}/include
    ${pbcopper_ThirdPartyDir}/googletest
)

# grab unit test source files
include(files.cmake)
set(SOURCES
    ${pbcopperTest_H}
    ${pbcopperTest_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${pbcopper_CXX_FLAGS}")

# define unit test exe
add_executable(test_pbcopper EXCLUDE_FROM_ALL ${SOURCES})
set_target_properties(test_pbcopper
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)
target_link_libraries(test_pbcopper
    ${pbcopper_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

# add unit tests exe to CTest framework
add_test(
    NAME PbcopperUnitTests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
    COMMAND test_pbcopper
)

# add 'make check' for easy 'build & run tests'
set(verbose_ctest_command ctest -VV)
add_custom_target(check COMMAND ${verbose_ctest_command})
add_dependencies(check
    # list all explicitly to get re-compile on any changes to lib or tests
    test_pbcopper
    pbcopper
)

# add 'make check-quiet' for 'make check' w/ simple pass/fail summary
add_custom_target(check-quiet COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check-quiet
    # list all explicitly to get re-compile on any changes to lib or tests
    test_pbcopper
    pbcopper
)

# add 'make check-xunit' for 'make check' w/ XUnit output (on Bamboo)
add_custom_target(check-xunit
    COMMAND test_pbcopper --gtest_output=xml:${CMAKE_BINARY_DIR}/pbcopper-unit.xml
)
add_dependencies(check-xunit
    # list all explicitly to get re-compile on any changes to lib or tests
    test_pbcopper
    pbcopper
)
