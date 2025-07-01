SET(CMAKE_CXX_STANDARD 17)

include(GoogleTest)

SET(RESOURCE_URL https://proposal.app.tu-dortmund.de/resources)
SET(TEST_FILES TestFiles-2023-01-05.tar.gz)

# download the test files
add_custom_command(OUTPUT TestFiles
    COMMAND curl --silent ${RESOURCE_URL}/${TEST_FILES} -o TestFiles.tar.gz
    COMMAND ${CMAKE_COMMAND} -E tar -xzf TestFiles.tar.gz
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Downloading and extracting test files."
    )

add_custom_target(testfiles DEPENDS TestFiles)

# small lib with utilities for PROPOSAL testing
add_library(TestUtilities INTERFACE)
add_library(PROPOSAL::TestUtilities ALIAS TestUtilities)
target_include_directories(TestUtilities INTERFACE PROPOSALTestUtilities)

find_package(GTest REQUIRED)
find_package(Boost REQUIRED filesystem)
# creat test binary and link against proposal and gtest
macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    add_dependencies(${TESTNAME} testfiles)
    target_link_libraries(${TESTNAME} PROPOSAL::PROPOSAL PROPOSAL::TestUtilities
        GTest::GTest Boost::filesystem)
    gtest_discover_tests(${TESTNAME})
endmacro()

package_add_test(UnitTest_Epairproduction Epairproduction_TEST.cxx)

