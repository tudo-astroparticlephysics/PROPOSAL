cmake_minimum_required(VERSION 3.9)

set(PROPOSAL_VERSION_MAJOR 7)
set(PROPOSAL_VERSION_MINOR 0)
set(PROPOSAL_VERSION_PATCH 1)
set(PROPOSAL_VERSION ${PROPOSAL_VERSION_MAJOR}.${PROPOSAL_VERSION_MINOR}.${PROPOSAL_VERSION_PATCH})

project(PROPOSAL
    VERSION ${PROPOSAL_VERSION}
    LANGUAGES CXX
    DESCRIPTION "The very best photon and lepton propagator."
    )

if (EXISTS ${CMAKE_BINARY_DIR}/conan_paths.cmake)
    include(${CMAKE_BINARY_DIR}/conan_paths.cmake)
endif()

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type has been specified. Using default system build type. You may want to change this to"
  	"'Release' if you are using PROPOSAL for production purposes.")
else()
	if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  		message(WARNING "Build type set to ${CMAKE_BUILD_TYPE}. This build type should only be used for development"
  			"and debugging purposes and not for production since the runtime of PROPOSAL with be significantly higher!")
  	else()
		message(STATUS "Build type set to ${CMAKE_BUILD_TYPE}.")
  	endif()
endif()


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# removes tons of compiler warnings on clang (issue #141). there is probably a cleaner solution.
if(NOT BUILD_SHARED_LIBS)
	set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif()

option(BUILD_PYTHON "build python interface" OFF)
option(BUILD_EXAMPLE "build example" OFF)
option(BUILD_DOCUMENTATION "build documentation" OFF)
option(BUILD_TESTING "build testing" OFF)

add_subdirectory(src)

if(BUILD_EXAMPLE)
    add_subdirectory(example)
endif()

if(BUILD_DOCUMENTATION)
    add_subdirectory(docs)
endif()

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
