cmake_minimum_required(VERSION 3.9)

project(PROPOSAL VERSION 7.0.0 LANGUAGES CXX)

project(PROPOSAL LANGUAGES CXX)

set(PROPOSAL_VERSION_MAJOR 7)
set(PROPOSAL_VERSION_MINOR 0)
set(PROPOSAL_VERSION_PATCH 0)
set(PROPOSAL_VERSION ${PROPOSAL_VERSION_MAJOR}.${PROPOSAL_VERSION_MINOR}.${PROPOSAL_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

option(BUILD_PYTHON "build python interface" OFF)
option(BUILD_EXAMPLE "build example" OFF)
option(BUILD_DOCUMENTATION "build documentation" OFF)

add_subdirectory(src)

if(BUILD_EXAMPLE)
    add_subdirectory(example)
endif()

if(BUILD_DOCUMENTATION)
    add_subdirectory(docs)
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
