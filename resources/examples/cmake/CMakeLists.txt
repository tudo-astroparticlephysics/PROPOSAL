cmake_minimum_required(VERSION 3.8)

project(proposal_cmake_example VERSION 0.1.0 LANGUAGES CXX)

find_package(PROPOSAL REQUIRED)

add_executable(example src/example.cpp)
target_link_libraries(example PRIVATE PROPOSAL::PROPOSAL)
