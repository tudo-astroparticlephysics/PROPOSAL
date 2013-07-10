# Load some basic macros which are needed later on
include(FindROOT.cmake)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


# Create a library called "libEvent" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)
 
link_directories( ${LINK_DIRECTORIES})

# List of source files
set(Examples_SRC
plot_range.cxx
)

# fill list of header files 
#set(Examples_HEADERS 
#Event.h
#)


############### build the library #####################
#add_library(Event ${Event_SRCS})

add_executable(Plot_range ${Examples_SRC})

target_link_libraries(Plot_range ${ROOT_LIBRARIES})
target_link_libraries(Plot_range PROPOSAL)



