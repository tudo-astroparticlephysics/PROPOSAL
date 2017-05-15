CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

ADD_DEFINITIONS(-DPROPOSAL_STANDALONE=1)

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

SET(DEBUG OFF)
OPTION (DEBUG "DEBUG" OFF)
IF(DEBUG)
        SET( CMAKE_CXX_FLAGS "-g -pg -O0" )
ENDIF()

SET(ICETRAY_INCLUDE_PATH "/home/koehne/Simulation/icesim4_candidate/V04-00-01-RC/icetray/public/")


#################################################################
#################           boost       #########################
#################################################################

FIND_PACKAGE( Boost COMPONENTS program_options REQUIRED )
SET(LIBRARYS_TO_LINK ${Boost_LIBRARIES})

#################################################################
#################           python      #########################
#################################################################

# SET(ADD_PYTHON ON CACHE BOOL "Choose to compile the python wrapper library")
OPTION (ADD_PYTHON "Choose to compile the python wrapper library" ON)

IF(ADD_PYTHON)
	MESSAGE(STATUS "Enabled to build the python wrapper library.")
	FIND_PACKAGE( PythonLibs 2.7 REQUIRED )

	IF(PYTHONLIBS_FOUND)
		include_directories( ${PYTHON_INCLUDE_DIRS} )
	ENDIF(PYTHONLIBS_FOUND)
ELSE(ADD_PYTHON)
	MESSAGE(STATUS "No python wrapper library will be build.")
ENDIF(ADD_PYTHON)

#################################################################
#################           ROOT        #########################
#################################################################

# Load some basic macros which are needed later on
INCLUDE(FindROOT.cmake)

#if ROOT is found ROOT files with ROOT trees can be written
IF(ROOT_FOUND)
    ADD_DEFINITIONS(-DROOT_SUPPORT=1)

    SET(INCLUDE_DIRECTORIES
    ${ROOT_INCLUDE_DIR}
    )

    INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES})

    SET(LINK_DIRECTORIES
    ${ROOT_LIBRARY_DIR}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    )

    LINK_DIRECTORIES( ${LINK_DIRECTORIES})
    SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${ROOT_LIBRARIES})

ELSE(ROOT_FOUND)

    ADD_DEFINITIONS(-DROOT_SUPPORT=0)

    MESSAGE(STATUS "ROOT not found...")
    MESSAGE(STATUS "ROOT examples will not be builded.")
    MESSAGE(STATUS "No ROOT Output is available.")
    MESSAGE(STATUS "Make sure you have ROOT installed and ROOTSYS is set.")

ENDIF(ROOT_FOUND)


#################################################################
#################           log4cplus     #######################
#################################################################

INCLUDE(resources/FindLog4cplus.cmake)

IF(LOG4CPLUS_FOUND)
    ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=1)
	ADD_DEFINITIONS(-DLOG4CPLUS_CONFIG=\"${PROJECT_SOURCE_DIR}/resources/log4cplus.conf\")
    SET(LIBRARYS_TO_LINK ${LIBRARYS_TO_LINK} ${LOG4CPLUS_LIBRARIES})
ELSE(LOG4CPLUS_FOUND)
    ADD_DEFINITIONS(-DLOG4CPLUS_SUPPORT=0)
    MESSAGE(STATUS "No logging will be done.")
ENDIF(LOG4CPLUS_FOUND)


#################################################################
#################           GTest       #########################
#################################################################

FIND_PACKAGE(GTest)

if (GTEST_FOUND)
	INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
 	ENABLE_TESTING()
else (GTEST_FOUND)
	SET(DO_TESTING FALSE)
    MESSAGE(STATUS "No tests will be build.")
endif (GTEST_FOUND)

#################################################################
#################           Libraries    ########################
#################################################################

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/public" "${PROJECT_SOURCE_DIR}" ${LOG4CPLUS_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

FILE(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/private/PROPOSAL/*.cxx)
ADD_LIBRARY(PROPOSAL SHARED
		${SRC_FILES}
        # private/PROPOSAL/Integral.cxx
        # private/PROPOSAL/methods.cxx
        # private/PROPOSAL/MathModel.cxx
        # private/PROPOSAL/Bremsstrahlung.cxx
        # private/PROPOSAL/CrossSections.cxx
        # private/PROPOSAL/Decay.cxx
        # private/PROPOSAL/Epairproduction.cxx
        # private/PROPOSAL/Ionization.cxx
        # private/PROPOSAL/Photonuclear.cxx
        # private/PROPOSAL/Medium.cxx
        # private/PROPOSAL/PROPOSALParticle.cxx
        # private/PROPOSAL/EnergyCutSettings.cxx
		# private/PROPOSAL/Interpolant.cxx
        # private/PROPOSAL/StandardNormal.cxx
        # private/PROPOSAL/RootFinder.cxx
        # private/PROPOSAL/ProcessCollection.cxx
        # private/PROPOSAL/Propagator.cxx
        # private/PROPOSAL/ContinuousRandomization.cxx
        # private/PROPOSAL/Geometry.cxx
        # private/PROPOSAL/Scattering.cxx
        # private/PROPOSAL/ScatteringFirstOrder.cxx
        # private/PROPOSAL/ScatteringMoliere.cxx
        # private/PROPOSAL/Output.cxx
	)

SET_TARGET_PROPERTIES(PROPOSAL PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O2")# -Wextra -pedantic")
TARGET_LINK_LIBRARIES(PROPOSAL ${LIBRARYS_TO_LINK})

ADD_EXECUTABLE(PROPOSALtest
        private/test/PROPOSAL.cxx
)

SET_TARGET_PROPERTIES(PROPOSALtest PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(PROPOSALtest PROPOSAL)

EXECUTE_PROCESS(COMMAND ln -s ${CMAKE_SOURCE_DIR}/resources ${CMAKE_BINARY_DIR}/resources)

ADD_EXECUTABLE(WriteSectorsFromDomList
        private/test/WriteSectorsFromDomList.cxx
)

TARGET_LINK_LIBRARIES(WriteSectorsFromDomList PROPOSAL)

INSTALL(TARGETS PROPOSAL DESTINATION lib)

FILE(GLOB INC_FILES ${PROJECT_SOURCE_DIR}/public/PROPOSAL/*.h)
INSTALL(FILES ${INC_FILES} DESTINATION include/PROPOSAL)

#################################################################
#################           Tests        ########################
#################################################################

IF(DO_TESTING)
  EXECUTE_PROCESS(COMMAND mkdir ${PROPOSAL_BINARY_DIR}/bin/ OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

  #create tar directory with "tar -czvf TestFiles.tar.Z TestFiles/" and put it in Test directory
  EXECUTE_PROCESS(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/Test/TestFiles.tar.Z -C ${PROPOSAL_BINARY_DIR}/bin/
                    OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)
  EXECUTE_PROCESS(COMMAND  tar -xvf ${PROJECT_SOURCE_DIR}/Test/TestFiles2.tar.Z -C ${PROPOSAL_BINARY_DIR}/bin/
                    OUTPUT_VARIABLE _output OUTPUT_STRIP_TRAILING_WHITESPACE)

  ADD_EXECUTABLE(UnitTest_Scattering Test/Scattering_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_StandardNormal Test/StandardNormal_Test.cxx)
  ADD_EXECUTABLE(UnitTest_Photonuclear Test/Photonuclear_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Integral Test/Integral_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Interpolant Test/Interpolant_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Bremsstrahlung Test/Bremsstrahlung_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Epairproduction Test/Epairproduction_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Ionization Test/Ionization_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_RootFinder Test/RootFinder_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Medium Test/Medium_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Particle Test/Particle_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_EnergyCutSettings Test/EnergyCutSettings_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Decay Test/Decay_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_ProcessCollection Test/ProcessCollection_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_ContinuousRandomization Test/ContinuousRandomization_TEST.cxx)
  ADD_EXECUTABLE(UnitTest_Geometry Test/Geometry_TEST.cxx)

  TARGET_LINK_LIBRARIES(UnitTest_Scattering PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_StandardNormal PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Integral PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Interpolant PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Ionization PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Bremsstrahlung PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Epairproduction PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Photonuclear PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_RootFinder PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Medium PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Particle PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_EnergyCutSettings PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Decay PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_ProcessCollection PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_ContinuousRandomization PROPOSAL gtest gtest_main)
  TARGET_LINK_LIBRARIES(UnitTest_Geometry PROPOSAL gtest gtest_main)

  ADD_TEST(UnitTest_Scattering bin/UnitTest_Scattering)
  ADD_TEST(UnitTest_ContinuousRandomization bin/UnitTest_ContinuousRandomization)
  ADD_TEST(UnitTest_ProcessCollection bin/UnitTest_ProcessCollection)
  ADD_TEST(UnitTest_Decay bin/UnitTest_Decay)
  ADD_TEST(UnitTest_RootFinder bin/UnitTest_RootFinder)
  ADD_TEST(UnitTest_StandardNormal bin/UnitTest_StandardNormal)
  ADD_TEST(UnitTest_Integral bin/UnitTest_Integral)
  ADD_TEST(UnitTest_Medium bin/UnitTest_Medium)
  ADD_TEST(UnitTest_Particle bin/UnitTest_Particle)
  ADD_TEST(UnitTest_EnergyCutSettings bin/UnitTest_EnergyCutSettings)
  ADD_TEST(UnitTest_Interpolant bin/UnitTest_Interpolant)
  ADD_TEST(UnitTest_Epairproduction bin/UnitTest_Epairproduction)
  ADD_TEST(UnitTest_Ionization bin/UnitTest_Ionization)
  ADD_TEST(UnitTest_Bremsstrahlung bin/UnitTest_Bremsstrahlung)
  ADD_TEST(UnitTest_Photonuclear bin/UnitTest_Photonuclear)
  ADD_TEST(UnitTest_Geometry bin/UnitTest_Geometry)

ENDIF()

ADD_SUBDIRECTORY( doc )

IF(ROOT_FOUND)
    ADD_SUBDIRECTORY( root_examples )
ENDIF(ROOT_FOUND)

IF(ADD_PYTHON)
	IF(PYTHONLIBS_FOUND)
		ADD_LIBRARY(pyPROPOSAL SHARED private/python/pybindings.cxx)
		SET_TARGET_PROPERTIES(pyPROPOSAL PROPERTIES PREFIX "" COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
		TARGET_LINK_LIBRARIES(pyPROPOSAL boost_python ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} PROPOSAL)
		INSTALL(TARGETS pyPROPOSAL DESTINATION lib)
	ENDIF(PYTHONLIBS_FOUND)
ENDIF(ADD_PYTHON)

