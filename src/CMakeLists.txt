# ===== src directory =====

# define variable PROPOSAL_Sources for the PROPOSAL sources

set(PROPOSAL_Sources Integral DFunctionOfx Interpolate MathModel
                FunctionOfx PhysicsModel Medium
                methods Output  Ionizationloss IonizStochastic
                IonizContinuous  Epairproduction EpairContinuous
                EpairStochastic Photonuclear PhotoStochastic
                PhotoContinuous Bremsstrahlung BremsStochastic
                BremsContinuous Decay FindRoot Energy2Loss
                Energy2LossE Energy2LossX StandardNormal Scattering
                CrossSections PROPOSALParticle Propagate MMCexception
		Mpairproduction MpairContinuous MpairStochastic
		Amanda)
# define sources for executable
set(EXE_Sources TestMain)


# create and fill the build/include directory
foreach(source ${PROPOSAL_Sources})
        configure_file(${PROPOSAL_SOURCE_DIR}/src/${source}.h ${PROPOSAL_BINARY_DIR}/public/${source}.h COPYONLY)
endforeach(source ${PROPOSAL_Sources})

# include files (.h files) directories
include_directories( ${PROPOSAL_BINARY_DIR}/public )

# Copy medidef-file to the build directory
configure_file(${PROPOSAL_SOURCE_DIR}/src/mediadef ${PROPOSAL_BINARY_DIR}/src/mediadef COPYONLY)

# Create a library called "PROPOSALCore"
# which includes the source files "*.cxx"


add_library (PROPOSAL ${PROPOSAL_Sources})

# Link the library PROPOSALCore to the libraries
# ROOT version 5.26
#target_link_libraries(PROPOSALCore Core Hist Thread Tree Gpad RIO Cint MathCore Matrix Net Graf dl)

# Create executable

add_executable(PROPOSAL_EXE ${EXE_Sources} )

# Link the library PROPOSALCore to the executable
target_link_libraries(PROPOSAL_EXE PROPOSAL )
 
