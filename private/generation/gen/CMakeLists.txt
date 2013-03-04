
i3_add_library(PROPOSALGeneration
        Atmosphere.cxx
        CosCorr.cxx
	CteqPDF.cxx
	EarthModel.cxx
	ExpCorr.cxx
	IntFlux.cxx
	NeutrinoInt.cxx
	NeutrinoTot.cxx
	AtmFlux.cxx
	USE_PROJECTS PROPOSAL
	)

target_link_libraries(PROPOSALGeneration ${I3_BUILD}/lib/libCteq61Pdf.so)
target_link_libraries(PROPOSALGeneration libgfortran.so )

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(GENERATION TestMain.cxx )
target_link_libraries(GENERATION PROPOSALGeneration)


