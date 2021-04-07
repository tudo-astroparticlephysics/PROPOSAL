add_library(PROPOSALTestUtilities INTERFACE)


set(PROPOSALTestUtilities_HEADER
    "TestFilesHandling.h"
    )

set_target_properties(PROPOSALTestUtilities PROPERTIES
    PUBLIC_HEADER ${PROPOSALTestUtilities_HEADER})
