#!/usr/bin/env python

import os

from icecube.dataclasses import I3Particle
from icecube.PROPOSAL import I3PropagatorServicePROPOSAL


# FIXME - The table path I3PropagatorServicePROPOSAL settles
#         on behind the scenes is not always going to be 
#         $I3_BUILD/PROPOSAL/resources/tables/
#         It would be nice then to be able to instantiate an
#         I3PropagatorServicePROPOSAL object, figure out what
#         path settled on so we can make these checkes here.
#         1) Add a 'generate_tables' method to I3PropagatorServicePROPOSAL
#         2) Don't include expensive, blocking processes in constructors.
#         3) Expose the path it expects to find tables in, so
#            we can instantiate multiple services in different processes
#            without them wrestling with each other.

# This should have generated all data files.
# Generate a marker file to tell the build system this task
# has been successfully completed.
file_lock_path = os.path.expandvars("$I3_BUILD/PROPOSAL/resources/tables/.tables.auto_generated")
if not os.path.exists(file_lock_path):
    with open(file_lock_path, 'w') as f:
        try:
            I3PropagatorServicePROPOSAL()
        except:
            # something went wrong and table generation was not successful
            os.remove(file_lock_path)
        f.write('success') # need to be able to communicate to other objects instantiated
        # in different processes whether the tables it needs are good to go.  if not, it
        # should probably emit an error message "tables are not ready yet."
else:
    print("%s exists already, indicating the tables are already (or in the process of being) generated." % file_lock_path)
