
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
	load_project("PROPOSAL", false);
	
	using namespace boost::python;
	import("icecube.sim_services");
	
	class_<I3PropagatorServicePROPOSAL, shared_ptr<I3PropagatorServicePROPOSAL>,
	    bases<I3PropagatorService>, boost::noncopyable>(
	    "I3PropagatorServicePROPOSAL", init<  const std::string& , optional<bool> > () )
          .def("propagate", &I3PropagatorServicePROPOSAL::Propagate)
	;
}
