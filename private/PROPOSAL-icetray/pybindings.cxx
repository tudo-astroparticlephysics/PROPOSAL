
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL-icetray/SimplePropagator.h"
#include "icetray/load_project.h"

I3_PYTHON_MODULE(PROPOSAL)
{
    load_project("PROPOSAL", false);

    using namespace PROPOSAL;
    using namespace boost::python;
    import("icecube.dataclasses");
    import("icecube.sim_services");

    class_<I3PropagatorServicePROPOSAL,
           boost::shared_ptr<I3PropagatorServicePROPOSAL>,
           bases<I3PropagatorService>,
           boost::noncopyable>(
        "I3PropagatorServicePROPOSAL",
        init<std::string>((arg("config_file") = I3PropagatorServicePROPOSAL::GetDefaultConfigFile()),
                          ":param config_file: Path to the config file\n"))
        .def("register_particletype", &I3PropagatorServicePROPOSAL::RegisterParticleType);

    class_<SimplePropagator, boost::shared_ptr<SimplePropagator>, boost::noncopyable>(
        "SimplePropagator",
        init<I3Particle::ParticleType, std::string, double, double, double>((arg("type")   = I3Particle::MuMinus,
                                                                             arg("medium") = "ice",
                                                                             arg("ecut")   = -1.0,
                                                                             arg("vcut")   = -1.0,
                                                                             arg("rho")    = -1.0)))
        .def("set_seed", &SimplePropagator::SetSeed)
        .def("propagate",
             &SimplePropagator::propagate,
             (args("p"), arg("distance"), arg("secondaries") = boost::shared_ptr<std::vector<I3Particle> >()));
}
