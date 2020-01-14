
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/Secondaries.h"
#include "PROPOSAL/decay/LeptonicDecayChannel.h"
#include "PROPOSAL/decay/ManyBodyPhaseSpace.h"
#include "PROPOSAL/decay/StableChannel.h"
#include "PROPOSAL/decay/DecayChannel.h"
#include "PROPOSAL/decay/TwoBodyPhaseSpace.h"
#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_decay(py::module& m) {
    py::module m_sub = m.def_submodule("decay");

    py::class_<DecayChannel, std::shared_ptr<DecayChannel>>(m_sub,
                                                            "DecayChannel")
        .def("__str__", &py_print<DecayChannel>)
        .def("__eq__", &DecayChannel::operator==)
        .def("__ne__", &DecayChannel::operator!=)
        .def("decay", &DecayChannel::Decay, "Decay the given particle")
        .def_static("boost", overload_cast_<DynamicData&, const Vector3D&, double, double>()(&DecayChannel::Boost))
        .def_static("boost", overload_cast_<Secondaries&, const Vector3D&, double, double>()(&DecayChannel::Boost));

    py::class_<LeptonicDecayChannelApprox,
               std::shared_ptr<LeptonicDecayChannelApprox>, DecayChannel>(
        m_sub, "LeptonicDecayChannelApprox")
        .def(py::init<const ParticleDef&, const ParticleDef&,
                      const ParticleDef&>());

    py::class_<LeptonicDecayChannel, std::shared_ptr<LeptonicDecayChannel>,
               DecayChannel>(m_sub, "LeptonicDecayChannel")
        .def(py::init<const ParticleDef&, const ParticleDef&,
                      const ParticleDef&>());

    py::class_<TwoBodyPhaseSpace, std::shared_ptr<TwoBodyPhaseSpace>,
               DecayChannel>(m_sub, "TwoBodyPhaseSpace")
        .def(py::init<ParticleDef, ParticleDef>());

    py::class_<ManyBodyPhaseSpace, std::shared_ptr<ManyBodyPhaseSpace>,
               DecayChannel>(m_sub, "ManyBodyPhaseSpace")
        .def(py::init<std::vector<const ParticleDef*>,
                      PROPOSAL::ManyBodyPhaseSpace::MatrixElementFunction>(),
             py::arg("particle_defs"), py::arg("matrix_element") = nullptr)
        .def_static("default_evaluate", &ManyBodyPhaseSpace::DefaultEvaluate,
                    "Return the default matrix element (default 1)")
        .def("evaluate", &ManyBodyPhaseSpace::Evaluate,
             "Return the matrix element (default 1)")
        .def("set_uniform_sampling", &DecayChannel::SetUniformSampling,
             "Decide to use uniform phase space sampling");

    py::class_<StableChannel, std::shared_ptr<StableChannel>, DecayChannel>(
        m_sub, "StableChannel")
        .def(py::init<>());

    py::class_<DecayTable, std::shared_ptr<DecayTable>>(m_sub, "DecayTable")
        .def(py::init<>())
        .def(py::init<const DecayTable&>())
        .def("__str__", &py_print<DecayTable>)
        .def("add_channel", &DecayTable::addChannel, "Add an decay channel")
        .def("select_channel", &DecayTable::SelectChannel, py::arg("rnd"),
             "Select an decay channel according to given branching ratios")
        .def("set_stable", &DecayTable::SetStable,
             "Define decay table for stable particles")
        .def("set_uniform_sampling", &DecayTable::SetUniformSampling,
             "Set whether to sample many body decays uniform in phase "
             "space");
}
