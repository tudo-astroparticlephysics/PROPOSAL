#include "PROPOSAL/PROPOSAL.h"
#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_scattering(py::module& m) {
    py::module m_sub = m.def_submodule("scattering");

    py::class_<Scattering, std::shared_ptr<Scattering>>(m_sub, "Scattering")
        .def("scatter", &Scattering::Scatter)
        .def_property_readonly("particle", &Scattering::GetParticle);

    py::class_<ScatteringMoliere, std::shared_ptr<ScatteringMoliere>,
               Scattering>(m_sub, "Moliere")
        .def(py::init<Particle&, const Medium&>());

    py::class_<ScatteringHighlandIntegral,
               std::shared_ptr<ScatteringHighlandIntegral>, Scattering>(
        m_sub, "HighlandIntegral")
        .def(py::init<Particle&, Utility&, InterpolationDef>());

    py::class_<ScatteringHighland, std::shared_ptr<ScatteringHighland>,
               Scattering>(m_sub, "Highland")
        .def(py::init<Particle&, const Medium&>());

    py::class_<ScatteringNoScattering, std::shared_ptr<ScatteringNoScattering>,
               Scattering>(m_sub, "NoScattering")
        .def(py::init<Particle&, const Medium&>());

    py::enum_<ScatteringFactory::Enum>(m_sub, "ScatteringModel")
        .value("HighlandIntegral", ScatteringFactory::HighlandIntegral)
        .value("Moliere", ScatteringFactory::Moliere)
        .value("Highland", ScatteringFactory::Highland)
        .value("NoScattering", ScatteringFactory::NoScattering);
}
