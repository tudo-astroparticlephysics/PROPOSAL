
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/scattering/ScatteringFactory.h"
#include "PROPOSAL/scattering/ScatteringHighland.h"
#include "PROPOSAL/scattering/ScatteringHighlandIntegral.h"
#include "PROPOSAL/scattering/ScatteringMoliere.h"
#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_scattering(py::module& m) {
    py::module m_sub = m.def_submodule("scattering");

    py::class_<Scattering, std::shared_ptr<Scattering>>(m_sub, "Scattering")
        .def("scatter", &Scattering::Scatter,
                py::arg("grammage"), py::arg("e_i"), py::arg("e_f"), py::arg("position"),  py::arg("random_numbers"),
            R"pbdoc(
                Calculate a random averaged scatterangle `u` alonge the given
                displacement`dr` and the particle directions after distance
                `n_i`.

                Args:
                    dr(double): displacement of particle
                    ei(double): inital energy
                    ef(double): final energy
            )pbdoc");
}
