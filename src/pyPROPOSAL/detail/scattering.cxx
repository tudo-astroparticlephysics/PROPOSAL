
#include "PROPOSAL/scattering/Scattering.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/scattering/ScatteringMultiplier.h"
#include "PROPOSAL/scattering/multiple_scattering/Highland.h"
#include "PROPOSAL/scattering/multiple_scattering/HighlandIntegral.h"
#include "PROPOSAL/scattering/multiple_scattering/Moliere.h"
#include "PROPOSAL/scattering/multiple_scattering/ScatteringFactory.h"

#include "PROPOSAL/scattering/stochastic_deflection/bremsstrahlung/TsaiApproximationBremsstrahlung.h"

#include "pyPROPOSAL/pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_scattering(py::module& m)
{
    py::module m_sub = m.def_submodule("scattering");

    py::class_<multiple_scattering::Parametrization,
        std::shared_ptr<multiple_scattering::Parametrization>>(
        m_sub, "MultipleScattering")
        .def("scatter",
            &multiple_scattering::Parametrization::CalculateRandomAngle,
            py::arg("grammage"), py::arg("e_i"), py::arg("e_f"),
            py::arg("random_numbers"),
            R"pbdoc(
                Calculate a random averaged scatterangle `u` in cartesian coordinates.

                Args:
                    dr(double): displacement of particle
                    ei(double): inital energy
                    ef(double): final energy
            )pbdoc");

    py::class_<multiple_scattering::ScatteringOffset>(m_sub, "scattering_offset")
            .def_readwrite("sx", &multiple_scattering::ScatteringOffset::sx)
            .def_readwrite("sy", &multiple_scattering::ScatteringOffset::sy)
            .def_readwrite("tx", &multiple_scattering::ScatteringOffset::tx)
            .def_readwrite("ty", &multiple_scattering::ScatteringOffset::ty);

    m_sub.def("scatter_initial_direction", &multiple_scattering::ScatterInitialDirection);

    auto scattering_doc
        = R"pbdoc(Factory method for a MultipleScattering object.

                Args:
                    name(string): name of the scattering method
                    particle_def(ParticleDef): particle related constants
                    target(Medium): medium related constants

                Returns: MultipleScattering
    )pbdoc";

    m.def(
        "make_multiple_scattering",
        [](std::string const& n, ParticleDef const& p, Medium const& m,
            crosssection_list_t c, bool i) {
            return std::shared_ptr<multiple_scattering::Parametrization>(
                make_multiple_scattering(n, p, m, c, i));
        },
        py::arg("name"), py::arg("particle_def"), py::arg("target"),
        py::arg("cross"), py::arg("interpolate"), scattering_doc);

    m.def(
        "make_multiple_scattering",
        [](std::string const& n, ParticleDef const& p, Medium const& m) {
            return std::shared_ptr<multiple_scattering::Parametrization>(
                make_multiple_scattering(n, p, m));
        },
        py::arg("name"), py::arg("particle_def"), py::arg("target"),
        scattering_doc);

    py::class_<stochastic_deflection::Parametrization,
        std::shared_ptr<stochastic_deflection::Parametrization>>(
        m_sub, "StochasticDeflection")
        .def("n_rnd",
            &stochastic_deflection::Parametrization::RequiredRandomNumbers,
            R"pbdoc(Required random numbers for a stochastic deflection calculation)pbdoc")
        .def("type",
            &stochastic_deflection::Parametrization::GetInteractionType,
            R"pbdoc(Interaction type which causes the stochastic deflection calculation)pbdoc")
        .def("stochastic_deflection",
            &stochastic_deflection::Parametrization::
                CalculateStochasticDeflection,
            py::arg("initial_energy"), py::arg("loss_energy"),
            py::arg("random_numbers"),
            R"pbdoc(TODO: Doc is missing because it's not clear if second argument
            should be the lost energy or the final energy. Please contact the
            maintainers if required.)pbdoc");

    using deflect_ptr = std::shared_ptr<stochastic_deflection::Parametrization>;
    using deflect_list_t = std::vector<deflect_ptr>;

    m.def("make_stochastic_deflection",
            [](std::string const& n, ParticleDef const& p, Medium const& m) {
                return std::shared_ptr<stochastic_deflection::Parametrization>(
                        make_stochastic_deflection(n, p, m));
            },
            py::arg("name"), py::arg("particle_def"), py::arg("target"),
            scattering_doc);

    m.def(
        "make_default_stochastic_deflection",
        [](std::vector<InteractionType> t, ParticleDef const& p,
            Medium const& m) {
            auto v_shared = deflect_list_t();
            for (auto& v_i : make_default_stochastic_deflection(t, p, m))
                v_shared.emplace_back(v_i->clone());
            return v_shared;
        },
        py::arg("type"), py::arg("particle_def"), py::arg("target"));

    using multiple_scatter_t = multiple_scattering::Parametrization;

    py::class_<Scattering, std::shared_ptr<Scattering>>(m_sub, "Scattering")
        .def(py::init([](multiple_scatter_t const& s, deflect_list_t const& d) {
            return std::make_shared<Scattering>(s, d);
        }),
            py::arg("multiple_scatter"), py::arg("stochastic_deflection"))
        .def(py::init([](multiple_scatter_t const& s) {
            return std::make_shared<Scattering>(s, nullptr);
        }),
            py::arg("multiple_scatter"))
        .def(py::init([](deflect_list_t const& d) {
            return std::make_shared<Scattering>(nullptr, d);
        }),
            py::arg("stochastic_deflection"))
        .def("n_rnd_mulitple_scatter",
            &Scattering::MultipleScatteringRandomNumbers,
            R"pbdoc(Required number of random numbers to sample multiple scattering.
            Returns:
                int: required rnd numbers)pbdoc")
        .def("n_rnd_stochastic_deflect",
            &Scattering::StochasticDeflectionRandomNumbers,
            R"pbdoc(Required number of random numbers to sample a stochastic interaction for given type.
            Args:
                type (Interaction_Type): type of interaction
            Returns:
                int: required rnd numbers)pbdoc")
        .def("stochastic_deflection",
            &Scattering::CalculateStochasticDeflection<double, double,
                std::vector<double> const&>,
            py::arg("type"), py::arg("initial_energy"), py::arg("final_energy"),
            py::arg("rnd"),
            R"pbdoc(Sample stochastic defleciton angles in radians.
            Args:
                type (Interaction_Type): type of stochastic loss
                initial_energy (double): energy before continuous loss
                final_energy (double): energy after continuous loss
                rnd (list(double)): container of random numbers
            Returns:
                list(double): deflection angles)pbdoc")
        .def("multiple_scattering",
            &Scattering::CalculateMultipleScattering<double, double, double,
                const std::array<double, 4>&>,
            py::arg("grammage"), py::arg("inital_energy"),
            py::arg("final_energy"), py::arg("rnd"),
            R"pbdoc(Sample multiple scattering angles in cartesian coordinates.
            Args:
                initial_energy (double): energy before continuous loss
                final_energy (double): energy after continuous loss
                rnd (list(double)): container of random numbers
            Returns:
                list(double): averaged deflection over continuous loss and final direction)pbdoc");

    py::class_<ScatteringMultiplier, Scattering,
        std::shared_ptr<ScatteringMultiplier>>(m_sub, "ScatteringMultiplier")
        .def(py::init([](multiple_scatter_t const& s, deflect_list_t const& d,
                          double mm,
                          std::vector<std::pair<InteractionType, double>> dm) {
            return std::make_shared<ScatteringMultiplier>(s, d, mm, dm);
        }),
            py::arg("multiple_scatter"), py::arg("stochastic_deflection"),
            py::arg("multiple_scatter_multiplier"),
            py::arg("stochastic_deflection_multiplier"),
            R"pbdoc(Sample multiple scattering angles in cartesian coordinates.
            Args:
                multiple_scatter (MultipleScattering): multiple scattering calculator
                stochastic_deflection (StochasticDeflection): stochastic deflection calculator
                multiple_scatter_multiplier (double): multiple scatter factor
                stochastic_deflection_multiplier ([tuple(Interaction_Type, double)]): interaction types with corresponding multiplier )pbdoc")
        .def(py::init([](multiple_scatter_t const& s, double mm) {
            auto dm = std::vector<std::pair<InteractionType, double>>();
            return std::make_shared<ScatteringMultiplier>(s, nullptr, mm, dm);
        }),
            py::arg("multiple_scatter"), py::arg("multiple_scatter_multiplier"))
        .def(py::init([](deflect_list_t const& d,
                          std::vector<std::pair<InteractionType, double>> dm) {
            return std::make_shared<ScatteringMultiplier>(nullptr, d, 1., dm);
        }),
            py::arg("stochastic_deflection"),
            py::arg("stochastic_deflection_multiplier"));
}
