
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "pyPROPOSAL/pyBindings.h"

#define MEDIUM_DEF(module, cls)                                 \
    py::class_<cls, Medium, std::shared_ptr<cls>>(module, #cls) \
        .def(py::init<>());


namespace py = pybind11;
using namespace PROPOSAL;

void init_medium(py::module& m) {
    py::module m_sub = m.def_submodule("medium");

    m_sub.doc() = R"pbdoc(
            A medium is a object which contains serveral physical constants
            about the material. Media are constant objects and can not be
            modified once there are initalized, with the exception of their
            density distribution. There are several preimplemented media.

            +------------+------+----------+------------------+--------------------+
            | Water      | Ice  | Salt     | CalciumCarbonate | StandardRock       |
            +------------+------+----------+------------------+--------------------+
            | FrejusRock | Iron | Hydrogen | Lead             | Copper             |
            +------------+------+----------+------------------+--------------------+
            | Uranium    | Air  | Paraffin | AntaresWater     | CascadiaBasinWater |
            +------------+------+----------+------------------+--------------------+

            By default a media will be assumed as homogen distributed. There
            is the possibility to create inhomogen densities as listed below.
            Split densities are faster than multiple sectors and are therefore
            recommended when using different density profiles.

            +---------------------+--------------------+-----------------+
            | Density_exponential | Density_polynomial | Density_splines |
            +---------------------+--------------------+-----------------+

            There is currently an issue in the calculation of the LPM effect.
            The calculation of the LPM effect depends on the
            bremsstrahlungcrosssection which again depends on the depth.
            When generating the interpolation table, the depth is assumed to
            be constant.
            )pbdoc";

    py::class_<Medium, std::shared_ptr<Medium>>(m_sub, "Medium", R"pbdoc(
			Each Medium is a container of physical constants. There are several
			preimplemented media. They can be scaled by a correction factor if
			a homogen density is used. Otherwise the density_correction take
			care of scaling tasks.

            Example:
				Creating a homogeneous air medium with a linear correction
				factor.

				>>> correction_factor = 0.9
				>>> air = proposal.medium.Air(correction_factor)
				)pbdoc")
        .def("__str__", &py_print<Medium>)
        // .def(py::init<>())
        .def(py::init<
                 std::string,  double, double, double, double, double,
                 double, double, double,
                 const std::vector<Component>&>(),
             py::arg("name"),  py::arg("I"), py::arg("C"),
             py::arg("a"), py::arg("m"), py::arg("X0"), py::arg("X1"),
             py::arg("d0"), py::arg("massDensity"), py::arg("components"))
        .def_property_readonly("sum_charge", &Medium::GetSumCharge,
                               R"pbdoc(Sum of charges of all nuclei.)pbdoc")
        .def_property_readonly("ratio_ZA", &Medium::GetZA,
                               R"pbdoc(Ratio of atomic and mass number.)pbdoc")
        .def_property_readonly("ionization_potential", &Medium::GetI,
                               R"pbdoc(Ionization potential in [eV])pbdoc")
        .def_property_readonly(
            "radiation_length",
            (double (Medium::*)() const) & Medium::GetRadiationLength,
            R"pbdoc(Radiation length of the medium in [].)pbdoc")
        .def_property_readonly(
            "mass_density", &Medium::GetMassDensity,
            R"pbdoc(Mass density of the medium in [g/cm3].)pbdoc")
        .def_property_readonly(
            "mol_density", &Medium::GetMolDensity,
            R"pbdoc(Mol density of the medium in [#/cm3].)pbdoc")
        .def_property_readonly(
            "average_nucleon_weigth", &Medium::GetMM,
            R"pbdoc(Averaged component nucleon weights in [MeV])pbdoc")
        .def_property_readonly("sum_nucleons", &Medium::GetSumNucleons,
                               R"pbdoc(Sum of nucleons.)pbdoc")
        .def_property_readonly(
            "num_components", &Medium::GetNumComponents,
            R"pbdoc(Number of components preserved in the medium.)pbdoc")
        .def_property_readonly(
            "components", &Medium::GetComponents,
            R"pbdoc(List of components preserved in the medium.)pbdoc")
        .def_property_readonly(
            "name", &Medium::GetName,
            R"pbdoc(Internal name of the particle. Should be unique to prevent
				undefined behaviour.)pbdoc")
        .def_property_readonly("hash", &Medium::GetHash,
            R"pbdoc(Get medium hash.)pbdoc");

    MEDIUM_DEF(m_sub, Water)
    MEDIUM_DEF(m_sub, Ice)
    MEDIUM_DEF(m_sub, Salt)
    MEDIUM_DEF(m_sub, CalciumCarbonate)
    MEDIUM_DEF(m_sub, StandardRock)
    MEDIUM_DEF(m_sub, FrejusRock)
    MEDIUM_DEF(m_sub, Iron)
    MEDIUM_DEF(m_sub, Hydrogen)
    MEDIUM_DEF(m_sub, Lead)
    MEDIUM_DEF(m_sub, Copper)
    MEDIUM_DEF(m_sub, Uranium)
    MEDIUM_DEF(m_sub, Air)
    MEDIUM_DEF(m_sub, Paraffin)
    MEDIUM_DEF(m_sub, AntaresWater)
    MEDIUM_DEF(m_sub, CascadiaBasinWater)
}
