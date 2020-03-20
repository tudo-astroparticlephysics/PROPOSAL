
#include "PROPOSAL/medium/Components.h"
#include "pyBindings.h"

#define COMPONENT_DEF(module, cls)                             \
    py::class_<Components::cls, Components::Component,         \
               std::shared_ptr<Components::cls>>(module, #cls) \
        .def(py::init<double>(), py::arg("atom_in_molecule") = 1.0);

namespace py = pybind11;
using namespace PROPOSAL;

void init_components(py::module& m) {
    py::module m_sub = m.def_submodule("component");

    m_sub.doc() = R"pbdoc(
        You could create a new component or select one of the implemented.
        Components are used to define a medium as part of a sector.

        A existing component can be called for example with:

        >>> hydro = proposal.component.Hydrogen()
        >>> hydro.atomic_number
        1.00794

        There listed components can be created without initialization:

        * Hydrogen
        * Carbon
        * Nitrogen
        * Oxygen
        * Sodium
        * Magnesium
        * Sulfur
        * Argon
        * Potassium
        * Calcium
        * Iron
        * Copper
        * Lead
        * Uranium
        * StandardRock
        * FrejusRock

        Otherwise you have to initalize the component yourself.
    )pbdoc";

    py::class_<Components::Component, std::shared_ptr<Components::Component>>(
        m_sub, "Component")
        .def(py::init<std::string, double, double, double>(), py::arg("name"),
             py::arg("charge"), py::arg("atomic_num"),
             py::arg("atom_in_molecule"), R"pbdoc(
            Creating a new static component.

            Args:
               name (str):  The name of component.
               charge (float):  Charge in units of Coulomb.
               atomic_num (float): Atom number in periodic table.
               atomic_in_molecule (float): Number of atoms in molecule.
        )pbdoc")
        .def("__str__", &py_print<Components::Component>)
        .def_property_readonly("name", &Components::Component::GetName,
                               R"pbdoc(
                    Get name of component.

                    Returns:
                        str: Name of component
                )pbdoc")
        .def_property_readonly("nuclear_charge",
                               &Components::Component::GetNucCharge,
                               R"pbdoc(
                    Get nuclear charge of component.

                    Returns:
                        float: Nuclear charge of component
                )pbdoc")
        .def_property_readonly("atomic_number",
                               &Components::Component::GetAtomicNum,
                               R"pbdoc(
                    Get atomic number of component.

                    Returns:
                        float: Atomic number of component
                )pbdoc")
        .def_property_readonly("atoms_in_molecule",
                               &Components::Component::GetAtomInMolecule,
                               R"pbdoc(
                    Get number of atoms in one molecule.

                    Returns:
                        float: number of atoms in molecule
                )pbdoc")
        .def_property_readonly("log_constant",
                               &Components::Component::GetLogConstant,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("bprime", &Components::Component::GetBPrime,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("average_nucleon_weight",
                               &Components::Component::GetAverageNucleonWeight,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("wood_saxon", &Components::Component::GetWoodSaxon,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc");

    COMPONENT_DEF(m_sub, Hydrogen)
    COMPONENT_DEF(m_sub, Carbon)
    COMPONENT_DEF(m_sub, Nitrogen)
    COMPONENT_DEF(m_sub, Oxygen)
    COMPONENT_DEF(m_sub, Sodium)
    COMPONENT_DEF(m_sub, Magnesium)
    COMPONENT_DEF(m_sub, Sulfur)
    COMPONENT_DEF(m_sub, Argon)
    COMPONENT_DEF(m_sub, Potassium)
    COMPONENT_DEF(m_sub, Calcium)
    COMPONENT_DEF(m_sub, Iron)
    COMPONENT_DEF(m_sub, Copper)
    COMPONENT_DEF(m_sub, Lead)
    COMPONENT_DEF(m_sub, Uranium)
    COMPONENT_DEF(m_sub, StandardRock)
    COMPONENT_DEF(m_sub, FrejusRock)
}
