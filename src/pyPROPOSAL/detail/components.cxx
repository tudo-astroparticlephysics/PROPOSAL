
#include "PROPOSAL/medium/Components.h"
#include "pyPROPOSAL/pyBindings.h"

#define COMPONENT_DEF(module, cls)                             \
    py::class_<Components::cls, Component,         \
               std::shared_ptr<Components::cls>>(module, #cls) \
        .def(py::init<double>(), py::arg("atom_in_molecule") = 1.0);

namespace py = pybind11;
using namespace PROPOSAL;

void init_components(py::module& m) {
    py::module m_sub = m.def_submodule("component");

    m_sub.def("component_map", []() {return *Component::component_map;});

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

    py::class_<Component, std::shared_ptr<Component>>(
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
        .def("__str__", &py_print<Component>)
        .def_property_readonly("name", &Component::GetName,
                               R"pbdoc(
                    Get name of component.

                    Returns:
                        str: Name of component
                )pbdoc")
        .def_property_readonly("nuclear_charge",
                               &Component::GetNucCharge,
                               R"pbdoc(
                    Get nuclear charge of component.

                    Returns:
                        float: Nuclear charge of component
                )pbdoc")
        .def_property_readonly("atomic_number",
                               &Component::GetAtomicNum,
                               R"pbdoc(
                    Get atomic number of component.

                    Returns:
                        float: Atomic number of component
                )pbdoc")
        .def_property_readonly("atoms_in_molecule",
                               &Component::GetAtomInMolecule,
                               R"pbdoc(
                    Get number of atoms in one molecule.

                    Returns:
                        float: number of atoms in molecule
                )pbdoc")
        .def_property_readonly("log_constant",
                               &Component::GetLogConstant,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("bprime", &Component::GetBPrime,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("average_nucleon_weight",
                               &Component::GetAverageNucleonWeight,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("wood_saxon", &Component::GetWoodSaxon,
                               R"pbdoc(
                    Explanation still has to be added.
                )pbdoc")
        .def_property_readonly("hash", &Component::GetHash, R"pbdoc(
                    Hash value of component.
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
