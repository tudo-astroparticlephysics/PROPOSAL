
#include "PROPOSAL/medium/Components.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/density_distr/density_exponential.h"
#include "PROPOSAL/medium/density_distr/density_homogeneous.h"
#include "PROPOSAL/medium/density_distr/density_polynomial.h"
#include "PROPOSAL/medium/density_distr/density_splines.h"
#include "pyBindings.h"

#define MEDIUM_DEF(module, cls)                                 \
    py::class_<cls, Medium, std::shared_ptr<cls>>(module, #cls) \
        .def(py::init<double>(), py::arg("density_correction") = 1.0);

#define AXIS_DEF(module, cls)                                 \
    py::class_<cls, Axis, std::shared_ptr<cls>>(module, #cls) \
        .def(py::init<Vector3D, Vector3D>(), py::arg("axis"), \
             py::arg("reference_point"));

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
        .def(py::init<>())
        .def(py::init<
                 std::string, double, double, double, double, double, double,
                 double, double, double,
                 const std::vector<Components::Component>&>(),
             py::arg("name"), py::arg("rho"), py::arg("I"), py::arg("C"),
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
        .def_property("density_distribution", &Medium::GetDensityDistribution,
                      &Medium::SetDensityDistribution,
                      R"pbdoc(Density distribution of medium.)pbdoc");

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

    py::class_<Density_distr, std::shared_ptr<Density_distr>>(
        m_sub, "density_distribution", R"pbdoc(
		Density dirstribution is a factor for the medium dependend of the place
		of the particle.  There has to be an axis :math:`\vec{a}` , a reference
		point :math:`\vec{r}_p` and a function :math:`\rho` to calculate the
		density correction.

		Note:
			Partical locations :math:`\vec{x}` will be written in dependent of
			propagated distance :math:`l` along particle direction
			:math:`\vec{d}` since last interaction point :math:`\vec{x}_r` and.

			.. math::

				\vec{x} = \vec{x}_i + l \cdot \vec{d}

		Density function will be evaluated at the projection of the particle
 		position onto the density axis. If corrected displacement is in the
		dector, it will be returned, otherwise a failure will be raised and
		catched by the sector class to propagate the particle continously to
		the sector border.

            )pbdoc")
        .def("correct", &Density_distr::Correct, py::arg("xi"),
             py::arg("direction"), py::arg("displacement"),
             py::arg("distance_to_border"),
             R"pbdoc(
			Calculate displacemet :math:`l` equivalent by inverting the equation
			below. The left hand side represent the displacement in reference
			medium which will be translated with density correction.

			.. math::
				\int_{E_i}^{E_r} \frac{1}{f(E)} dE= \int^{\vec{x}_i + l}_{\vec{x}_i} \rho(\vec{a} \vec{x}) dx

			If inverting the integral is not possible the newton method will
			be used. When the displacement equivalent is outside the sector an
			error will be raised, because it imply undefined behaviour.

			Parameters:
				xi (Vector3D): particle position
				direction (Vector3D): particle direction
				displacement (float): distance since last stochastical loss
				distance_to_border (float): distance to sector border

			Return:
				float: displacement equivalent
            )pbdoc")
        .def("evaluate", &Density_distr::Evaluate, py::arg("xi"),
             R"pbdoc(
			Evaluate the density correction at given point.

			.. math::
				\rho(\vec{a} \vec{x})

			Parameters:
				xi (Vector3D): particle position

			Return:
				float: density

            )pbdoc")
        .def("calculate", &Density_distr::Calculate, py::arg("xi"),
             py::arg("direction"), py::arg("displacement"),
             R"pbdoc(
			The depth will be calculated. Therefore, the dense integral is
			solved from the last interaction point to the next interaction
			point.

			.. math::
				\int^{\vec{x}_i + l \cdot \vec{d}}_{\vec{x}_i} \rho(\vec{a} \vec{x}) dx

			Is equivalent to the correction factor to consider density
			distribution for this place.

			Parameters:
				xi (Vector3D): particle position
				direction (Vector3D): particle direction
				displacement (float): distance since last stochastical loss

			Return:
				float: density correction
            )pbdoc");

    py::class_<Density_homogeneous, Density_distr,
               std::shared_ptr<Density_homogeneous>>(m_sub,
                                                     "density_homogeneous")
        .def(py::init<>());

    py::class_<Density_exponential, Density_distr,
               std::shared_ptr<Density_exponential>>(m_sub,
                                                     "density_exponential")
        .def(py::init<const Axis&, double>(), py::arg("density_axis"),
             py::arg("sigma"));

    py::class_<Density_polynomial, Density_distr,
               std::shared_ptr<Density_polynomial>>(m_sub, "density_polynomial")
        .def(py::init<const Axis&, const Polynom&>(), py::arg("density_axis"),
             py::arg("polynom"));

    py::class_<Density_splines, Density_distr,
               std::shared_ptr<Density_splines>>(m_sub, "density_splines")
        .def(py::init<const Axis&, const Spline&>(), py::arg("density_axis"),
             py::arg("splines"));

    py::class_<Axis, std::shared_ptr<Axis>>(m_sub, "Density_axis")
        .def_property_readonly("fAxis", &Axis::GetAxis)
        .def_property_readonly("refernce_point", &Axis::GetFp0)
        .def("depth", &Axis::GetDepth, py::arg("position"),
             R"pbdoc(
                Calculates in dependence of the particle position the depthcorrection.

                Parameters:
                    position (Vector3D): particle position

                Return:
                    float: depthcorrection
            )pbdoc")
        .def("depth", &Axis::GetEffectiveDistance, py::arg("position"),
             py::arg("direction"),
             R"pbdoc(
                Calculates in dependence of the particle position the effective Distance.

                Parameters:
                    position (Vector3D): particle position
                    direction (Vector3D): direction position

                Return:
                    float: effective Distance
            )pbdoc");

    AXIS_DEF(m_sub, RadialAxis);
    AXIS_DEF(m_sub, CartesianAxis);
}
