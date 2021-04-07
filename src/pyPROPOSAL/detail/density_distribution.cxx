#include "PROPOSAL/density_distr/density_exponential.h"
#include "PROPOSAL/density_distr/density_homogeneous.h"
#include "PROPOSAL/density_distr/density_polynomial.h"
#include "PROPOSAL/density_distr/density_splines.h"
#include "pyPROPOSAL/pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_density_distribution(py::module& m) {
    py::module m_sub = m.def_submodule("density_distribution");

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
            .def(py::init<double>(), py::arg("mass_density"));

    py::class_<Density_exponential, Density_distr,
            std::shared_ptr<Density_exponential>>(m_sub,
                                                  "density_exponential")
            .def(py::init<const Axis&, double, double, double>(), py::arg("density_axis"),
                 py::arg("sigma"), py::arg("d0"), py::arg("mass_density"));

    py::class_<Density_polynomial, Density_distr,
            std::shared_ptr<Density_polynomial>>(m_sub, "density_polynomial")
            .def(py::init<const Axis&, const Polynom&, double>(), py::arg("density_axis"),
                 py::arg("polynom"), py::arg("mass_density"));

    py::class_<Density_splines, Density_distr,
            std::shared_ptr<Density_splines>>(m_sub, "density_splines")
            .def(py::init<const Axis&, const Spline&, double>(), py::arg("density_axis"),
                 py::arg("splines"), py::arg("mass_density"));

    py::class_<Axis, std::shared_ptr<Axis>>(m_sub, "Density_axis")
            .def_property_readonly("reference_point", &Axis::GetFp0)
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

    py::class_<CartesianAxis, Axis, std::shared_ptr<CartesianAxis>>(m_sub, "cartesian_axis")
        .def(py::init<const Vector3D&, const Vector3D&>(), py::arg("axis"),
             py::arg("reference_point"));

    py::class_<RadialAxis, Axis, std::shared_ptr<RadialAxis>>(m_sub, "radial_axis")
            .def(py::init<const Vector3D&>(), py::arg("reference_point"));
}
