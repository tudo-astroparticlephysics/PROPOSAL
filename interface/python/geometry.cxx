
#include "PROPOSAL/geometry/GeometryFactory.h"
#include "PROPOSAL/geometry/Box.h"
#include "PROPOSAL/geometry/Cylinder.h"
#include "PROPOSAL/geometry/Sphere.h"
#include "pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_geometry(py::module& m) {
    py::module m_sub = m.def_submodule("geometry");

    m_sub.doc() = R"pbdoc(
        Every sector is defined by a specific medium and a geometry.
        There are three different classes defined to build a mathematical
        body. All of them are a object of type :meth:`Geometry`.
        Besides the information of the shape an object of the class
        geometry contains the position relativ to the coordinate origin.

        Based on the position of the geometry, distances of the propagated
        particle to geometry sizes can be determined.
    )pbdoc";

    py::enum_<Geometry_Type>(m_sub, "Shape")
        .value("Sphere", Geometry_Type::SPHERE)
        .value("Box", Geometry_Type::BOX)
        .value("Cylinder", Geometry_Type::CYLINDER);

    py::class_<Geometry, std::shared_ptr<Geometry>>(m_sub, "Geometry")
        .def("__str__", &py_print<Geometry>)
        .def("is_infront", &Geometry::IsInfront,
             R"pbdoc(
            Check if particle is in fron of the geometry.

            Parameters:
                arg0 (Vector3D): particle position
                arg1 (Vector3D): particle direction

            Return:
                bool: Is particle in front of the geometry?
        )pbdoc")
        .def("is_inside", &Geometry::IsInside,
             R"pbdoc(
            Check if particle is in the geometry.

            Parameters:
                arg0 (Vector3D): particle position
                arg1 (Vector3D): particle direction

            Return:
                bool: Is particle in the geometry?
        )pbdoc")
        .def("is_behind", &Geometry::IsBehind,
             R"pbdoc(
            Check if particle has passed the geometry.

            Parameters:
                arg0 (Vector3D): particle position
                arg1 (Vector3D): particle direction

            Return:
                bool: Is particle in front of the geometry?
        )pbdoc")
        .def("distance_to_border", &Geometry::DistanceToBorder,
             py::arg("position"), py::arg("direction"),
             R"pbdoc(
            Calculates in dependence of the particle position and direction
            the distance to the next border

            Parameters:
                position (Vector3D): particle position
                direction (Vector3D): particle direction

            Return:
                float: distance to border
        )pbdoc")
        .def("distance_to_closet_approach",
             &Geometry::DistanceToClosestApproach,
             R"pbdoc(
            Calculates in dependence of the particle position and direction
            the distance where the particle pass the geometry center with
            minimal distance.

            Parameters:
                arg0 (Vector3D): particle position
                arg1 (Vector3D): particle direction

            Return:
                float: distance to closest approach
        )pbdoc")
        .def_property_readonly("name", &Geometry::GetName,
                               R"pbdoc(
            name of the geometry
        )pbdoc")
        .def_property("position", &Geometry::GetPosition,
                      &Geometry::SetPosition,
                      R"pbdoc(
            position of the geometry center.
        )pbdoc")
        .def_property("hierarchy", &Geometry::GetHierarchy,
                      &Geometry::SetHierarchy,
                      R"pbdoc(
            hierachy of the geometry. If sectors overlap, the sector
            with the highest hierachy will be selected.
        )pbdoc");

    py::class_<Sphere, std::shared_ptr<Sphere>, Geometry>(m_sub, "Sphere")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double>())
        .def(py::init<const Sphere&>())
        .def_property("inner_radius", &Sphere::GetInnerRadius,
                      &Sphere::SetInnerRadius,
                      R"pbdoc(
                inner radius of the sphere
            )pbdoc")
        .def_property("radius", &Sphere::GetRadius, &Sphere::SetRadius,
                      R"pbdoc(
                outer radius of the sphere
            )pbdoc");

    py::class_<Box, std::shared_ptr<Box>, Geometry>(m_sub, "Box")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Box&>())
        .def_property("width", &Box::GetX, &Box::SetX,
                      R"pbdoc(
                width of the box (x-axis)
            )pbdoc")
        .def_property("height", &Box::GetY, &Box::SetY,
                      R"pbdoc(
                height of the box (y-axis)
            )pbdoc")
        .def_property("depth", &Box::GetZ, &Box::SetZ,
                      R"pbdoc(
                depth of the box (z-axis)
            )pbdoc");

    py::class_<Cylinder, std::shared_ptr<Cylinder>, Geometry>(m_sub, "Cylinder",
                                                              R"pbdoc(
                A cylinder can be created as a hollow cylinder.
                For this purpose, a corresponding radius must be
                selected for the bore along the main axis. A
                cylinder without a bore is equal to a bore radius
                equal to zero.
            )pbdoc")
        .def(py::init<>())
        .def(py::init<Vector3D, double, double, double>())
        .def(py::init<const Cylinder&>())
        .def_property("inner_radius", &Cylinder::GetInnerRadius,
                      &Cylinder::SetInnerRadius,
                      R"pbdoc(
                the inner radius of the bore through the main axis
                of the cylinder
            )pbdoc")
        .def_property("radius", &Cylinder::GetRadius, &Cylinder::SetRadius,
                      R"pbdoc(
                radius of outer shell of the cylinder
            )pbdoc")
        .def_property("height", &Cylinder::GetZ, &Cylinder::SetZ,
                      R"pbdoc(
                height of the cylinder
            )pbdoc");
}
