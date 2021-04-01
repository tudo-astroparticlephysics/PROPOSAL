
#include "PROPOSAL/math/InterpolantBuilder.h"
#include "PROPOSAL/math/Function.h"
#include "PROPOSAL/math/Spline.h"
#include "pyPROPOSAL/pyBindings.h"

namespace py = pybind11;
using namespace PROPOSAL;

void init_math(py::module& m) {
    py::module m_sub = m.def_submodule("math");

    py::class_<Interpolant1DBuilder::Definition, std::shared_ptr<Interpolant1DBuilder::Definition>>(m_sub, "Interpolant1DBuilderDefinition")
        .def(py::init<>())
        .def_readwrite("function1d", &Interpolant1DBuilder::Definition::function1d)
        .def_readwrite("nodes", &Interpolant1DBuilder::Definition::nodes)
        .def_readwrite("xmin", &Interpolant1DBuilder::Definition::xmin)
        .def_readwrite("xmax", &Interpolant1DBuilder::Definition::xmax)
        .def_readwrite("romberg", &Interpolant1DBuilder::Definition::romberg)
        .def_readwrite("rational", &Interpolant1DBuilder::Definition::rational)
        .def_readwrite("relative", &Interpolant1DBuilder::Definition::relative)
        .def_readwrite("isLog", &Interpolant1DBuilder::Definition::isLog)
        .def_readwrite("rombergY", &Interpolant1DBuilder::Definition::rombergY)
        .def_readwrite("rationalY", &Interpolant1DBuilder::Definition::rationalY)
        .def_readwrite("relativeY", &Interpolant1DBuilder::Definition::relativeY)
        .def_readwrite("logSubst", &Interpolant1DBuilder::Definition::logSubst) ;

    py::class_<Interpolant2DBuilder::Definition, Interpolant1DBuilder::Definition, std::shared_ptr<Interpolant2DBuilder::Definition>>(m_sub, "Interpolant2DBuilderDefinition")
        .def(py::init<>())
        .def_readwrite("function2d", &Interpolant2DBuilder::Definition::function2d)
        .def_readwrite("nodes2", &Interpolant2DBuilder::Definition::nodes2)
        .def_readwrite("x2min", &Interpolant2DBuilder::Definition::x2min)
        .def_readwrite("x2max", &Interpolant2DBuilder::Definition::x2max)
        .def_readwrite("romberg2", &Interpolant2DBuilder::Definition::romberg2)
        .def_readwrite("rational2", &Interpolant2DBuilder::Definition::rational2)
        .def_readwrite("relative2", &Interpolant2DBuilder::Definition::relative2)
        .def_readwrite("isLog2", &Interpolant2DBuilder::Definition::isLog2);

    py::class_<Polynom, std::shared_ptr<Polynom>>(m_sub, "Polynom")
        .def("__call__", py::vectorize(&Polynom::evaluate), py::arg("argument"))
        .def("__str__", &py_print<Polynom>)
        .def(py::init<std::vector<double>>())
        .def("derivate", &Polynom::GetDerivative)
        .def("antiderivative", &Polynom::GetAntiderivative, py::arg("constant"))
        .def_property_readonly("coeff", &Polynom::GetCoefficient);

    py::class_<Spline, std::shared_ptr<Spline>>(m_sub, "Spline", R"pbdoc(
            Splines are widely used in PROPOSAL to interpolate between nodes.
            This allows data to be compressed by storing only the sampling
            points and efficiently evaluating the functions using splines.
        )pbdoc")
        .def("__call__", py::vectorize(&Spline::evaluate), py::arg("argument"),
             R"pbdoc(
             Evaluating Spline at given argument.

             Parameters:
                argument (double): Argument of the function.

             Return:
                double: Value of the function.
             )pbdoc")
        .def("__str__", &py_print<Spline>)
        .def("derivate", &Spline::Derivative, R"pbdoc(
            Calculate the derivation of each subintervall and save it inplace.
        )pbdoc")
        .def("antiderivative", &Spline::Antiderivative, py::arg("constant"),
             R"pbdoc(
            Calculate the antiderivation of each subintervall and save it inplace.

            Parameters:
                constant (double): Integration constant for each subintervall.

        )pbdoc")
        .def("save", &Spline::save, py::arg("path"), py::arg("binary"), R"pbdoc(
            Save current spline in file. Calculation of splineintervalls can be
            omitted by cost of reading the file.  Currently only textfile is
            supported.

            Parameters:
                path (string): name of saved file
                binary (bool): use of binary type

            Return:
                bool: filediscriptor if save was successfully
            )pbdoc");

    py::class_<Linear_Spline, std::shared_ptr<Linear_Spline>, Spline>(
        m_sub, "Linear_spline")
        .def(py::init<std::vector<double>, std::vector<double>>(), py::arg("x"),
             py::arg("y"))
        .def(py::init<std::vector<Polynom>, std::vector<double>>(),
             py::arg("Polynom"), py::arg("definition_area"));

    py::class_<Cubic_Spline, std::shared_ptr<Cubic_Spline>, Spline>(
        m_sub, "Cubic_spline")
        .def(py::init<std::vector<double>, std::vector<double>>(), py::arg("x"),
             py::arg("y"))
        .def(py::init<std::vector<Polynom>, std::vector<double>>(),
             py::arg("Polynom"), py::arg("definition_area"))
        .def(py::init<std::string, bool>(), py::arg("Filename"),
             py::arg("binary"));
}
