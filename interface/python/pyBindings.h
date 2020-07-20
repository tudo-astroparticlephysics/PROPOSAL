#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


template <class T>
std::string py_print(const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;
