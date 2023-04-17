import os
from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy, get, rmdir
from conan.tools.scm import Version

required_conan_version = ">=1.53.0"


class PROPOSALConan(ConanFile):
    name = "proposal"
    homepage = "https://github.com/tudo-astroparticlephysics/PROPOSAL"
    license = "LGPL-3.0"
    package_type = "library"
    description = "monte Carlo based lepton and photon propagator"
    topics = ("propagator", "lepton", "photon", "stochastic")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_python": [True, False],
        "with_testing": [True, False],
        "with_documentation": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_python": False,
        "with_testing": False,
        "with_documentation": False,
    }

    def config_options(self):
        if self.settings.os == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self, src_folder="src")

    def requirements(self):
        # cubicinterpolation: headers are transitively included, and function calls are made
        # from implementation in headers (templates)
        self.requires("cubicinterpolation/0.1.5", transitive_headers=True, transitive_libs=True)
        # spdlog: requires transitive_libs due to direct calls to functionality from headers
        self.requires("spdlog/1.11.0", transitive_headers=True, transitive_libs=True)
        # nlohmann_json: public headers include json.hpp and json_fwd.hpp
        self.requires("nlohmann_json/3.11.2", transitive_headers=True)
        if self.options.with_python:
            self.requires("pybind11/2.10.1")
        if self.options.with_testing:
            self.requires("boost/1.75.0")
            self.requires("gtest/1.11.0")
        if self.options.with_documentation:
            self.requires("doxygen/1.8.20")

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.install()
        if self.options.with_testing:
            cmake.test()

    def package(self):
        copy(self, "LICENSE.md", self.source_folder, os.path.join(self.package_folder, "licenses"))
        cmake = CMake(self)
        cmake.install()
        rmdir(self, os.path.join(self.package_folder, "lib", "cmake"))

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["BUILD_TESTING"] = self.options.with_testing
        tc.cache_variables["BUILD_PYTHON"] = self.options.with_python
        tc.cache_variables["BUILD_DOCUMENTATION"] = self.options.with_documentation
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "PROPOSAL")
        self.cpp_info.set_property("cmake_target_name", "PROPOSAL::PROPOSAL")
        self.cpp_info.libs = ["PROPOSAL"]

        # TODO: to remove in conan v2 once cmake_find_package_* generators removed
        self.cpp_info.names["cmake_find_package"] = "PROPOSAL"
        self.cpp_info.names["cmake_find_package_multi"] = "PROPOSAL"