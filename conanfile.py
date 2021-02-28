from conans import ConanFile, CMake, tools
import os


class PROPOSALConan(ConanFile):
    name = "PROPOSAL"
    homepage = "https://github.com/tudo-astroparticlephysics/PROPOSAL"
    license = "GNU Lesser General Public License v3.0"
    description = "the very best lepton and photon propagator"
    topics = ("propagator", "lepton", "photon", "stochastic")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "build_testing": [True, False],
        "build_python": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "build_testing": False,
        "build_python": False,
    }
    generators = "cmake"
    _cmake = None

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def requirements(self):
        self.requires("cubic-interpolation/0.1.1")
        self.requires("spdlog/1.8.2")
        if self.options.build_python:
            self.requires("pybind11/2.6.2")
        if self.options.build_testing:
            self.requires("gtest/1.10.0")

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake
        self._cmake = CMake(self)
        self._cmake.definitions["BUILD_TESTING"] = self.options.build_testing
        self._cmake.definitions["BUILD_PYTHON"] = self.options.build_python
        self._cmake.configure()
        return self._cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy("LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)
